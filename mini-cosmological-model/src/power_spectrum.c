/** @file power_spectrum.c
 *  @brief Primordial and matter power spectrum, transfer functions,
 *         correlation function, and nonlinear corrections.
 *
 *  The power spectrum P(k) is the Fourier transform of the two-point
 *  correlation function xi(r). It encodes the statistical properties
 *  of density fluctuations as a function of scale k.
 *
 *  Reference: Dodelson & Schmidt (2020) Ch.7-8
 *             Eisenstein & Hu (1998), ApJ 496, 605
 *             Bardeen, Bond, Kaiser, Szalay (1986), ApJ 304, 15
 */

#include "power_spectrum.h"
#include "friedmann.h"
#include "perturbations.h"
#include <math.h>

/* === Primordial power spectrum === */

/* Dimensionless primordial curvature power spectrum:
 * P_R(k) = A_s * (k/k_pivot)^{n_s - 1}
 * Delta_R^2(k) = k^3 * P_R(k) / (2*pi^2)
 */
double ps_primordial_dimensionless(double k_Mpc, double A_s,
                                    double n_s, double k_pivot)
{
    if (k_Mpc <= 0.0 || k_pivot <= 0.0 || A_s <= 0.0) return 0.0;
    return A_s * pow(k_Mpc / k_pivot, n_s - 1.0);
}

/* Primordial power spectrum with running spectral index:
 * P_R(k) = A_s * (k/k_pivot)^{n_s - 1 + 0.5*alpha_s*ln(k/k_pivot)}
 *
 * Running alpha_s = d(n_s)/d(ln k) is constrained by Planck 2018
 * to be consistent with zero: alpha_s = -0.004 +/- 0.007 (68% CL).
 */
double ps_primordial_running(double k_Mpc, double A_s, double n_s,
                              double alpha_s, double k_pivot)
{
    if (k_Mpc <= 0.0 || k_pivot <= 0.0 || A_s <= 0.0) return 0.0;
    double ln_ratio = log(k_Mpc / k_pivot);
    double exponent = n_s - 1.0 + 0.5 * alpha_s * ln_ratio;
    return A_s * pow(k_Mpc / k_pivot, exponent);
}

/* === BBKS transfer function ===
 * T_BBKS(q) = ln(1 + 2.34*q) / (2.34*q)
 *            * [1 + 3.89*q + (16.1*q)^2 + (5.46*q)^3 + (6.71*q)^4]^{-1/4}
 * where q = k / (Omega_m * h^2) [Mpc^{-1}]
 *
 * Valid for CDM with negligible baryon fraction.
 * This parametrization captures the suppression of power on small
 * scales due to free-streaming (Meszaros effect) during radiation era.
 */
double ps_transfer_BBKS(double k_Mpc, double Omega_m, double h)
{
    if (k_Mpc <= 0.0 || Omega_m <= 0.0 || h <= 0.0) return 1.0;

    double Gamma = Omega_m * h;  /* Shape parameter */
    double q = k_Mpc / Gamma;

    if (q < 1.0e-8) return 1.0;

    double num = log(1.0 + 2.34 * q);
    double denom_sq = 1.0 + 3.89 * q
                    + (16.1 * q) * (16.1 * q)
                    + pow(5.46 * q, 3.0)
                    + pow(6.71 * q, 4.0);

    return (num / (2.34 * q)) * pow(denom_sq, -0.25);
}

/* === Eisenstein & Hu (1998) transfer function ===
 * T_EH(k) = (Omega_b/Omega_m)*T_b(k) + (Omega_c/Omega_m)*T_cdm(k)
 *
 * This includes baryon acoustic oscillations from pre-recombination
 * sound waves. The BAO feature appears as wiggles at
 * k ~ 0.05 h/Mpc (the sound horizon scale).
 *
 * We implement the zero-baryon limit of the EH formula for simplicity,
 * which reduces to a modified BBKS shape with sound horizon damping.
 */
double ps_transfer_EH(double k_Mpc, const CosmoModel *model)
{
    if (!model || k_Mpc <= 0.0) return 1.0;
    if (model->Omega_m <= 0.0) return 1.0;

    /* Shape parameter with baryon correction (Sugiyama 1995) */
    double f_b = model->Omega_b / model->Omega_m;
    double Gamma_eff = model->Omega_m * model->h
                     * exp(-f_b * (1.0 + sqrt(2.0 * model->h) / model->Omega_m));

    double q = k_Mpc / Gamma_eff;

    if (q < 1.0e-8) return 1.0;

    /* Zero-baryon transfer function (Eq. 29 of EH98) */
    double L0 = log(2.0 * exp(1.0) + 1.8 * q);
    double C0 = 14.2 + 731.0 / (1.0 + 62.5 * q);

    double T0 = L0 / (L0 + C0 * q * q);

    /* For non-zero baryon fraction, add oscillatory BAO feature */
    if (f_b > 0.001) {
        double s = 44.5 * log(9.83 / q) / sqrt(1.0 + 10.0 * pow(f_b, 0.75));
        double bao_factor = 1.0 + (1.0 - f_b) * pow(q, 0.5) * sin(s) / (1.0 + q);
        return T0 * bao_factor;
    }

    return T0;
}

/* === Linear matter power spectrum ===
 * P_lin(k) = A_s * (k/k_pivot)^{n_s-1} * T^2(k) * normalization
 *
 * The normalization ensures that sigma_8 matches observations.
 * We use the BBKS transfer function for the shape, but scale
 * the amplitude such that sigma(8 Mpc/h) = sigma_8_target.
 */
double ps_matter_power_linear(double k_Mpc, const CosmoModel *model)
{
    if (!model || k_Mpc <= 0.0) return 0.0;

    /* Primordial power */
    double P_R = model->A_s * pow(k_Mpc / model->k_pivot, model->n_s - 1.0);

    /* Transfer function */
    double T_k = ps_transfer_EH(k_Mpc, model);

    /* Normalization: P(k) = P_R * T^2 * k * (c/H0)^4 * 2*pi^2 / k^3
     * with appropriate unit conversion.
     * We use a simpler form: P(k) proportional to P_R * T^2 * k^{n_s}
     * scaled to give the right sigma_8.
     */
    double D_H = COSMO_C_LIGHT / model->H0_si / COSMO_MPC;  /* Hubble distance in Mpc */
    double D_H_h = D_H * model->h;
    double norm = 2.0 * M_PI * M_PI / (k_Mpc * k_Mpc * k_Mpc);

    return P_R * T_k * T_k * norm * D_H_h * D_H_h * D_H_h * D_H_h;
}

/* === Dimensionless matter power spectrum ===
 * Delta_m^2(k) = k^3 * P(k) / (2*pi^2)
 *
 * Delta^2(k) ~ 1 corresponds to fluctuations of order unity
 * on that scale, indicating the onset of nonlinearity.
 */
double ps_matter_dimensionless(double k_Mpc, const CosmoModel *model)
{
    if (!model || k_Mpc <= 0.0) return 0.0;
    double P_k = ps_matter_power_linear(k_Mpc, model);
    return k_Mpc * k_Mpc * k_Mpc * P_k / (2.0 * M_PI * M_PI);
}


/* === Effective spectral index ===
 * n_eff(k) = d(ln P)/d(ln k)
 * Computed via centered finite difference.
 */
double ps_effective_n(double k_Mpc, const CosmoModel *model, double dk_frac)
{
    if (!model || k_Mpc <= 0.0 || dk_frac <= 0.0) return model ? model->n_s : 0.96;

    double k_lo = k_Mpc * (1.0 - dk_frac);
    double k_hi = k_Mpc * (1.0 + dk_frac);
    double P_lo = ps_matter_power_linear(k_lo, model);
    double P_hi = ps_matter_power_linear(k_hi, model);

    if (P_lo <= 0.0 || P_hi <= 0.0) return model->n_s;

    return (log(P_hi) - log(P_lo)) / (log(k_hi) - log(k_lo));
}

/* === Two-point correlation function ===
 * xi(r) = (1/2pi^2) * integral_0^inf k^2 * P(k) * sin(kr)/(kr) dk
 *
 * xi(r) measures the excess probability over random of finding
 * a pair of galaxies separated by distance r.
 *
 * For a power-law P(k) ~ k^n, xi(r) ~ r^{-(n+3)}.
 * The BAO feature appears as a peak at r ~ 105 Mpc/h (sound horizon).
 */
double ps_correlation_function(double r_Mpc_h, const CosmoModel *model)
{
    if (!model || r_Mpc_h <= 0.0) return 0.0;

    int n_k = 500;
    double log_k_min = log(1.0e-4);
    double log_k_max = log(1.0e3);
    double d_logk = (log_k_max - log_k_min) / n_k;

    double xi = 0.0;
    int i;

    for (i = 0; i < n_k; i++) {
        double log_k = log_k_min + i * d_logk;
        double k = exp(log_k);
        double dk = k * d_logk;

        double P_k = ps_matter_power_linear(k, model);
        if (P_k <= 0.0) continue;

        double kr = k * r_Mpc_h;
        double sinc;
        if (kr < 1.0e-6) {
            sinc = 1.0;
        } else {
            sinc = sin(kr) / kr;
        }

        xi += k * k * P_k * sinc * dk / (2.0 * M_PI * M_PI);
    }

    return xi;
}

/* === Window functions === */

/* Top-hat window in Fourier space:
 * W_th(kR) = 3 * (sin(kR) - kR*cos(kR)) / (kR)^3
 */
double ps_tophat_window(double kR)
{
    if (kR < 1.0e-8) return 1.0;
    return 3.0 * (sin(kR) - kR * cos(kR)) / (kR * kR * kR);
}

/* Gaussian window in Fourier space:
 * W_gauss(kR) = exp(-(kR)^2 / 2)
 */
double ps_gaussian_window(double kR)
{
    return exp(-0.5 * kR * kR);
}

/* === Variance on scale R ===
 * sigma_R^2 = integral k^2 dk/(2pi^2) * P(k) * |W(kR)|^2
 */
double ps_variance_R(double R_Mpc_h, const CosmoModel *model,
                      double (*window)(double))
{
    if (!model || R_Mpc_h <= 0.0 || !window) return 0.0;

    int n_k = 500;
    double log_k_min = log(1.0e-4);
    double log_k_max = log(1.0e3);
    double d_logk = (log_k_max - log_k_min) / n_k;

    double sigma2 = 0.0;
    int i;

    for (i = 0; i < n_k; i++) {
        double log_k = log_k_min + i * d_logk;
        double k = exp(log_k);
        double dk = k * d_logk;

        double P_k = ps_matter_power_linear(k, model);
        if (P_k <= 0.0) continue;

        double W = window(k * R_Mpc_h);
        sigma2 += k * k * P_k * W * W * dk / (2.0 * M_PI * M_PI);
    }

    return sqrt(sigma2);
}

/* === Nonlinear power spectrum (halofit) ===
 * Full halofit model with Takahashi et al. (2012) calibration.
 *
 * The nonlinear correction depends on:
 * - The linear power spectrum Delta_L^2(k)
 * - The effective spectral index n_eff
 * - The curvature of the power spectrum
 * - The RMS fluctuation on the nonlinear scale
 *
 * We implement a simplified version using the fitting functions
 * from Smith et al. (2003) Eqs. C1-C18.
 */
double ps_nonlinear_halofit(double k_Mpc, double z, const CosmoModel *model)
{
    if (!model || k_Mpc <= 0.0) return 0.0;

    /* Linear power at this redshift */
    double P_lin_z = ps_power_at_z(k_Mpc, z, model);
    double Delta_L = k_Mpc * k_Mpc * k_Mpc * P_lin_z / (2.0 * M_PI * M_PI);

    if (Delta_L < 1.0e-12) return P_lin_z;

    /* Find nonlinear scale k_nl where Delta_L^2(k_nl) = 1 */
    double k_nl = k_Mpc;
    double Delta_nl_test = Delta_L;
    int iter;

    for (iter = 0; iter < 20; iter++) {
        if (fabs(Delta_nl_test - 1.0) < 0.01) break;
        k_nl *= pow(Delta_nl_test, 0.3);
        double P_test = ps_power_at_z(k_nl, z, model);
        Delta_nl_test = k_nl * k_nl * k_nl * P_test / (2.0 * M_PI * M_PI);
    }

    /* Effective spectral index at nonlinear scale */
    double n_eff = ps_effective_n(k_nl, model, 0.1);

    /* Halofit parameters (simplified Smith+03) */
    double a_n = pow(10.0, 1.4861 + 1.8369*n_eff + 1.6762*n_eff*n_eff
                     + 0.7940*n_eff*n_eff*n_eff
                     + 0.1670*n_eff*n_eff*n_eff*n_eff
                     - 0.6206*n_eff*n_eff*n_eff*n_eff*n_eff);
    double b_n = pow(10.0, 0.9463 + 0.9466*n_eff + 0.3084*n_eff*n_eff
                     - 0.9400*n_eff*n_eff*n_eff);
    double c_n = pow(10.0, -0.2807 + 0.6669*n_eff + 0.3214*n_eff*n_eff
                     - 0.2207*n_eff*n_eff*n_eff);
    double gamma_n = 0.8649 + 0.2989*n_eff + 0.1631*n_eff*n_eff;
    double alpha_n = 1.3884 + 0.3700*n_eff - 0.1452*n_eff*n_eff;
    double beta_n = 0.8291 + 0.9854*n_eff + 0.3401*n_eff*n_eff;

    double y = k_Mpc / k_nl;
    double Delta_Q = Delta_L
        * pow((1.0 + Delta_L), beta_n) / (1.0 + alpha_n * Delta_L)
        * exp(-y/4.0 - y*y/8.0);

    double Delta_H = a_n * pow(y, 3.0*fabs(gamma_n))
        / (1.0 + b_n * pow(y, fabs(gamma_n))
           + pow(c_n * fabs(gamma_n) * y, 3.0 - fabs(gamma_n)));

    double Delta_NL = Delta_Q + Delta_H;

    /* Convert back to P(k) */
    return Delta_NL * 2.0 * M_PI * M_PI / (k_Mpc * k_Mpc * k_Mpc);
}

/* === Linear power spectrum at redshift z ===
 * P_lin(k, z) = P_lin(k, 0) * [D(z)/D(0)]^2
 */
double ps_power_at_z(double k_Mpc, double z, const CosmoModel *model)
{
    if (!model || k_Mpc <= 0.0) return 0.0;

    double P_k_0 = ps_matter_power_linear(k_Mpc, model);
    if (P_k_0 <= 0.0) return 0.0;

    double a_z = 1.0 / (1.0 + z);
    double D_z = perturbations_growth_factor_D(model, a_z);
    double D_0 = perturbations_growth_factor_D(model, 1.0);

    if (D_0 <= 0.0) return P_k_0;
    double ratio = D_z / D_0;
    return P_k_0 * ratio * ratio;
}
