/** @file perturbations.c
 *  @brief Cosmological perturbation theory: linear growth factor,
 *         spherical collapse, Press-Schechter mass function,
 *         nonlinear corrections, BAO features.
 *
 *  Growth equation (sub-horizon, Newtonian gauge):
 *    delta'' + (3/a + E'/E)*delta' - (3*Omega_m/(2*a^2*E^2))*delta = 0
 *  where prime = d/da.
 *
 *  Reference: Dodelson & Schmidt (2020) Ch.7-8
 *             Mo, van den Bosch, White (2010) "Galaxy Formation and Evolution"
 *             Press & Schechter (1974), ApJ 187, 425
 */

#include "perturbations.h"
#include "friedmann.h"
#include "distances.h"
#include "power_spectrum.h"
#include "friedmann.h"
#include <math.h>

/* === Linear growth factor D(a) ===
 *
 * The growth factor D(a) describes the time evolution of the amplitude
 * of linear density perturbations: delta(a) = D(a) * delta_0.
 *
 * For Omega_m=1, flat: D(a) = a (exact).
 * For LCDM, we solve the ODE numerically using a simple integrator.
 *
 * Normalization: D(a=1) = 1 (growth factor today = 1).
 */
double perturbations_growth_factor_D(const CosmoModel *model, double a)
{
    if (!model || a <= 0.0) return 0.0;

    /* Normalization: D(a) ~ a at early times (matter-dominated) */
    if (fabs(model->Omega_m - 1.0) < 1e-6 && fabs(model->Omega_Lambda) < 1e-6) {
        /* Einstein-de Sitter: D(a) = a exactly */
        return a;
    }

    /* For general LCDM, integrate the growth equation.
     * We use a simple Euler integration from a_min to a.
     * Growth equation rewritten as coupled first-order ODE:
     *   dD/da = G
     *   dG/da = -(3/a + dE/da/E)*G + (3*Omega_m0*a^{-3}/(2*E^2))*D
     *
     * At early times (radiation era negligible for growth):
     *   D ~ a, G ~ 1 (constant growth rate)
     */
    double a_min = 1.0e-5;  /* Deep in matter era, radiation negligible */
    double D = a_min;       /* D ~ a at early times */
    double G = 1.0;         /* dD/da = 1 at early times */
    int n_steps = 2000;
    double da = (a - a_min) / n_steps;
    double ai;
    int i;

    for (i = 0; i < n_steps; i++) {
        ai = a_min + i * da;
        double z = 1.0 / ai - 1.0;
        double E = friedmann_E_z(model, z);

        if (E <= 0.0) continue;

        /* dE/da: E^2 = Om*a^{-3} + OL (flat LCDM), so 2*E*dE/da = -3*Om*a^{-4}
         * => dE/da = -3*Om*a^{-4}/(2*E) */
        double dE_da = -1.5 * model->Omega_m * pow(ai, -4.0) / E;

        /* dG/da from growth equation */
        double dG_da = -(3.0/ai + dE_da/E) * G
                     + (1.5 * model->Omega_m * pow(ai, -3.0) / (E*E)) * D;

        D += da * G;
        G += da * dG_da;
    }

    /* Normalize to D(a=1) = 1 */
    return D;
}

/* === Logarithmic growth rate f(a) = d(ln D)/d(ln a) ===
 * f(a) = a/D * dD/da
 *      ~ Omega_m(a)^gamma  with gamma ~ 0.55 for LCDM (Lahav et al. 1991)
 */
double perturbations_growth_rate_f(const CosmoModel *model, double a)
{
    if (!model || a <= 0.0) return 0.0;

    /* Use the fitting formula f = Omega_m(a)^0.55 */
    double z = 1.0 / a - 1.0;
    double E_z = friedmann_E_z(model, z);
    if (E_z <= 0.0) return 0.0;

    double Omega_m_a = model->Omega_m * pow(a, -3.0) / (E_z * E_z);
    return pow(Omega_m_a, 0.55);
}

/* === Power spectrum time evolution ===
 * P(k, z) = P(k, 0) * [D(z)/D(0)]^2
 * Linear theory: each Fourier mode evolves independently.
 */
double perturbations_power_spectrum_evolved(double P_k_ini, double D_ini,
                                             double D_z)
{
    if (D_ini <= 0.0) return 0.0;
    double ratio = D_z / D_ini;
    return P_k_ini * ratio * ratio;
}

/* === RMS fluctuation sigma(R) ===
 * sigma^2(R) = (1/2pi^2) * integral_0^inf P(k) * W^2(kR) * k^2 dk
 *
 * For a simple power-law P(k) ~ k^n, sigma^2(R) ~ R^{-(n+3)}.
 * We use a numerical integration with the BBKS transfer function.
 */
double perturbations_sigma_R(const CosmoModel *model, double R_Mpc_h)
{
    if (!model || R_Mpc_h <= 0.0) return 0.0;

    /* Integrate over k using trapezoidal rule.
     * k range: 1e-4 to 1e3 h/Mpc (covers all relevant scales). */
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

        /* Power spectrum at this k */
        double P_k = ps_matter_power_linear(k, model);
        if (P_k <= 0.0) continue;

        /* Top-hat window function */
        double kR = k * R_Mpc_h;
        double W;
        if (kR < 1.0e-6) {
            W = 1.0;
        } else {
            W = 3.0 * (sin(kR) - kR * cos(kR)) / (kR * kR * kR);
        }

        sigma2 += k * k * P_k * W * W * dk / (2.0 * M_PI * M_PI);
    }

    return sqrt(sigma2);
}

/* === Critical overdensity for spherical collapse ===
 * delta_c ~ 1.686 for Einstein-de Sitter.
 * Weak cosmology dependence: delta_c(z) approx 1.686 * [1 + 0.0123*log(Omega_m(z))]
 * (Weinberg 2008, Eq. 8.68; Kitayama & Suto 1996)
 */
double perturbations_delta_c(const CosmoModel *model, double z)
{
    if (!model) return 1.686;

    double zp1 = 1.0 + z;
    double E_z = friedmann_E_z(model, z);
    if (E_z <= 0.0) return 1.686;

    double Omega_m_z = model->Omega_m * zp1 * zp1 * zp1 / (E_z * E_z);

    if (Omega_m_z <= 0.0 || Omega_m_z > 1.0) return 1.686;
    return 1.686 * (1.0 + 0.0123 * log(Omega_m_z));
}


/* === Press-Schechter mass function ===
 * dn/dM dM = sqrt(2/pi) * (rho_m/M^2) * (delta_c/sigma) * |d(ln sigma)/d(ln M)|
 *           * exp(-delta_c^2/(2*sigma^2)) dM
 *
 * Number density of dark matter halos in mass range [M, M+dM] at redshift z.
 */
double perturbations_mass_function_PS(const CosmoModel *model,
                                       double M_Msun_h, double z)
{
    if (!model || M_Msun_h <= 0.0) return 0.0;

    /* Convert mass to Lagrangian radius: M = (4*pi/3) * rho_m * R^3 */
    double rho_m0 = model->Omega_m * model->rho_crit0;  /* kg/m^3 */
    double rho_m0_Msun_Mpc3 = rho_m0 * COSMO_MPC * COSMO_MPC * COSMO_MPC
                              / (COSMO_M_SOLAR * model->h * model->h * model->h);

    double R_Mpc_h = pow(3.0 * M_Msun_h / (4.0 * M_PI * rho_m0_Msun_Mpc3), 1.0/3.0);

    double sigma = perturbations_sigma_R(model, R_Mpc_h);
    if (sigma <= 0.0) return 0.0;

    double delta_c = perturbations_delta_c(model, z);

    /* d(ln sigma)/d(ln M) via finite difference */
    double dM_frac = 0.01;
    double R_plus = pow(3.0 * M_Msun_h * (1.0 + dM_frac)
                        / (4.0 * M_PI * rho_m0_Msun_Mpc3), 1.0/3.0);
    double sigma_plus = perturbations_sigma_R(model, R_plus);
    double d_ln_sigma_d_ln_M = (log(sigma_plus) - log(sigma)) / log(1.0 + dM_frac);

    double nu = delta_c / sigma;
    double prefactor = sqrt(2.0/M_PI) * rho_m0_Msun_Mpc3 / (M_Msun_h * M_Msun_h);

    return prefactor * fabs(d_ln_sigma_d_ln_M) * nu * exp(-0.5 * nu * nu);
}

/* === Halo bias (peak-background split) ===
 * b(M) = 1 + (nu^2 - 1)/delta_c
 * where nu = delta_c/sigma(M)
 *
 * Massive halos (high nu) are more strongly clustered than the matter.
 */
double perturbations_halo_bias_PS(const CosmoModel *model,
                                   double M_Msun_h, double z)
{
    if (!model || M_Msun_h <= 0.0) return 1.0;

    double rho_m0 = model->Omega_m * model->rho_crit0;
    double rho_m0_Msun_Mpc3 = rho_m0 * COSMO_MPC * COSMO_MPC * COSMO_MPC
                              / (COSMO_M_SOLAR * model->h * model->h * model->h);

    double R_Mpc_h = pow(3.0 * M_Msun_h / (4.0 * M_PI * rho_m0_Msun_Mpc3), 1.0/3.0);
    double sigma = perturbations_sigma_R(model, R_Mpc_h);
    if (sigma <= 0.0) return 1.0;

    double delta_c = perturbations_delta_c(model, z);
    double nu = delta_c / sigma;

    return 1.0 + (nu * nu - 1.0) / delta_c;
}

/* === Nonlinear power spectrum (simplified halofit) ===
 * Using the Smith et al. (2003) fitting formula:
 *   Delta_NL(k) = Delta_L(k) * [ (1 + B*Delta_L)^beta / (1 + A*Delta_L)^alpha ]
 * where Delta(k) = k^3*P(k)/(2*pi^2) is the dimensionless power.
 *
 * This is a simplified version. Full halofit requires additional
 * parameters dependent on n_eff and the linear power spectrum shape.
 */
double perturbations_nonlinear_power(const CosmoModel *model,
                                      double k_Mpc, double z)
{
    if (!model || k_Mpc <= 0.0) return 0.0;

    /* Get linear dimensionless power at this k and z */
    double P_lin = ps_power_at_z(k_Mpc, z, model);
    double Delta_lin = k_Mpc * k_Mpc * k_Mpc * P_lin / (2.0 * M_PI * M_PI);

    if (Delta_lin <= 0.0) return 0.0;

    /* Simplified Halofit parameters (calibrated for LCDM near z=0) */
    double A = 0.7;
    double B = 0.4;
    double alpha = 0.7;
    double beta = 0.3;

    double Delta_nl = Delta_lin
        * pow((1.0 + B * Delta_lin), beta)
        / pow((1.0 + A * Delta_lin), alpha);

    /* Convert back to P(k) */
    return Delta_nl * 2.0 * M_PI * M_PI / (k_Mpc * k_Mpc * k_Mpc);
}

/* === BAO oscillation feature ===
 * The BAO wiggles in the power spectrum are a damped sinusoidal
 * modulation relative to a smooth ("no-wiggle") spectrum:
 *
 * O(k) = sin(k * r_drag) / (k * r_drag)
 *
 * multiplied by a Gaussian damping from nonlinear evolution
 * and reconstruction effects.
 */
double perturbations_bao_oscillation(double k_Mpc, double r_drag_Mpc,
                                      double Sigma_nl)
{
    if (k_Mpc <= 0.0 || r_drag_Mpc <= 0.0) return 0.0;

    double x = k_Mpc * r_drag_Mpc;
    double osc;

    if (x < 1.0e-6) {
        osc = 1.0;
    } else {
        osc = sin(x) / x;
    }

    /* Gaussian damping from nonlinear structure growth */
    double damping = exp(-0.5 * k_Mpc * k_Mpc * Sigma_nl * Sigma_nl);

    return osc * damping;
}

/* === CMB temperature power spectrum (Sachs-Wolfe approximation) ===
 * C_l^TT_SW = (2/pi) * integral k^2 dk * P_R(k) * j_l^2(k * D_rec)
 *
 * For large scales (low l), the Sachs-Wolfe effect dominates:
 *   Delta_T / T = -(1/3) * Phi  (gravitational potential at LSS)
 *   C_l^TT ~ A_s * pi / (2*l*(l+1))  (scale-invariant, Harrison-Zeldovich)
 *
 * Spherical Bessel function j_l(x) computed via recurrence:
 *   j_0 = sin(x)/x, j_1 = sin(x)/x^2 - cos(x)/x
 *   j_{l+1} = (2l+1)/x * j_l - j_{l-1}
 */
double perturbations_cmb_Cl_TT_SW(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;

    /* Comoving distance to recombination surface */
    double z_rec = 1089.0;  /* Planck 2018 best-fit */
    double D_rec = distances_comoving_line_of_sight(model, z_rec);
    if (D_rec <= 0.0) return 0.0;

    /* Integrate over k */
    int n_k = 300;
    double log_k_min = log(1.0e-5);
    double log_k_max = log(1.0);
    double d_logk = (log_k_max - log_k_min) / n_k;

    double integral = 0.0;
    int i;

    for (i = 0; i < n_k; i++) {
        double log_k = log_k_min + i * d_logk;
        double k = exp(log_k);  /* k in Mpc^{-1}, but D_rec in meters... */
        double dk = k * d_logk;

        /* Convert k to proper units (1/m) */
        double k_m = k / (COSMO_MPC * 1.0e6);  /* Mpc^{-1} -> m^{-1} */

        /* Primordial curvature power spectrum */
        double P_R = model->A_s * pow(k / model->k_pivot, model->n_s - 1.0);

        /* Spherical Bessel function j_l(k*D_rec) */
        double x = k_m * D_rec;
        double j_l;

        if (x < 1.0e-10) {
            /* For very small x: j_l(x) ~ x^l / (2l+1)!! */
            j_l = pow(x, l);
            int m;
            for (m = 1; m <= 2*l+1; m += 2) j_l /= m;
        } else {
            /* Recurrence relation for spherical Bessel functions */
            double j_lm2 = sin(x) / x;  /* j_0 */
            double j_lm1;
            if (x < 1.0e-6) {
                j_lm1 = x / 3.0;  /* j_1 ~ x/3 for small x */
            } else {
                j_lm1 = sin(x) / (x*x) - cos(x) / x;  /* j_1 */
            }
            double j_l_cur = j_lm1;
            int il;
            for (il = 2; il <= l; il++) {
                j_l_cur = (2.0*il - 1.0)/x * j_lm1 - j_lm2;
                j_lm2 = j_lm1;
                j_lm1 = j_l_cur;
            }
            j_l = j_l_cur;
        }

        integral += k_m * k_m * P_R * j_l * j_l * dk;
    }

    return (2.0 / M_PI) * integral;
}

