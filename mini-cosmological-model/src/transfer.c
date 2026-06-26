/** @file transfer.c
 *  @brief Cosmological transfer functions: CMB anisotropy, Sachs-Wolfe
 *         effect, integrated Sachs-Wolfe, Doppler contribution,
 *         polarization, and lensing potential.
 *
 *  The CMB angular power spectrum C_l receives contributions from
 *  multiple physical effects on the last-scattering surface (LSS)
 *  and along the line of sight:
 *
 *  C_l = C_l^SW + C_l^ISW + C_l^Doppler + C_l^pol + ...
 *
 *  Reference: Dodelson & Schmidt (2020) Ch.9-10
 *             Hu & Dodelson (2002), ARAA 40, 171
 *             Seljak & Zaldarriaga (1996), ApJ 469, 437
 */

#include "cosmo_model.h"
#include "friedmann.h"
#include "distances.h"
#include "perturbations.h"
#include <math.h>

/* =========================================================================
 * Key cosmological scales
 * ========================================================================= */

/* Sound horizon at recombination:
 * r_s(z_rec) = integral_{z_rec}^{inf} c_s(z) dz / H(z)
 * where c_s = c / sqrt(3*(1+R)) is the sound speed,
 * R = 3*rho_b/(4*rho_gamma) = 3*Omega_b/(4*Omega_gamma) * a.
 *
 * This sets the physical scale of the BAO feature (~147 Mpc for Planck 2018).
 */
double transfer_sound_horizon(const CosmoModel *model, double z_rec)
{
    if (!model || z_rec < 0.0) return -1.0;

    int n = model->n_z_steps;
    if (n < 10) n = 1000;

    double z_max = 1.0e7;
    double log_zp1_min = log(1.0 + z_rec);
    double log_zp1_max = log(1.0 + z_max);
    double d_logzp1 = (log_zp1_max - log_zp1_min) / n;

    double R_eq = 3.0 * model->Omega_b / (4.0 * COSMO_OMEGA_R0);
    double integral = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double log_zp1_i = log_zp1_min + i * d_logzp1;
        double log_zp1_mid = log_zp1_i + 0.5 * d_logzp1;
        double log_zp1_next = log_zp1_i + d_logzp1;

        double zp1_i   = exp(log_zp1_i);
        double zp1_mid = exp(log_zp1_mid);
        double zp1_next = exp(log_zp1_next);

        double z_i   = zp1_i   - 1.0;
        double z_mid = zp1_mid - 1.0;
        double z_next = zp1_next - 1.0;

        /* Sound speed c_s^2 = 1 / [3*(1 + R*a)] where a = 1/(1+z) */
        double a_i   = 1.0 / zp1_i;
        double a_mid = 1.0 / zp1_mid;
        double a_next = 1.0 / zp1_next;

        double cs_i   = 1.0 / sqrt(3.0 * (1.0 + R_eq * a_i));
        double cs_mid = 1.0 / sqrt(3.0 * (1.0 + R_eq * a_mid));
        double cs_next = 1.0 / sqrt(3.0 * (1.0 + R_eq * a_next));

        double H_i   = friedmann_H_z_si(model, z_i);
        double H_mid = friedmann_H_z_si(model, z_mid);
        double H_next = friedmann_H_z_si(model, z_next);

        double f_i   = (H_i   > 0.0) ? cs_i   * COSMO_C_LIGHT / H_i   : 0.0;
        double f_mid = (H_mid > 0.0) ? cs_mid * COSMO_C_LIGHT / H_mid : 0.0;
        double f_next = (H_next > 0.0) ? cs_next * COSMO_C_LIGHT / H_next : 0.0;

        double dz = z_next - z_i;
        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return integral;  /* meters */
}

/* Damping scale at recombination (Silk damping):
 * k_D^{-2} = integral c_s^2 * (16/15 + R^2/(1+R)) / (H * a^2 * tau_c) da
 * where tau_c = 1/(n_e * sigma_T * c) is the Thomson scattering time.
 *
 * Simplified form: k_D ~ 0.14 * Omega_m^0.15 * Omega_b^0.30 * h^0.70 Mpc^{-1}
 */
double transfer_damping_scale(const CosmoModel *model)
{
    if (!model) return -1.0;

    /* Fitting formula from Hu & White (1997) */
    double Om_m = model->Omega_m;
    double Om_b = model->Omega_b;
    double h = model->h;

    double k_D = 0.14 * pow(Om_m, 0.15) * pow(Om_b, 0.30) * pow(h, 0.70);

    return k_D;  /* Mpc^{-1} */
}

/* =========================================================================
 * Sachs-Wolfe (SW) effect
 * ========================================================================= */

/* The Sachs-Wolfe effect relates CMB temperature fluctuations to
 * gravitational potential perturbations on the last scattering surface.
 *
 * Delta_T/T = -(1/3) * Phi  (for adiabatic perturbations)
 *
 * On large scales (l < 100), the SW effect dominates.
 * C_l^SW ~ (2/pi) * (A_s/25) * (Gamma((3-n_s)/2) * Gamma(l+(n_s-1)/2))
 *         / (Gamma((5-n_s)/2) * Gamma(l+2-(n_s-1)/2))
 *
 * For a scale-invariant spectrum (n_s=1):
 *   l*(l+1)*C_l^SW / (2*pi) = A_s / 25 = constant
 */

/* Compute the Sachs-Wolfe contribution to C_l^TT for scale-invariant spectrum.
 * l*(l+1)*C_l/(2*pi) = A_s * pi^2 / 25 * (Gamma(3/2) * Gamma(l+1/2))
 *                     / (Gamma(1) * Gamma(l+3/2))
 *
 * For large l: C_l ~ A_s * pi / (25 * l*(l+1)) -> constant l(l+1)C_l.
 */
double transfer_Cl_SW(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;

    double A_s = model->A_s;
    double n_s = model->n_s;

    /* Gamma function ratio using Lanczos approximation */
    /* For scale-invariant (n_s=1), use the closed form: */
    if (fabs(n_s - 1.0) < 1.0e-4) {
        return A_s * M_PI / (25.0 * l * (l + 1.0));
    }

    /* For tilted spectrum (n_s != 1), the analytical expression involves
     * Gamma functions. We compute a simplified form. */
    double gamma_ratio = tgamma(0.5 * (3.0 - n_s)) * tgamma(l + 0.5*(n_s - 1.0))
                       / (tgamma(0.5 * (5.0 - n_s)) * tgamma(l + 2.0 - 0.5*(n_s - 1.0)));

    return (2.0/M_PI) * (A_s/25.0) * gamma_ratio;
}

/* Compute the temperature fluctuation amplitude for mode l.
 * l(l+1)*C_l^TT/(2*pi) [microKelvin^2] */
double transfer_Cl_TT_SW_uK2(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;
    double Cl = transfer_Cl_SW(model, l);
    /* Convert to microKelvin^2: multiply by T_cmb^2 */
    double T_cmb_uK = model->T_cmb0 * 1.0e6;
    return l * (l + 1.0) * Cl * T_cmb_uK * T_cmb_uK / (2.0 * M_PI);
}

/* =========================================================================
 * Integrated Sachs-Wolfe (ISW) effect
 * ========================================================================= */

/* The ISW effect arises from the time evolution of gravitational potentials
 * as CMB photons traverse evolving large-scale structure.
 *
 * Delta_T/T = 2 * integral (dPhi/dt) dt  (along the line of sight)
 *
 * ISW is significant at late times (dark energy domination) and
 * at very early times (radiation-matter transition: early ISW).
 *
 * The ISW contribution peaks at low l (l < 30) and is correlated
 * with large-scale structure tracers (ISW-galaxy cross-correlation).
 */

/* Compute ISW kernel: the weight function for the ISW integral.
 * W_ISW(z) = 3*Omega_m*H0^2/(c^3*k^2) * d(D/a)/dz
 * where D is the growth factor and a is the scale factor.
 */
double transfer_isw_kernel(const CosmoModel *model, double z)
{
    if (!model) return 0.0;

    double a = 1.0 / (1.0 + z);
    double D = perturbations_growth_factor_D(model, a);
    if (D <= 0.0) return 0.0;

    double D_over_a = D / a;

    /* Derivative d(D/a)/dz via finite difference */
    double dz = 0.01;
    double zp = z + dz;
    double ap = 1.0 / (1.0 + zp);
    double Dp = perturbations_growth_factor_D(model, ap);
    double D_over_a_p = Dp / ap;

    double d_Dovera_dz = (D_over_a_p - D_over_a) / dz;

    double prefactor = 3.0 * model->Omega_m * model->H0_si * model->H0_si
                      / (COSMO_C_LIGHT * COSMO_C_LIGHT * COSMO_C_LIGHT);

    return prefactor * fabs(d_Dovera_dz);
}

/* ISW contribution to C_l^TT:
 * C_l^ISW ~ integral dz * W_ISW^2(z) * P_Phi(k=l/chi(z), z)
 * where chi(z) is the comoving distance to redshift z.
 *
 * Simplified: C_l^ISW ~ A_ISW * exp(-l/l_ISW)
 */
double transfer_Cl_ISW(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;

    /* ISW amplitude scales with the amount of dark energy */
    double A_ISW = model->Omega_Lambda * model->Omega_Lambda
                 * (model->A_s / 25.0) * M_PI;

    /* ISW decays exponentially above l ~ 30 */
    double l_ISW = 30.0;

    return A_ISW * exp(-(double)l / l_ISW) / (l * (l + 1.0));
}

/* =========================================================================
 * Doppler effect
 * ========================================================================= */

/* The Doppler contribution from baryon velocity at recombination
 * produces peaks in the CMB power spectrum at l ~ 200-500.
 *
 * C_l^Doppler ~ (Omega_b*h^2)^2 * sin^2(k*r_s) * exp(-k^2/k_D^2)
 *                 * j_l'(k*D_rec)^2
 *
 * The Doppler peaks are out of phase with the acoustic (compression)
 * peaks by pi/2, causing the alternating peak heights observed in C_l.
 */

/* Approximate Doppler contribution for mode l using a semi-analytic form.
 * The Doppler term is subdominant to the acoustic term for most l,
 * except at the trough between the first and second acoustic peaks. */
double transfer_Cl_Doppler(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;

    /* Doppler amplitude scales with baryon density */
    double A_Dop = model->Omega_b * model->Omega_b * model->h * model->h
                 * (model->A_s / 25.0);

    /* Doppler roughly constant in l(l+1)C_l across the acoustic peaks */
    double k_eff = (double)l / 100.0;  /* Effective wavenumber */

    double r_s = COSMO_R_DRAG;
    double k_D = transfer_damping_scale(model);
    if (k_D <= 0.0) k_D = 0.14;

    double damping = exp(-k_eff * k_eff / (k_D * k_D));
    double osc = sin(k_eff * r_s) * sin(k_eff * r_s);

    return A_Dop * damping * osc / (l * (l + 1.0));
}


/* =========================================================================
 * Polarization contributions
 * ========================================================================= */

/* CMB polarization is generated by Thomson scattering of anisotropic
 * radiation at recombination. The E-mode polarization power spectrum
 * C_l^EE has peaks at the same acoustic scales as C_l^TT, but the
 * first peak is at l ~ 300 and polarization peaks are out of phase
 * with temperature peaks.
 *
 * B-mode polarization is generated only by tensor perturbations
 * (primordial gravitational waves) and gravitational lensing,
 * making its detection a key probe of inflation.
 */

/* Simplified E-mode polarization spectrum.
 * In the tight-coupling approximation:
 * C_l^EE ~ C_l^TT * (Delta_tau / tau_rec)^2 * (l/l_pol)^2
 * where Delta_tau is the width of the last-scattering surface.
 */
double transfer_Cl_EE(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;

    /* Polarization fraction at recombination */
    double pol_fraction = 0.01;   /* ~1% polarization at recombination */
    double l_pol_peak = 300.0;    /* First polarization peak */

    /* Quadratic rise at low l */
    double l_factor = (l / l_pol_peak) * (l / l_pol_peak);

    /* Get the temperature spectrum for the acoustic structure */
    double Cl_TT = transfer_Cl_SW(model, l);

    return pol_fraction * Cl_TT * l_factor;
}

/* =========================================================================
 * CMB lensing potential
 * ========================================================================= */

/* Gravitational lensing of the CMB by large-scale structure smooths
 * the acoustic peaks and generates B-mode polarization from E-modes.
 *
 * The lensing potential power spectrum C_l^phiphi is:
 * C_l^phiphi = 8*pi^2 * l^{-4}
 *   * integral dz (chi/chi_rec) * (1 - chi/chi_rec)^2
 *   * P_Phi(k=l/chi, z) * dchi/dz
 *
 * where chi is comoving distance, chi_rec is distance to LSS,
 * and P_Phi is the power spectrum of the gravitational potential.
 */

/* Compute the CMB lensing convergence power spectrum.
 * C_l^kk = l^2*(l+1)^2 * C_l^phiphi / 4
 *
 * Simplified form calibrated against Planck 2018.
 */
double transfer_Cl_lensing_kk(const CosmoModel *model, int l)
{
    if (!model || l < 2) return 0.0;

    /* Lensing amplitude (Planck 2018: A_L = 1.180 +/- 0.065) */
    double A_L = 1.18;

    /* Lensing kernel peaks at z ~ 2 */
    double z_lens = 2.0;
    double D_lens = distances_comoving_line_of_sight(model, z_lens);
    double D_rec  = distances_comoving_line_of_sight(model, 1089.0);
    if (D_lens <= 0.0 || D_rec <= 0.0) return 0.0;

    /* Lens efficiency: chi/chi_rec * (1 - chi/chi_rec) */
    double chi_ratio = D_lens / D_rec;
    double efficiency = chi_ratio * (1.0 - chi_ratio);

    /* Approximate: C_l^kk scales as l^2 * C_l^phiphi */
    double Cl_phiphi = A_L * model->A_s * M_PI * M_PI
                     * pow(model->H0_si * D_rec / COSMO_C_LIGHT, 2.0)
                     * efficiency * efficiency
                     / (l * l * l * l);

    /* Convert to kappa: C_l^kk = l^2(l+1)^2/4 * C_l^phiphi */
    double l_factor = l * (l + 1.0);
    return l_factor * l_factor / 4.0 * Cl_phiphi;
}

/* =========================================================================
 * Matter transfer function (Poisson equation)
 * ========================================================================= */

/* The gravitational potential Phi is related to the matter density
 * contrast delta via the Poisson equation:
 *   k^2 * Phi = -(3/2) * Omega_m * H0^2 * delta / a
 *
 * The transfer function T(k) relates the primordial curvature
 * perturbation R(k) to the late-time density contrast:
 *   delta(k, a) = (2/5) * (k^2 * T(k) * R(k)) / (Omega_m * H0^2)
 *               * D(a) / a
 *
 * This relationship is encoded in the ratio between the
 * potential and matter transfer functions.
 */

/* Compute P(k) from the gravitational potential power spectrum P_Phi(k):
 * P_delta(k) = (4/9) * k^4 / (Omega_m^2 * H0^4) * P_Phi(k) * T^2(k)
 */
double transfer_potential_to_density(double k_Mpc, double P_Phi_k,
                                      double Omega_m, double H0)
{
    if (k_Mpc <= 0.0 || Omega_m <= 0.0 || H0 <= 0.0) return 0.0;

    double H0_Mpc = H0 / COSMO_MPC;  /* H0 in Mpc^{-1} (but H0 is km/s/Mpc) */
    double prefactor = 4.0 / 9.0 * k_Mpc * k_Mpc * k_Mpc * k_Mpc
                     / (Omega_m * Omega_m * H0_Mpc * H0_Mpc * H0_Mpc * H0_Mpc);

    return prefactor * P_Phi_k;
}

/* =========================================================================
 * Growth factor transfer function
 * ========================================================================= */

/* The ratio of the growth factor today to the growth factor
 * during matter domination encodes the suppression of growth
 * by dark energy.
 *
 * g(z) = (1+z) * D(z)
 *
 * During matter domination, g(z) = constant.
 * During Lambda domination, g(z) decreases as growth is suppressed.
 * The suppression factor at z=0 relative to z_high is:
 *   g_0 = D_0 / a_0 = D_0 (since a_0 = 1)
 *
 * For LCDM: g_0 ~ 0.78 (growth suppressed ~22% relative to EdS)
 */

/* Compute the growth suppression factor g(a) = D(a)/a.
 * For Einstein-de Sitter: g = 1 always.
 * For LCDM: g < 1 at late times due to dark energy. */
double transfer_growth_suppression(const CosmoModel *model, double a)
{
    if (!model || a <= 0.0) return 1.0;

    double D_a = perturbations_growth_factor_D(model, a);
    if (D_a <= 0.0) return 1.0;

    return D_a / a;
}

/* =========================================================================
 * f*sigma_8: the product of growth rate and clustering amplitude.
 * This is a key observable from redshift-space distortions (RSD).
 *
 * f*sigma_8(z) = f(z) * sigma_8 * D(z)/D(0)
 *
 * Planck 2018: f*sigma_8(z=0) = 0.421 +/- 0.039
 *               sigma_8 = 0.811 +/- 0.006
 */

/* Compute f*sigma_8 at redshift z */
double transfer_f_sigma8(const CosmoModel *model, double z)
{
    if (!model) return -1.0;

    double a = 1.0 / (1.0 + z);
    double f = perturbations_growth_rate_f(model, a);
    double D_z = perturbations_growth_factor_D(model, a);
    double D_0 = perturbations_growth_factor_D(model, 1.0);

    if (D_0 <= 0.0) return -1.0;

    /* sigma_8 is roughly sigma(8 Mpc/h) */
    double sigma_8 = perturbations_sigma_R(model, 8.0);

    return f * sigma_8 * D_z / D_0;
}

/* =========================================================================
 * S8 parameter: sigma_8 * sqrt(Omega_m/0.3)
 *
 * S8 is the parameter combination best constrained by weak lensing surveys.
 * It is less degenerate than sigma_8 alone.
 *
 * Planck 2018: S8 = 0.832 +/- 0.013
 * KiDS-1000:   S8 = 0.759 +/- 0.024  (2-3 sigma tension with Planck)
 */

/* Compute S8 parameter */
double transfer_S8(const CosmoModel *model)
{
    if (!model) return -1.0;

    double sigma_8 = perturbations_sigma_R(model, 8.0);
    if (sigma_8 <= 0.0) return -1.0;

    return sigma_8 * sqrt(model->Omega_m / 0.3);
}

/* =========================================================================
 * Neutrino free-streaming scale
 * ========================================================================= */

/* Massive neutrinos suppress the matter power spectrum on scales
 * smaller than their free-streaming scale:
 *
 * k_fs(z) ~ 0.8 * sqrt(Omega_m*(1+z)) / (1+z)
 *           * (m_nu/1eV) * h Mpc^{-1}
 *
 * The suppression is ~ -8*f_nu in power, where f_nu = Omega_nu/Omega_m
 * is the neutrino fraction.
 */

/* Compute the neutrino free-streaming wavenumber */
double transfer_neutrino_free_streaming(const CosmoModel *model, double z)
{
    if (!model) return -1.0;

    double Omega_nu = model->sum_m_nu / 93.14 / (model->h * model->h);
    if (Omega_nu <= 0.0) return INFINITY;

    double k_fs = 0.8 * sqrt(model->Omega_m * (1.0 + z))
                / (1.0 + z)
                * (model->sum_m_nu / 1.0)
                * model->h;

    return k_fs;  /* Mpc^{-1} */
}

/* =========================================================================
 * Reionization optical depth effect on CMB
 * ========================================================================= */

/* Thomson scattering by free electrons after reionization suppresses
 * the CMB temperature power spectrum at all l by a factor e^{-2*tau}.
 * The reionization bump at very low l (l < 10) is generated by
 * large-scale polarization generated during reionization.
 *
 * C_l^TT -> C_l^TT * e^{-2*tau}  (for l >> 10)
 * C_l^TT -> additional contribution from reionization at very low l
 */

/* Apply reionization suppression to a raw C_l value */
double transfer_apply_reionization(double Cl_raw, double tau)
{
    return Cl_raw * exp(-2.0 * tau);
}

/* Reionization contribution at very low l (l < 10).
 * This is an additional source of large-angle polarization. */
double transfer_Cl_reionization_bump(const CosmoModel *model, int l)
{
    if (!model || l > 20) return 0.0;

    double tau = model->tau_reio;
    double z_reion = 7.5;  /* Approximate reionization redshift */
    (void)z_reion;  /* Used in full CMB calculation with reionization modeling */

    /* Amplitude proportional to tau^2 and the optical depth */
    double amp = tau * tau * model->A_s * M_PI / 4.0;

    /* Only significant at l < 10 */
    return amp * exp(-0.5 * l * l / 4.0) / (l * (l + 1.0));
}

