/**
 * @file    gw_stochastic.c
 * @brief   Stochastic GW Background — Omega_GW, Cross-Correlation, Source Models
 *
 * L2 — Fractional energy density Omega_GW(f)
 * L4 — Cosmological background models (inflation, phase transitions, strings)
 * L6 — Astrophysical CBC background
 * L7 — Overlap reduction function, SGWB detection SNR
 */

#include "gw_stochastic.h"
#include "gw_detection.h"  /* for gw_psd_aligo_design, gw_psd_lisa */
#include <math.h>

/* ================================================================
 * L2 — Omega_GW Energy Density
 *
 * Omega_GW(f) = (1/rho_c) * d(rho_GW)/d(ln f)
 *
 * Power-law parameterization:
 *   Omega_GW(f) = Omega_ref * (f / f_ref)^alpha
 * ================================================================ */

double gw_omega_gw(double f, const GwStochasticParams *params) {
    if (f <= 0.0 || params->f_ref <= 0.0) return 0.0;
    return params->Omega_ref * pow(f / params->f_ref, params->alpha);
}

double gw_energy_density_from_omega(double Omega,
    double f_low, double f_high, double H0) {
    double rho_c = gw_critical_density(H0);
    if (f_low <= 0.0 || f_high <= f_low) return 0.0;
    double ln_ratio = log(f_high / f_low);
    return rho_c * Omega * ln_ratio;
}

double gw_critical_density(double H0) {
    /* rho_c = 3 H_0^2 / (8 pi G)
     * H_0 in [km/s/Mpc] -> convert to [1/s] first */
    double H0_si = gw_H0_to_per_sec(H0);
    return 3.0 * H0_si * H0_si / (8.0 * GW_PI * GW_G);
}

double gw_H0_to_per_sec(double H0_km_s_Mpc) {
    /* 1 km/s/Mpc = 1000 / 3.0857e22 = 3.2408e-20 s^{-1} */
    return H0_km_s_Mpc * 1000.0 / GW_MPC;
}

/* ================================================================
 * L6 — Astrophysical CBC Background
 *
 * Omega_CBC(f) = Omega_ref * (f / 25 Hz)^{2/3}
 *
 * The 2/3 spectral index comes from inspiral evolution:
 * dE/df ~ f^{-1/3} => Omega(f) ~ f^{2/3}
 * ================================================================ */

double gw_omega_cbc_background(double f, double Omega_ref) {
    double f_ref = 25.0;  /* reference frequency [Hz] */
    if (f <= 0.0) return 0.0;
    return Omega_ref * pow(f / f_ref, 2.0/3.0);
}

double gw_omega_cbc_integrated(double f, double A_ref) {
    return gw_omega_cbc_background(f, A_ref);
}

/* ================================================================
 * L6 — Cosmological Background Models
 * ================================================================ */

double gw_omega_inflation(double f, double r, double A_s) {
    /* Inflationary GW background.
     *
     * Omega_GW(f) ~ (r A_s / 12) * Omega_m * (f/f_eq)^{n_t}
     *                                / (1 + c f/f_eq)
     *
     * Simplified: Omega_GW ~ 1e-15 * (r / 0.1) at LIGO frequencies
     * (unobservably small for r < 1).
     */

    if (f <= 0.0) return 0.0;

    double n_t = -r / 8.0;  /* consistency relation */
    double f_pivot = 0.05 * 1.0e6 * GW_MPC / (3.0857e22);  /* 0.05 Mpc^{-1} in Hz ~ 4.9e-19 Hz */
    /* Actually for LIGO band: inflationary background is ~ flat with tiny amplitude */
    double Omega_0 = 1.0e-15 * (r / 0.1) * (A_s / 2.1e-9);

    return Omega_0 * pow(f / 1.0e-16, n_t);
}

double gw_omega_phase_transition(double f,
    double Omega_peak, double f_peak,
    double a, double b, double c) {
    /* Broken power-law spectral shape for first-order phase transition.
     *
     * S(x) = (a+b)^c / (b x^{-a/c} + a x^{b/c})^c
     *
     * where x = f / f_peak.
     */
    if (f <= 0.0 || f_peak <= 0.0) return 0.0;

    double x = f / f_peak;

    /* Guard against degenerate x values */
    if (x < 1e-100) x = 1e-100;
    if (x > 1e100)  x = 1e100;

    double x_ma = pow(x, -a/c);
    double x_pb = pow(x,  b/c);

    double num = pow(a + b, c);
    double den = pow(b * x_ma + a * x_pb, c);

    return Omega_peak * num / den;
}

double gw_omega_cosmic_strings(double f, double G_mu) {
    /* Cosmic string GW background.
     *
     * Omega_cs(f) ~ (128 pi / 9) * G_mu * Omega_r^{1/2}
     *               * sum_n n^{-4/3} * ...
     *
     * Approximate for LIGO band: nearly flat.
     * Omega_cs ~ few × 10^{-10} (G_mu / 10^{-8}) at f ~ 25 Hz.
     *
     * Current LVK bound: G_mu < 4 × 10^{-15} (LIGO O3).
     */

    if (G_mu <= 0.0) return 0.0;

    /* Approximate flat spectrum in the LIGO band */
    double Omega_ref = 2.0e-10 * (G_mu / 1.0e-8);
    return Omega_ref;
}

/* ================================================================
 * L7 — Overlap Reduction Function for LIGO Hanford-Livingston
 *
 * gamma(f) encodes the geometric correlation between two detectors.
 *
 * For LIGO H-L:
 *   gamma(f) ~ -0.9 * sinc(2 pi f d / c) * [approx]
 *
 * where d ~ 3002 km is the separation.
 * ================================================================ */

double gw_overlap_reduction_hl(double f) {
    /* LIGO Hanford-Livingston baseline ~ 3002 km */
    double d = 3002000.0;
    double x = 2.0 * GW_PI * f * d / GW_C;
    if (fabs(x) < 1e-10) return 1.0;

    /* Approximate gamma(f) for LIGO H-L pair.
     * Actual gamma oscillates and decays with frequency.
     * Simplified model: gamma(f) ~ -cos(alpha) * sinc(x) * decay factor */
    double gamma = sin(x) / x;  /* sinc(x) */
    gamma *= -0.9;  /* approximate alignment factor */

    return gamma;
}

double gw_overlap_reduction(double f,
    const double dx[3],
    double lat1, double lon1, double lat2, double lon2,
    double arm1_az1, double arm1_az2,
    double arm2_az1, double arm2_az2) {
    /* Simplified: use the baseline length for phase factor,
     * and approximate the angular integral by averaging.
     *
     * Full calculation involves integrating over the sky
     * with antenna patterns — computationally expensive.
     *
     * Simplified: gamma(f) = sinc(2 pi f d / c) * (5/8pi) * integral...
     */
    (void)lat1; (void)lon1; (void)lat2; (void)lon2;
    (void)arm1_az1; (void)arm1_az2; (void)arm2_az1; (void)arm2_az2;

    double d = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (d < 1.0) return 1.0;

    double x = 2.0 * GW_PI * f * d / GW_C;
    if (fabs(x) < 1e-10) return 1.0;

    return sin(x) / x;
}

/* ================================================================
 * L7 — SGWB Detection SNR (Cross-Correlation)
 *
 * SNR^2 = 2 T_obs integral_0^inf df
 *         [ Omega_GW^2(f) * gamma^2(f) ] / [ f^6 S_n1(f) S_n2(f) ]
 *
 * For power-law Omega_GW(f) = Omega_alpha * (f/f_ref)^alpha
 * ================================================================ */

double gw_stochastic_snr(double T_obs, double Omega_alpha,
    double alpha, double f_min, double f_max) {
    /* Mid-frequency approximation: evaluate integrand at
     * a set of frequency points and sum.
     *
     * Integrand ~ Omega_GW^2(f) * gamma^2(f) / (f^6 * S_n1 * S_n2)
     */

    int nf = 100;
    double df = (f_max - f_min) / (nf - 1.0);
    double f_ref = 25.0;
    double sum = 0.0;

    for (int i = 0; i < nf; i++) {
        double f = f_min + i * df;
        if (f <= 0.0) continue;

        double Omega_f = Omega_alpha * pow(f / f_ref, alpha);
        double gamma_f = gw_overlap_reduction_hl(f);
        double Sn1 = gw_psd_aligo_design(f);
        double Sn2 = Sn1;  /* assume identical detectors */

        if (Sn1 <= 0.0) continue;

        double integrand = Omega_f * Omega_f * gamma_f * gamma_f
                         / (pow(f, 6.0) * Sn1 * Sn2);

        sum += integrand * df;
    }

    double snr_sq = 2.0 * T_obs * sum * (3.0 * 1000.0 * 1000.0 * 1000.0)
                    / (10.0 * GW_PI * GW_PI);
    /* The normalization factor (3 H_0^2 / 10 pi^2)^2 is included
     * in the definition of Omega_GW — simplified here. */

    return sqrt(snr_sq);
}
