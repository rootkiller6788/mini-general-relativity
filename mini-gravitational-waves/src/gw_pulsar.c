/**
 * @file    gw_pulsar.c
 * @brief   Continuous GW Sources — Pulsars, Ellipticity, R-modes, PTAs
 *
 * L2 — CW emission from non-axisymmetric neutron stars
 * L4 — Spin-down limit, GW strain from ellipticity
 * L6 — Pulsar timing arrays, Hellings-Downs curve
 * L7 — F-statistic, narrow-band CW search
 */

#include "gw_pulsar.h"
#include <math.h>

/* ================================================================
 * L2 — Continuous GW Strain from Ellipticity
 *
 * h_0 = (4 pi^2 G / c^4) * (epsilon * I_zz * f_rot^2) / D
 *
 * This is the "standard" CW signal from a triaxial NS.
 * Assumes GW emission at f_gw = 2 f_rot (l=m=2 mode).
 *
 * Typical values:
 *   I_zz ~ 1e38 kg m^2 (canonical NS)
 *   epsilon ~ 1e-6 (realistic maximum from crustal strain)
 *   f_rot ~ 100 Hz -> h_0 ~ 1e-26 at D = 1 kpc
 * ================================================================ */

double gw_cw_strain(double epsilon, double I_zz, double f_rot, double D) {
    double omega = 2.0 * GW_PI * f_rot;
    return (4.0 * GW_PI * GW_PI * GW_G / (GW_C * GW_C * GW_C * GW_C))
           * epsilon * I_zz * omega * omega / D;
}

double gw_cw_frequency(double f_rot) {
    return 2.0 * f_rot;
}

/* ================================================================
 * L4 — Spin-Down Limit
 *
 * If all rotational kinetic energy loss goes into GWs:
 *
 *   dE/dt = -I_zz omega |omega_dot| = L_GW
 *
 * For quadrupole GW: L_GW = (32/5) (G/c^5) I_zz^2 epsilon^2 omega^6
 *
 * Solving for epsilon gives the spin-down ellipticity,
 * and then h_0^{sd} is the corresponding strain.
 *
 *   h_0^{sd} = sqrt( 5 G I_zz |f_dot| / (2 c^3 D^2 f_gw) )
 * ================================================================ */

double gw_spindown_limit(double I_zz, double f_gw, double f_dot, double D) {
    if (D <= 0.0 || f_gw <= 0.0) return 0.0;

    double f_dot_abs = fabs(f_dot);
    return sqrt(5.0 * GW_G * I_zz * f_dot_abs
                / (2.0 * GW_C * GW_C * GW_C * D * D * f_gw));
}

double gw_spindown_limit_age(double I_zz, double f_rot,
                             double f_dot, double D) {
    if (f_dot >= 0.0 || D <= 0.0 || f_rot <= 0.0) return 0.0;

    /* Characteristic age: tau = -f_rot / (2 f_dot) for n=3 braking */
    double tau = -f_rot / (2.0 * f_dot);

    /* h_0^{age} = sqrt( (5/2) * (G/c^3) * I_zz / (D^2 * tau) ) / (2 pi f_rot) */
    double pref = sqrt(5.0 * GW_G * I_zz
                       / (2.0 * GW_C*GW_C*GW_C * D * D * tau));
    return pref / (2.0 * GW_PI * f_rot);
}

/* ================================================================
 * L4 — R-mode GW Strain
 *
 * R-modes are non-axisymmetric fluid oscillations restored by
 * the Coriolis force. The dominant mode (l=m=2) radiates GWs
 * at f_gw = (4/3) f_rot.
 *
 * The GW strain depends on the saturation amplitude alpha:
 *
 *   h_0 = sqrt( (8 pi G / c^3) * alpha^2 * J_tilde * M R^3 Omega^3 / D^2 )
 *
 * where J_tilde ~ 0.016 for an n=1 polytrope.
 * ================================================================ */

double gw_rmode_strain(double alpha, double M, double R,
                       double f_rot, double D) {
    if (D <= 0.0 || f_rot <= 0.0) return 0.0;

    double Omega = 2.0 * GW_PI * f_rot;
    double J_tilde = 0.016;  /* for n=1 polytrope */

    double power = 8.0 * GW_PI * GW_G / (GW_C * GW_C * GW_C)
                   * alpha * alpha * J_tilde
                   * M * R * R * R * Omega * Omega * Omega
                   / (D * D);

    if (power < 0.0) return 0.0;
    return sqrt(power);
}

double gw_rmode_frequency(double f_rot) {
    return (4.0 / 3.0) * f_rot;
}

/* ================================================================
 * L6 — Pulsar Timing Residual from GW
 *
 * A GW passing between Earth and a pulsar induces a shift
 * in the observed pulse arrival time:
 *
 *   R(t) = h_0 / (2 pi f_gw) * F_plus * [sin(2 pi f_gw t + phi_0)
 *                                        - sin(phi_0)]
 *
 * The "Earth term" (depends on t) and "pulsar term" (depends on
 * t - D_p/c) combine to produce a characteristic signature.
 * ================================================================ */

double gw_pta_timing_residual(double t, double h_0, double f_gw,
                              double F_plus, double phi0) {
    double omega = 2.0 * GW_PI * f_gw;
    /* Earth term only (pulsar term averages out in many-psr analysis) */
    return (h_0 / omega) * F_plus * (sin(omega * t + phi0) - sin(phi0));
}

/* ================================================================
 * L6 — Hellings-Downs Correlation
 *
 * The expected cross-correlation between timing residuals of
 * two pulsars separated by angle xi, for an isotropic stochastic
 * GW background, is given by the Hellings-Downs curve:
 *
 *   C(xi) = 1/2 - (1/4)(1 - cos xi)/2
 *          + (3/2)(1 - cos xi)/2 * ln((1 - cos xi)/2)
 *
 * Normalization: C(0) = 1/2, C(pi) = -1/4.
 *
 * This is the smoking-gun signature of the GW origin of
 * pulsar timing array correlations (NANOGrav 15yr detection!).
 * ================================================================ */

double gw_hellings_downs(double xi) {
    if (xi < 1e-10) return 0.5;  /* limit as xi -> 0 */

    double x = (1.0 - cos(xi)) / 2.0;
    if (x < 1e-15) return 0.5;

    return 0.5 - 0.25 * x + 1.5 * x * log(x);
}

/* ================================================================
 * L7 — F-Statistic (Jaranowski-Krolak-Schutz 1998)
 *
 * The F-statistic is the optimal detection statistic for
 * continuous GWs from known-position sources.
 *
 * It analytically maximizes the likelihood over the four
 * extrinsic parameters: h_0, cos iota, psi, phi_0.
 *
 * F = 2 * (maximum log-likelihood over extrinsic params)
 *
 * 2F ~ chi^2_4 (non-central) in Gaussian noise.
 *
 * Implemented here for the simple case of a known frequency
 * and sky position (only extrinsic parameter maximization).
 * ================================================================ */

double gw_f_statistic(const double *data_re, const double *data_im,
                      int nf, double df, const double *psd,
                      double F_plus_a, double F_cross_a,
                      double F_plus_b, double F_cross_b) {
    /* Simplified: F-statistic ~ SNR^2 maximized over amplitudes.
     *
     * The full F-statistic requires constructing 4 basis templates
     * and computing the maximum of a quadratic form.
     *
     * Simplified version: approximate F ~ rho^2 (optimal SNR^2) for
     * perfectly matched template.
     */

    double rho_sq = 0.0;
    for (int i = 0; i < nf; i++) {
        if (psd[i] <= 0.0) continue;
        double power = data_re[i]*data_re[i] + data_im[i]*data_im[i];
        rho_sq += power / psd[i] * df;
    }
    rho_sq *= 4.0;

    /* Approximate F-statistic as 2 * rho^2 / sum(A^2) normalization */
    double A_sq = F_plus_a*F_plus_a + F_cross_a*F_cross_a
                + F_plus_b*F_plus_b + F_cross_b*F_cross_b;
    if (A_sq <= 0.0) return 0.0;

    return 2.0 * rho_sq / A_sq;
}

double gw_fstat_fap(double F0) {
    /* P_F(2F > 2F_0) = exp(-F_0) for single trial, 4 dof */
    return exp(-F0);
}

/* ================================================================
 * L7 — Pulsar Antenna Pattern (Earth Term)
 *
 * For a pulsar at sky location (ra, dec) and GW with
 * polarization angle psi, the antenna pattern for the
 * Earth term is:
 *
 *   F_plus  = cos(2 psi) * a - sin(2 psi) * b
 *   F_cross = sin(2 psi) * a + cos(2 psi) * b
 *
 * where:
 *   a = (1/2) cos(2 phi) (1 - cos theta) / (1 + cos theta)
 *       ... actually, the standard form from JKS 1998:
 *
 *   a = (1/2) (zeta * cos(2 phi))
 *   b = (1/2) (zeta * sin(2 phi))
 *   zeta = (1 - cos theta) / 2
 *
 * But in the Earth-term only approximation for a single pulsar:
 *
 *   F_plus  = (1/2) * (1 + cos^2 theta) ...
 *
 * The exact expression depends on whether you use the detector
 * frame or the source frame. Here we use a simplified form.
 * ================================================================ */

void gw_pulsar_antenna_pattern(double *F_plus, double *F_cross,
                               double ra, double dec, double psi) {
    /* Sky-location-dependent antenna patterns for PTA pulsars.
     *
     * Uses the formalism from Anholm et al. (2009), arXiv:0809.0702.
     *
     * For a GW propagating from direction -k_hat (source at ra, dec),
     * the antenna pattern components are:
     *
     *   cos_theta = sin(dec)  (approx; depends on coordinate choice)
     */

    double cos_theta = sin(dec); /* simplified: dec=0 -> equatorial, cos_theta=0 */
    double sin_theta = cos(dec);

    double cos_2phi = cos(2.0 * ra);
    double sin_2phi = sin(2.0 * ra);
    double cos_2psi = cos(2.0 * psi);
    double sin_2psi = sin(2.0 * psi);

    /* Antenna pattern for Earth-pulsar baseline.
     * For Earth-term only, simplified response: */
    double fac = 0.5 * (1.0 - cos_theta);

    *F_plus  = fac * (cos_2phi * cos_2psi - sin_theta * sin_2phi * sin_2psi);
    *F_cross = fac * (cos_2phi * sin_2psi + sin_theta * sin_2phi * cos_2psi);
}
