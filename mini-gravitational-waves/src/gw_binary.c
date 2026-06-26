/**
 * @file    gw_binary.c
 * @brief   Compact Binary Dynamics — Chirp, PN Evolution, Orbital Decay
 *
 * L1 — Chirp mass, mass ratio, reduced mass
 * L2 — Frequency evolution, time-to-coalescence
 * L4 — Kepler + quadrupole: binary inspiral under GW emission
 * L5 — Post-Newtonian (PN) phase approximants
 * L6 — GW150914, BNS, NSBH reference systems
 */

#include "gw_binary.h"
#include <math.h>

/* ================================================================
 * L1 — Binary Mass Parameters
 * ================================================================ */

double gw_chirp_mass(double m1, double m2) {
    double M = m1 + m2;
    if (M <= 0.0) return 0.0;
    double mu = m1 * m2 / M;
    return pow(mu, 0.6) * pow(M, 0.4);  /* = (m1*m2)^(3/5) / M^(1/5) */
}

double gw_symmetric_mass_ratio(double m1, double m2) {
    double M = m1 + m2;
    if (M <= 0.0) return 0.0;
    return m1 * m2 / (M * M);
}

double gw_reduced_mass(double m1, double m2) {
    double M = m1 + m2;
    if (M <= 0.0) return 0.0;
    return m1 * m2 / M;
}

void gw_binary_params_init(GwBinaryParams *p,
    double m1, double m2, double D_L, double iota,
    double f_ref, double phi_ref) {
    p->m1     = m1;
    p->m2     = m2;
    p->M_total = m1 + m2;
    p->chirp_mass = gw_chirp_mass(m1, m2);
    p->eta    = gw_symmetric_mass_ratio(m1, m2);
    p->D_L    = D_L;
    p->iota   = iota;
    p->f_ref  = f_ref;
    p->phi_ref = phi_ref;
    p->t_c    = 0.0;
    p->phi_c  = 0.0;
}

/* ================================================================
 * L2/L4 — Inspiral Frequency Evolution (Newtonian)
 *
 * Starting from:
 *   E_orb = - G M mu / (2a)          (Kepler energy)
 *   dE_orb/dt = - L_GW               (energy balance)
 *   L_GW = (32/5) G^4 mu^2 M^3 / (c^5 a^5)  (quadrupole)
 *   a = (G M)^{1/3} / (pi f)^{2/3}   (Kepler III)
 *
 * Chain rule gives:
 *   df/dt = (96/5) pi^{8/3} (G Mc / c^3)^{5/3} f^{11/3}
 * ================================================================ */

double gw_frequency_derivative(double Mc, double f) {
    if (f <= 0.0 || Mc <= 0.0) return 0.0;
    double GMc_over_c3 = GW_G * Mc / (GW_C * GW_C * GW_C);
    return (96.0 / 5.0) * GW_PI_8_3 * pow(GMc_over_c3, 5.0/3.0) * pow(f, 11.0/3.0);
}

double gw_time_to_coalescence(double Mc, double f) {
    if (f <= 0.0 || Mc <= 0.0) return 1e100;
    double GMc_over_c3 = GW_G * Mc / (GW_C * GW_C * GW_C);
    double f_pow = pow(f, -8.0/3.0);
    return (5.0 / 256.0) * pow(GW_PI, -8.0/3.0)
           * pow(GMc_over_c3, -5.0/3.0) * f_pow;
}

double gw_frequency_at_time(double Mc, double delta_t) {
    if (delta_t <= 0.0 || Mc <= 0.0) return 1e100;
    double GMc_over_c3 = GW_G * Mc / (GW_C * GW_C * GW_C);
    double fac = 256.0 / 5.0 * delta_t * pow(GMc_over_c3, 5.0/3.0);
    return (1.0 / GW_PI) * pow(fac, -3.0/8.0);
}

double gw_number_of_cycles(double Mc, double f_min, double f_max) {
    if (f_min >= f_max || Mc <= 0.0) return 0.0;
    double GMc_over_c3 = GW_G * Mc / (GW_C * GW_C * GW_C);
    return (1.0 / (32.0 * GW_PI_8_3))
           * pow(GMc_over_c3, -5.0/3.0)
           * (pow(f_min, -5.0/3.0) - pow(f_max, -5.0/3.0));
}

/* ================================================================
 * L5 — PN Phase Coefficients (up to 3.5PN)
 *
 * phi_k for the TaylorT2/F2 phase. These are functions of eta.
 *
 * phi_0 = 1
 * phi_1 = 0
 * phi_2 = 3715/756 + 55*eta/9
 * phi_3 = -16*pi
 * phi_4 = 15293365/508032 + 27145*eta/504 + 3085*eta^2/72
 * phi_5 = pi * (38645/756 - 65*eta/9) * [1 + log_terms]
 * phi_6 = 11583231236531/4694215680 - 640*pi^2/3 - 6848*gamma_E/21
 *         + (-15737765635/3048192 + 2255*pi^2/12)*eta
 *         + 76055*eta^2/1728 - 127825*eta^3/1296
 * phi_7 = pi * (77096675/254016 + 378515*eta/1512 - 74045*eta^2/756)
 *
 * Here we implement up to 3.5PN (phi_0..phi_7).
 * ================================================================ */

#define GW_EULER_GAMMA  0.5772156649015328606065120900824024310421

static void gw_pn_phase_coeffs(double coeffs[8], double eta) {
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;
    double pi   = GW_PI;

    coeffs[0] = 1.0;
    coeffs[1] = 0.0;
    coeffs[2] = 3715.0/756.0 + 55.0*eta/9.0;
    coeffs[3] = -16.0 * pi;
    coeffs[4] = 15293365.0/508032.0 + 27145.0*eta/504.0 + 3085.0*eta2/72.0;
    coeffs[5] = pi * (38645.0/756.0 - 65.0*eta/9.0);
    coeffs[6] = 11583231236531.0/4694215680.0
                - 640.0*pi*pi/3.0
                - 6848.0*GW_EULER_GAMMA/21.0
                + (-15737765635.0/3048192.0 + 2255.0*pi*pi/12.0)*eta
                + 76055.0*eta2/1728.0
                - 127825.0*eta3/1296.0;
    coeffs[7] = pi * (77096675.0/254016.0 + 378515.0*eta/1512.0
                      - 74045.0*eta2/756.0);
}

/* ================================================================
 * L5 — TaylorT2 Phase (Frequency Domain)
 *
 * Phi(f) = 2 pi f t_c - phi_c - pi/4
 *        + (3/128) (pi G M f / c^3)^{-5/3} sum_{k=0}^{7} phi_k x^{k/2}
 * ================================================================ */

double gw_taylor_t2_phase(double f, double M_total, double eta,
                          double t_c, double phi_c) {
    if (f <= 0.0 || M_total <= 0.0) return 0.0;

    double x = gw_pn_parameter(M_total, f);  /* x = v^2/c^2 */
    double x_half = sqrt(x);

    double coeffs[8];
    gw_pn_phase_coeffs(coeffs, eta);

    double sum = 0.0;
    double xk = 1.0;  /* x^{0/2} */
    for (int k = 0; k < 8; k++) {
        sum += coeffs[k] * xk;
        xk *= x_half;
    }

    /* Formula: Phi = 2 pi f t_c - phi_c - pi/4 + (3)/(128 eta) * x^{-5/2} * sum_{k} phi_k x^{k/2} */
    double phase = 2.0 * GW_PI * f * t_c - phi_c - GW_PI / 4.0
                   + (3.0 / (128.0 * eta)) * pow(x, -2.5) * sum;

    return phase;
}

/* ================================================================
 * L5 — TaylorT2 frequency-domain phase (TaylorF2)
 *
 * Psi(f) = 2 pi f t_c - phi_c - pi/4
 *        + (3/128 eta) * (pi G M f / c^3)^{-5/3}
 *        * sum_{k=0}^{7} psi_k x^{k/2}
 *
 * where psi_k = phi_k, same coefficients.
 * ================================================================ */

double gw_taylor_f2_phase(double f, double M_total, double eta,
                          double t_c, double phi_c) {
    /* TaylorF2 phase is identical to TaylorT2 in form.
     * The difference is only in how the phase is used:
     * T2: Phi(f) -> time-domain via inverse Fourier
     * F2: Psi(f) directly in frequency-domain waveform
     */
    return gw_taylor_t2_phase(f, M_total, eta, t_c, phi_c);
}

/* ================================================================
 * L5 — TaylorT4 RHS (Time Domain ODE)
 *
 * dx/dt = (64/5) nu x^5 * (1 + PN corrections) / M
 * dPhi/dt = x^{3/2} / M
 *
 * where x = (pi G M f / c^3)^{2/3}
 *       nu = eta
 *
 * We output df/dt and dPhi/dt for the user to integrate.
 * ================================================================ */

void gw_taylor_t4_rhs(double f, double M_total, double eta,
                      double *df_dt, double *dphi_dt) {
    if (f <= 0.0 || M_total <= 0.0) {
        *df_dt = 0.0;
        *dphi_dt = 0.0;
        return;
    }

    double x = gw_pn_parameter(M_total, f);

    /* PN corrections to energy flux (up to 3.5PN) */
    double x_half = sqrt(x);
    double x2 = x * x;
    double x3 = x2 * x;
    double x35 = x3 * x_half; /* x^{7/2} */

    double flux_corr = 1.0
        - (1247.0/336.0 + 35.0*eta/12.0) * x
        + 4.0 * GW_PI * x * x_half
        + (-44711.0/9072.0 + 9271.0*eta/504.0 + 65.0*eta*eta/18.0) * x2
        + (-8191.0/672.0 - 583.0*eta/24.0) * GW_PI * x2 * x_half
        + (6643739519.0/69854400.0 + 16.0*GW_PI*GW_PI/3.0
           - 1712.0*GW_EULER_GAMMA/105.0
           + (-134543.0/7776.0 + 41.0*GW_PI*GW_PI/48.0)*eta
           - 94403.0*eta*eta/3024.0 - 775.0*eta*eta*eta/324.0) * x3
        + (-16285.0/504.0 + 214745.0*eta/1728.0
           + 193385.0*eta*eta/3024.0) * GW_PI * x35;

    /* dx/dt from flux balance */
    double dx_dt = (64.0/5.0) * eta * x * x * x * x * x * flux_corr
                   * (GW_C * GW_C * GW_C) / (GW_G * M_total);

    /* Convert dx/dt to df/dt:
     * x = (pi G M f / c^3)^{2/3} => f = (c^3 / (pi G M)) x^{3/2}
     * df/dt = (3/2) (c^3 / (pi G M)) x^{1/2} dx/dt
     */
    double f_scale = GW_C * GW_C * GW_C / (GW_PI * GW_G * M_total);
    *df_dt = 1.5 * f_scale * x_half * dx_dt;
    *dphi_dt = 2.0 * GW_PI * f;  /* orbital phase = 2 pi f_orb, f = 2 f_orb */
}

/* ================================================================
 * L5 — Orbital Evolution ODE (Peters 1964)
 *
 * For eccentric binaries.
 * ================================================================ */

void gw_orbital_evolution_rhs(double dydt[3], const double y[3],
                              double M, double mu) {
    double a = y[0];
    double e = y[1];
    (void)y[2]; /* Phi stored for completeness */

    if (a <= 0.0 || e < 0.0 || e >= 1.0) {
        dydt[0] = 0.0;
        dydt[1] = 0.0;
        dydt[2] = 0.0;
        return;
    }

    double G3 = GW_G * GW_G * GW_G;
    double C5 = GW_C * GW_C * GW_C * GW_C * GW_C;
    double M2 = M * M;
    double muM2 = mu * M2;

    double e2 = e * e;
    double one_minus_e2 = 1.0 - e2;
    double f_e = (1.0 + 73.0*e2/24.0 + 37.0*e2*e2/96.0)
                 / pow(one_minus_e2, 3.5);
    double g_e = (1.0 + 121.0*e2/304.0)
                 / pow(one_minus_e2, 2.5);

    /* da/dt */
    dydt[0] = -(64.0/5.0) * (G3 / C5) * (muM2) / (a*a*a) * f_e;

    /* de/dt */
    dydt[1] = -(304.0/15.0) * (G3 / C5) * (muM2) / (a*a*a*a) * g_e * e;

    /* dPhi/dt = omega = sqrt(GM/a^3) */
    dydt[2] = sqrt(GW_G * M / (a*a*a));
}

/* ================================================================
 * L5 — Eccentricity Enhancement Functions (Peters & Mathews 1963)
 * ================================================================ */

double gw_eccentricity_enhancement_da(double e) {
    double e2 = e * e;
    double denom = pow(1.0 - e2, 3.5);
    return (1.0 + 73.0*e2/24.0 + 37.0*e2*e2/96.0) / denom;
}

double gw_eccentricity_enhancement_de(double e) {
    double e2 = e * e;
    double denom = pow(1.0 - e2, 2.5);
    return (1.0 + 121.0*e2/304.0) / denom;
}

double gw_peak_frequency_eccentric(double M, double a, double e) {
    double omega_orb = sqrt(GW_G * M / (a*a*a));
    double f_orb = omega_orb / (2.0 * GW_PI);
    double f_peak = f_orb * sqrt(1.0 + e) / pow(1.0 - e, 1.5);
    /* Peak GW frequency: 2 * f_orb * (1+e)^{1/2} / (1-e)^{3/2} for quadrupole */
    return 2.0 * f_peak;
}

/* ================================================================
 * L6 — Reference Systems
 * ================================================================ */

void gw_params_gw150914(GwBinaryParams *p) {
    double m1 = 36.0 * GW_MSUN;
    double m2 = 29.0 * GW_MSUN;
    double D_L = 410.0 * GW_MPC;
    gw_binary_params_init(p, m1, m2, D_L, 0.0, 20.0, 0.0);
    p->t_c   = 0.0;
    p->phi_c = 0.0;
}

void gw_params_bns(GwBinaryParams *p) {
    double m1 = 1.4 * GW_MSUN;
    double m2 = 1.4 * GW_MSUN;
    double D_L = 40.0 * GW_MPC;
    gw_binary_params_init(p, m1, m2, D_L, 0.0, 10.0, 0.0);
}

void gw_params_nsbh(GwBinaryParams *p) {
    double m1 = 10.0 * GW_MSUN;
    double m2 = 1.4 * GW_MSUN;
    double D_L = 100.0 * GW_MPC;
    gw_binary_params_init(p, m1, m2, D_L, 0.0, 15.0, 0.0);
}
