/**
 * @file    gw_memory.c
 * @brief   GW Memory Effect — Linear, Nonlinear, CBC Memory
 *
 * L2 — Linear memory from unbound systems
 * L4 — Nonlinear memory from GW energy flux
 * L7 — Memory detectability, stacking
 */

#include "gw_memory.h"
#include "gw_quadrupole.h"
#include <math.h>

/* ================================================================
 * L2 — Linear GW Memory
 *
 * Delta h_{ij}^{TT} = (2G / c^4 r) * Lambda_{ij,kl}
 *                   * [dI_{kl}/dt|_{t=-inf}^{t=+inf}]
 *
 * Arises when asymptotic velocities differ (hyperbolic encounters,
 * supernova kicks, etc.).
 * ================================================================ */

void gw_linear_memory(GwTensor3 *Dh_tt,
    const GwTensor3 *Idot_initial, const GwTensor3 *Idot_final,
    double r, double nx, double ny, double nz) {
    /* Delta I_dot */
    GwTensor3 Delta_Idot;
    gw_tensor_sub(&Delta_Idot, Idot_final, Idot_initial);

    /* Apply quadrupole formula scaling */
    GwTensor3 scaled;
    gw_quadrupole_strain(&scaled, &Delta_Idot, r);

    /* TT-project along n-hat */
    gw_tensor_tt_project(Dh_tt, &scaled, nx, ny, nz);
}

/* ================================================================
 * L4 — Nonlinear (Christodoulou) Memory
 *
 * Caused by GWs themselves sourcing additional GWs.
 * Memory amplitude proportional to total radiated energy:
 *
 *   h_mem ~ (G/c^4 r) * Delta_E_GW * sin^2 iota * (17 + cos^2 iota)
 *
 * For a CBC, Delta_E_GW ~ eta * M_total * c^2 * (radiated fraction).
 * ================================================================ */

void gw_nonlinear_memory(double *plus_mem, double *cross_mem,
                         double Delta_E_GW, double r, double iota) {
    double sin_iota = sin(iota);
    double cos_iota = cos(iota);

    /* Favata (2010) memory formula for CBC */
    double pref = GW_G / (GW_C * GW_C * GW_C * GW_C * r) * Delta_E_GW;

    /* Angular dependence:
     * h_mem^+ ~ sin^2 iota * (17 + cos^2 iota) / (something) */
    double ang_plus = sin_iota * sin_iota * (17.0 + cos_iota*cos_iota) / 16.0;

    *plus_mem  = pref * ang_plus;
    *cross_mem = 0.0;  /* Non-precessing CBC: memory is purely +-polarized */
}

/* ================================================================
 * L4 — CBC Memory Amplitude (Orbital Parameterization)
 * ================================================================ */

void gw_cbc_memory(double *h_plus_mem,
    double M_total, double mu, double a, double r, double iota) {
    /* Memory builds up roughly like h_mem ~ (G/c^4 r) * (mu M / a) * ...
     * This is a quasi-static memory estimate valid at a given separation. */

    double sin_iota = sin(iota);
    double cos_iota = cos(iota);

    double pref = GW_G * GW_G * M_total * mu
                / (GW_C * GW_C * GW_C * GW_C * r * a);

    *h_plus_mem = pref * sin_iota*sin_iota
                  * (17.0 + cos_iota*cos_iota) / 8.0;
}

/* ================================================================
 * L7 — Memory Growth During Inspiral
 *
 * During the inspiral, the memory grows roughly as:
 *
 *   h_mem(f) ~ h_mem_0 * log(f / f_0)
 *
 * This slow growth makes single-event detection challenging
 * in ground-based detectors (memory is DC/low-frequency),
 * but stacking many events enables detection.
 * ================================================================ */

void gw_memory_growth(double *memory_out,
    double Mc, double r, double iota, double f) {
    if (f <= 0.0) {
        *memory_out = 0.0;
        return;
    }

    /* Memory amplitude scales ~ (G Mc / c^2 r) * (G Mc f / c^3)^{2/3} */
    double GMc_c2 = GW_G * Mc / (GW_C * GW_C);
    double v2 = pow(GW_PI * GW_G * Mc * f / (GW_C*GW_C*GW_C), 2.0/3.0);
    double sin_iota = sin(iota);

    *memory_out = (GMc_c2 / r) * v2 * sin_iota*sin_iota * 0.2;
}

/* ================================================================
 * L7 — Memory SNR in aLIGO
 * ================================================================ */

double gw_memory_snr_aligo(double Mc, double D_L, double iota) {
    /* Memory is a DC offset in h(t). Ground-based detectors
     * are AC-coupled (seismic wall at ~10 Hz), so they are
     * insensitive to true DC memory.
     *
     * However, the "turn-on" of memory during the merger
     * produces power at higher frequencies.
     *
     * Approximate: SNR_mem ~ 0.05 * SNR_CBC for aLIGO
     * (Lasky et al. 2016, ApJ 825)
     */

    /* Rough estimate using GW150914-like parameters as baseline:
     * GW150914: Mc ~ 28 Msun, D_L ~ 410 Mpc, SNR ~ 24
     * Memory SNR ~ 0.5 for that event */

    double Mc_ref = 28.0 * GW_MSUN;
    double D_ref  = 410.0 * GW_MPC;
    double SNR_ref = 24.0;   /* CBC SNR for GW150914-like */
    double mem_SNR_ref = 0.5; /* ~2% of CBC SNR */

    /* Scale by Mc (roughly linear in amplitude), D_L (1/D_L) */
    double snr_cbc_est = SNR_ref * (Mc / Mc_ref) * (D_ref / D_L);
    double snr_mem_est = snr_cbc_est * 0.02;  /* ~2% of CBC SNR */

    return snr_mem_est;
}

double gw_memory_events_for_detection(double snr_single) {
    if (snr_single <= 0.0) return 1e100;
    double snr_target = 3.0;
    double N = (snr_target / snr_single) * (snr_target / snr_single);
    return N;
}
