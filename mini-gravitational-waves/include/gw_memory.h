/**
 * @file    gw_memory.h
 * @brief   Gravitational Wave Memory Effect
 *
 * Reference: Zel'dovich & Polnarev (1974), Christodoulou (1991),
 *            Favata (2010), Lasky et al. (2016)
 *
 * L2 — Linear and nonlinear GW memory
 * L4 — Memory formula: permanent strain offset after GW passage
 * L7 — Detectability with LIGO/Virgo, PTA
 */

#ifndef GW_MEMORY_H
#define GW_MEMORY_H

#include "gw_core.h"

/* ================================================================
 * L2 — Linear Memory
 * ================================================================ */

/**
 * Linear GW memory: permanent strain offset from unbound sources.
 *
 * Caused by change in second time derivative of quadrupole moment
 * between t = -inf and t = +inf.
 *
 *   Delta h_{ij}^{TT} = (2G / c^4 r) * Lambda_{ij,kl}
 *                     * [d I_{kl}/dt |_{-inf}^{+inf}]
 *
 * For a system that goes from unbound to unbound (hyperbolic).
 *
 * @param Dh_tt   Output memory tensor (TT gauge)
 * @param Idot_initial  dI_{ij}/dt at t = -inf [kg m^2/s]
 * @param Idot_final    dI_{ij}/dt at t = +inf [kg m^2/s]
 * @param r       Distance [m]
 * @param nx,ny,nz  Propagation direction
 */
void gw_linear_memory(GwTensor3 *Dh_tt,
    const GwTensor3 *Idot_initial, const GwTensor3 *Idot_final,
    double r, double nx, double ny, double nz);

/* ================================================================
 * L2 — Nonlinear (Christodoulou) Memory
 * ================================================================ */

/**
 * Nonlinear GW memory from GW energy flux radiated.
 *
 *   Delta h_{ij}^{nl}(theta,phi) =
 *     (4G / c^4 r) integral_{-inf}^{t} dt'
 *     integral dOmega' [ dE_GW/dt'/dOmega' * ... ]
 *
 * Simplified: memory amplitude ~ (G / c^4 r) * Delta E_GW / theta_jet
 *
 * For compact binary mergers, the nonlinear memory builds up
 * mainly around merger.
 *
 * @param plus_mem   Output +-polarization memory amplitude
 * @param cross_mem  Output x-polarization memory amplitude
 * @param Delta_E_GW Total energy radiated in GWs [J]
 * @param r          Distance [m]
 * @param iota       Inclination angle [rad]
 */
void gw_nonlinear_memory(double *plus_mem, double *cross_mem,
                         double Delta_E_GW, double r, double iota);

/* ================================================================
 * L4 — Memory from Compact Binary Coalescence
 * ================================================================ */

/**
 * Nonlinear memory for aligned-spin CBC.
 *
 * Uses Favata (2010) fitting formula:
 *
 *   h_mem(t) ~ (G/c^4 r) * (mu M / a(t)) * sin^2 iota * (17 + cos^2 iota)
 *
 * where a(t) is the orbital separation.
 *
 * @param h_plus_mem   Output memory strain (plus only for non-precessing)
 * @param M_total      Total mass [kg]
 * @param mu           Reduced mass [kg]
 * @param a            Orbital separation [m]
 * @param r            Distance [m]
 * @param iota         Inclination [rad]
 */
void gw_cbc_memory(double *h_plus_mem,
    double M_total, double mu, double a, double r, double iota);

/**
 * Late-time memory growth for CBC.
 * Memory grows like log(f) during inspiral.
 *
 * @param memory_out  Memory amplitude at frequency f
 * @param Mc          Chirp mass [kg]
 * @param r           Distance [m]
 * @param iota        Inclination [rad]
 * @param f           Current GW frequency [Hz]
 */
void gw_memory_growth(double *memory_out,
    double Mc, double r, double iota, double f);

/* ================================================================
 * L7 — Memory SNR Estimation
 * ================================================================ */

/**
 * Estimate SNR of the memory signal in aLIGO.
 *
 * Memory is DC/low-frequency; dominant SNR comes from stacking
 * multiple events or from low-frequency detectors (LISA, PTA).
 *
 * @param Mc       Chirp mass [kg]
 * @param D_L      Luminosity distance [m]
 * @param iota     Inclination [rad]
 * @return Approximate memory SNR in aLIGO (single event)
 */
double gw_memory_snr_aligo(double Mc, double D_L, double iota);

/**
 * Number of CBC events needed to detect memory at SNR=3,
 * stacking N events: SNR_N = sqrt(N) * SNR_1.
 */
double gw_memory_events_for_detection(double snr_single);

#endif /* GW_MEMORY_H */
