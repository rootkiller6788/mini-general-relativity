/**
 * @file    gw_stochastic.h
 * @brief   Stochastic Gravitational Wave Background (SGWB)
 *
 * Reference: Maggiore Vol.2 Ch.23, Allen & Romano (1999),
 *            NANOGrav 15-year results, LIGO-Virgo-KAGRA Stoch. analysis
 *
 * L2 — Stochastic background as superposition of unresolved sources
 * L4 — Omega_GW(f): fractional energy density in GWs
 * L6 — Astrophysical background (CBC), cosmological background (inflation)
 * L7 — Cross-correlation detection, overlap reduction function
 */

#ifndef GW_STOCHASTIC_H
#define GW_STOCHASTIC_H

#include "gw_core.h"

/* ================================================================
 * L2 — SGWB Energy Density
 * ================================================================ */

/**
 * Omega_GW(f): GW energy density per logarithmic frequency,
 * normalized by critical density.
 *
 *   Omega_GW(f) = (1 / rho_c) * d(rho_GW) / d(ln f)
 *
 * where rho_c = 3 H_0^2 c^2 / (8 pi G).
 *
 * @param f      Frequency [Hz]
 * @param params Stochastic background parameters
 * @return Omega_GW(f)
 */
double gw_omega_gw(double f, const GwStochasticParams *params);

/**
 * GW energy density rho_GW from Omega_GW.
 *
 * @param Omega  Omega_GW value
 * @param f_low  Lower bound [Hz]
 * @param f_high Upper bound [Hz]
 * @param H0     Hubble constant [km/s/Mpc]
 * @return rho_GW [J/m^3]
 */
double gw_energy_density_from_omega(double Omega,
    double f_low, double f_high, double H0);

/**
 * Critical density of the Universe.
 *
 *   rho_c = 3 H_0^2 / (8 pi G)
 *
 * @param H0  Hubble constant [km/s/Mpc]
 * @return rho_c [kg/m^3]
 */
double gw_critical_density(double H0);

/**
 * Convert H0 from [km/s/Mpc] to [1/s].
 */
double gw_H0_to_per_sec(double H0_km_s_Mpc);

/* ================================================================
 * L6 — Astrophysical Background Models
 * ================================================================ */

/**
 * Omega_GW from compact binary coalescence (CBC) background.
 *
 * Power-law model: Omega_CBC(f) = Omega_ref * (f / f_ref)^{2/3}
 *
 * (The 2/3 spectral index comes from inspiral-only GW spectrum
 *  integrated over cosmic history.)
 *
 * @param f       Frequency [Hz]
 * @param Omega_ref Reference amplitude at f_ref = 25 Hz
 * @return Omega_CBC(f)
 */
double gw_omega_cbc_background(double f, double Omega_ref);

/**
 * Integrated Omega_GW from CBC population model.
 *
 * Sum over all mergers in cosmic history, folding in star formation
 * rate, delay time distribution, and merger rate evolution.
 *
 * Simplified: Omega(f) = A_ref * (f / 25Hz)^{2/3}
 *
 * @param A_ref   Amplitude at 25 Hz (e.g. from LVK O3: A ~ 5e-9)
 * @param f       Frequency [Hz]
 * @return Omega_GW
 */
double gw_omega_cbc_integrated(double f, double A_ref);

/* ================================================================
 * L6 — Cosmological Background Models
 * ================================================================ */

/**
 * Omega_GW from slow-roll inflation (single-field).
 *
 *   Omega_GW(f) ~ (r * A_s / ... ) * (f / f_pivot)^{n_t}
 *
 * where r = tensor-to-scalar ratio, n_t ~ -r/8 (consistency relation).
 *
 * @param f       Frequency [Hz]
 * @param r       Tensor-to-scalar ratio
 * @param A_s     Scalar amplitude (~2.1e-9 at k_pivot = 0.05 Mpc^{-1})
 * @return Omega_GW from inflationary tensors
 */
double gw_omega_inflation(double f, double r, double A_s);

/**
 * Omega_GW from first-order phase transitions.
 *
 * Omega_pt(f) ~ Omega_peak * S(f/f_peak)
 *
 * where S is a broken power-law spectral shape:
 *
 *   S(x) = (a+b) / (b x^{-a/c} + a x^{b/c})^c
 *
 * @param f          Frequency [Hz]
 * @param Omega_peak Peak amplitude
 * @param f_peak     Peak frequency [Hz]
 * @param a          Low-frequency slope
 * @param b          High-frequency slope
 * @param c          Transition width
 * @return Omega_GW from phase transition
 */
double gw_omega_phase_transition(double f,
    double Omega_peak, double f_peak,
    double a, double b, double c);

/**
 * Omega_GW from cosmic strings (Nambu-Goto).
 *
 * Omega_cs(f) ~ G mu * sum_n C_n(f) / n^{q-1}
 *
 * Simplified power law over LIGO band:
 *   Omega(f) ~ Omega_cs_ref * (f / f_ref)^{0}
 *
 * (Nearly flat spectrum from loops in radiation era.)
 *
 * @param f       Frequency [Hz]
 * @param G_mu    Dimensionless string tension
 * @return Omega_GW from cosmic strings
 */
double gw_omega_cosmic_strings(double f, double G_mu);

/* ================================================================
 * L7 — Overlap Reduction Function
 * ================================================================ */

/**
 * Overlap reduction function gamma(f) for two GW detectors.
 *
 * Normalized: gamma(f) = (5/8pi) * integral dOmega
 *   sum_A F_1^A(n) F_2^A(n) exp(2 pi i f n·Delta_x / c)
 *
 * where F^A are antenna patterns, Delta_x is detector separation.
 *
 * For co-located and co-aligned detectors: gamma(f=0) = 1.
 *
 * This implementation computes gamma for the LIGO Hanford-Livingston
 * pair using the known baseline geometry.
 *
 * @param f  Frequency [Hz]
 * @return gamma(f) (real part)
 */
double gw_overlap_reduction_hl(double f);

/**
 * General overlap reduction function for two detectors with given
 * positions and orientations.
 *
 * @param f     Frequency [Hz]
 * @param dx    Baseline vector [m], length 3
 * @param lat1,lon1  Detector 1 location [rad]
 * @param lat2,lon2  Detector 2 location [rad]
 * @param arm1_az1  Detector 1 arm 1 azimuth [rad]
 * @param arm1_az2  Detector 1 arm 2 azimuth [rad]
 * @param arm2_az1  Detector 2 arm 1 azimuth [rad]
 * @param arm2_az2  Detector 2 arm 2 azimuth [rad]
 * @return gamma(f) (real part)
 */
double gw_overlap_reduction(double f,
    const double dx[3],
    double lat1, double lon1, double lat2, double lon2,
    double arm1_az1, double arm1_az2,
    double arm2_az1, double arm2_az2);

/**
 * SNR for stochastic background cross-correlation search.
 *
 *   SNR^2 = 2 T_obs integral_0^inf df
 *           [ Omega_GW^2(f) * gamma^2(f) ] / [ f^6 S_n1(f) S_n2(f) ]
 *
 * Simplified: mid-frequency approximation with band integration.
 *
 * @param T_obs       Observation time [s]
 * @param Omega_alpha Power-law amplitude
 * @param alpha       Power-law index
 * @param f_min, f_max Integration range [Hz]
 * @return Expected SNR
 */
double gw_stochastic_snr(double T_obs, double Omega_alpha,
    double alpha, double f_min, double f_max);

#endif /* GW_STOCHASTIC_H */
