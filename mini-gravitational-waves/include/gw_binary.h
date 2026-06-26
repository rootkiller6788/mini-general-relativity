/**
 * @file    gw_binary.h
 * @brief   Compact Binary Systems — Chirp, Inspiral, Merger
 *
 * Reference: Maggiore Vol.1 Ch.4, Blanchet (2014) Living Rev.Rel.,
 *            MIT 8.962 — Binary coalescence
 *
 * L1 — Chirp mass, symmetric mass ratio, PN parameters
 * L2 — Binary inspiral under GW emission
 * L4 — Kepler's third law + GW energy balance → frequency evolution
 * L5 — PN phase evolution (TaylorT2, TaylorT4 approximants)
 * L6 — GW150914-like binary, NS-NS mergers
 */

#ifndef GW_BINARY_H
#define GW_BINARY_H

#include "gw_core.h"

/* ================================================================
 * L1 — Binary Mass Parameters
 * ================================================================ */

/**
 * Chirp mass: M_c = (m1 m2)^{3/5} / (m1 + m2)^{1/5}
 *
 * Dominant parameter controlling inspiral GW phase evolution.
 */
double gw_chirp_mass(double m1, double m2);

/**
 * Symmetric mass ratio: eta = m1 * m2 / (m1 + m2)^2
 *
 * Dimensionless, 0 < eta <= 1/4 (max for equal masses).
 */
double gw_symmetric_mass_ratio(double m1, double m2);

/**
 * Reduced mass: mu = m1 * m2 / (m1 + m2)
 */
double gw_reduced_mass(double m1, double m2);

/**
 * Initialize GwBinaryParams from component masses and distance.
 */
void gw_binary_params_init(GwBinaryParams *p,
    double m1, double m2, double D_L, double iota,
    double f_ref, double phi_ref);

/* ================================================================
 * L2/L4 — Inspiral Frequency Evolution (Newtonian)
 * ================================================================ */

/**
 * GW frequency evolution df/dt (Newtonian quadrupole):
 *
 *   df/dt = (96/5) pi^{8/3} (G M_c / c^3)^{5/3} f^{11/3}
 *
 * Derived from: dE_orb/dt = -L_GW, with E_orb = -G M mu / (2a)
 * and Kepler: f = (1/pi) sqrt(GM/a^3).
 *
 * @param Mc   Chirp mass [kg]
 * @param f    GW frequency [Hz]
 * @return df/dt [Hz/s]
 */
double gw_frequency_derivative(double Mc, double f);

/**
 * Time to coalescence from GW frequency f (Newtonian):
 *
 *   tau(f) = (5/256) (pi f)^{-8/3} (G M_c / c^3)^{-5/3}
 *
 * @param Mc  Chirp mass [kg]
 * @param f   Current GW frequency [Hz]
 * @return Time to coalescence [s]
 */
double gw_time_to_coalescence(double Mc, double f);

/**
 * GW frequency at time t before coalescence:
 *
 *   f(t) = (1/pi) [ (5/256) (G M_c / c^3)^{-5/3} / (t_c - t) ]^{3/8}
 *
 * @param Mc     Chirp mass [kg]
 * @param delta_t Time before coalescence t_c - t [s]
 * @return GW frequency [Hz]
 */
double gw_frequency_at_time(double Mc, double delta_t);

/**
 * Number of GW cycles in a frequency band (Newtonian):
 *
 *   N_cyc = (1 / (32 pi^{8/3})) (G M_c / c^3)^{-5/3}
 *           * (f_min^{-5/3} - f_max^{-5/3})
 */
double gw_number_of_cycles(double Mc, double f_min, double f_max);

/* ================================================================
 * L5 — Post-Newtonian Phase Evolution
 * ================================================================ */

/**
 * TaylorT2 approximant: phase Phi(f) as PN expansion in frequency domain.
 *
 * Phi(f) = 2 pi f t_c - Phi_c - pi/4
 *        + (3/128) (pi G M f / c^3)^{-5/3} * sum_{k=0}^{7} phi_k x^{k/2}
 *
 * where x = (pi G M f / c^3)^{2/3} and phi_k are PN coefficients.
 *
 * @param f       GW frequency [Hz]
 * @param M_total Total mass [kg]
 * @param eta     Symmetric mass ratio
 * @param t_c     Coalescence time [s]
 * @param phi_c   Coalescence phase [rad]
 * @return Orbital phase Phi(f) [rad]
 */
double gw_taylor_t2_phase(double f, double M_total, double eta,
                          double t_c, double phi_c);

/**
 * TaylorT4 approximant: dPhi/dt and df/dt combined ODE right-hand side.
 *
 * Evolves (f, Phi) as:
 *   df/dt  = (96/5) pi^{8/3} Mc^{5/3} f^{11/3} * [1 + PN corrections]
 *   dPhi/dt = 2 pi f
 *
 * @param f       Current GW frequency [Hz]
 * @param M_total Total mass [kg]
 * @param eta     Symmetric mass ratio
 * @param df_dt   Output: df/dt [Hz/s]
 * @param dphi_dt Output: dPhi/dt [rad/s]
 */
void gw_taylor_t4_rhs(double f, double M_total, double eta,
                      double *df_dt, double *dphi_dt);

/**
 * TaylorF2 approximant: frequency-domain phase only.
 * Suitable for matched filtering.
 *
 * @param f       GW frequency [Hz]
 * @param M_total Total mass [kg]
 * @param eta     Symmetric mass ratio
 * @param t_c     Coalescence time [s]
 * @param phi_c   Coalescence phase [rad]
 * @return Frequency-domain phase Psi(f) [rad]
 */
double gw_taylor_f2_phase(double f, double M_total, double eta,
                          double t_c, double phi_c);

/* ================================================================
 * L5 — Orbital Evolution ODE System
 * ================================================================ */

/**
 * Right-hand side of the 2.5PN orbital evolution equations.
 *
 * Evolves state vector y = [a, e, Phi]:
 *   da/dt = - (64/5) (G^3/c^5) (mu M^2) / a^3 * f_e(e)
 *   de/dt = - (304/15) (G^3/c^5) (mu M^2) / a^4 * g_e(e)
 *   dPhi/dt = sqrt(GM/a^3)
 *
 * where f_e, g_e are enhancement functions for eccentricity.
 *
 * @param dydt   Output derivatives [da/dt, de/dt, dPhi/dt]
 * @param y      State [a (m), e, Phi (rad)]
 * @param M      Total mass [kg]
 * @param mu     Reduced mass [kg]
 */
void gw_orbital_evolution_rhs(double dydt[3], const double y[3],
                              double M, double mu);

/* ================================================================
 * L5 — Eccentricity Enhancement Functions
 * ================================================================ */

/**
 * Peters & Mathews (1963) eccentricity enhancement factor f(e)
 * for da/dt of an eccentric binary under GW emission.
 *
 * f(e) = (1 + 73e^2/24 + 37e^4/96) / (1-e^2)^{7/2}
 */
double gw_eccentricity_enhancement_da(double e);

/**
 * Eccentricity enhancement factor g(e) for de/dt.
 *
 * g(e) = (1 + 121e^2/304) / (1-e^2)^{5/2}
 */
double gw_eccentricity_enhancement_de(double e);

/**
 * Peak GW frequency for eccentric binary:
 *   f_peak = sqrt(GM/a^3) / pi * (1+e)^{1/2} / (1-e)^{3/2}
 */
double gw_peak_frequency_eccentric(double M, double a, double e);

/* ================================================================
 * L6 — GW150914 Reference Values
 * ================================================================ */

/**
 * Fill parameters for a GW150914-like system.
 *
 * Source-frame masses: m1 ~ 36 Msun, m2 ~ 29 Msun
 * Luminosity distance: ~410 Mpc
 */
void gw_params_gw150914(GwBinaryParams *p);

/**
 * Fill parameters for a typical BNS (binary neutron star) system.
 * m1 = m2 = 1.4 Msun, D_L = 40 Mpc (Virgo cluster distance)
 */
void gw_params_bns(GwBinaryParams *p);

/**
 * Fill parameters for a typical NS-BH system.
 * m1 = 10 Msun, m2 = 1.4 Msun, D_L = 100 Mpc
 */
void gw_params_nsbh(GwBinaryParams *p);

#endif /* GW_BINARY_H */
