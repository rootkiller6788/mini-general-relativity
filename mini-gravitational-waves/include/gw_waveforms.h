/**
 * @file    gw_waveforms.h
 * @brief   GW Waveform Models — Inspiral, IMR, Ringdown
 *
 * Reference: Maggiore Vol.1, Blanchet (2014), Khan et al. (2016) IMRPhenom,
 *            Bohé et al. (2017) SEOBNRv4
 *
 * L5 — Waveform generation (time + frequency domain)
 * L6 — IMR waveforms, ringdown modelling
 * L7 — Data analysis waveform templates
 */

#ifndef GW_WAVEFORMS_H
#define GW_WAVEFORMS_H

#include "gw_core.h"

/* ================================================================
 * L5 — Time-Domain Inspiral Waveform (Newtonian + PN)
 * ================================================================ */

/**
 * Generate time-domain inspiral-only waveform h(t).
 *
 * Uses TaylorT4 to evolve (f, Phi), then computes h_{+,x}(t).
 * Newtonian quadrupole amplitude, PN phase.
 *
 * @param s         Output strain series (must be pre-allocated)
 * @param params    Binary parameters
 * @param f_start   Starting GW frequency [Hz]
 * @param f_end     Ending GW frequency [Hz] (e.g. f_ISCO)
 * @param dt        Desired time step [s]
 * @param max_samples   Max number of samples (buffer limit)
 * @return Actual number of time samples generated, -1 on error.
 */
int gw_waveform_inspiral_td(GwStrainSeries *s,
    const GwBinaryParams *params,
    double f_start, double f_end,
    double dt, size_t max_samples);

/* ================================================================
 * L5 — Frequency-Domain Waveform (Stationary Phase)
 * ================================================================ */

/**
 * Generate frequency-domain inspiral waveform via Stationary Phase
 * Approximation (SPA):
 *
 *   h(f) = A(f) exp[ i Psi(f) ]
 *
 * where A(f) is the Newtonian amplitude and Psi(f) = TaylorF2 phase.
 *
 * @param h_re     Output real part [nf]
 * @param h_im     Output imag part [nf]
 * @param f_arr    Frequency array [Hz], length nf
 * @param nf       Number of frequencies
 * @param params   Binary parameters
 * @param t_c      Coalescence time [s]
 * @param phi_c    Coalescence phase [rad]
 */
void gw_waveform_fd_spa(double *h_re, double *h_im,
    const double *f_arr, int nf,
    const GwBinaryParams *params,
    double t_c, double phi_c);

/* ================================================================
 * L5 — IMR (Inspiral-Merger-Ringdown) Phenomenological Waveform
 * ================================================================ */

/**
 * IMRPhenomD-like frequency-domain waveform.
 *
 * Piecewise: inspiral (TaylorF2) | intermediate | merger-ringdown (Lorentzian)
 *
 * h(f) = A_IMR(f) exp[ i Phi_IMR(f) ]
 *
 * @param h_re, h_im  Output arrays [nf]
 * @param f_arr       Frequency array [Hz], length nf
 * @param nf          Number of frequency points
 * @param params      Binary parameters
 * @param t_c         Coalescence time [s]
 * @param phi_c       Coalescence phase [rad]
 */
void gw_waveform_imr_fd(double *h_re, double *h_im,
    const double *f_arr, int nf,
    const GwBinaryParams *params,
    double t_c, double phi_c);

/* ================================================================
 * L6 — Ringdown Waveform (Quasi-Normal Modes)
 * ================================================================ */

/**
 * Damped sinusoid ringdown waveform:
 *
 *   h(t) = A_rd * exp(-t / tau) * cos(2 pi f_rd t + phi_rd)
 *
 * for a Kerr BH with mass M_f and spin a_f.
 *
 * @param s       Output strain series (pre-allocated)
 * @param M_f     Final BH mass [kg]
 * @param a_f     Final BH dimensionless spin (0..1)
 * @param A_rd    Ringdown amplitude
 * @param phi_rd  Ringdown phase [rad]
 * @param dt      Time step [s]
 * @param n       Number of samples
 */
void gw_waveform_ringdown(GwStrainSeries *s,
    double M_f, double a_f, double A_rd, double phi_rd,
    double dt, size_t n);

/**
 * QNM frequency and damping time for Kerr BH (l=2, m=2, n=0 mode).
 *
 * Fitting formula from Berti, Cardoso & Will (2006):
 *
 *   M * omega_220 = f_22(a) + i / tau_22(a)
 *
 * @param f_rd    Output: ringdown frequency [Hz]
 * @param tau     Output: damping time [s]
 * @param M_f     Final mass [kg]
 * @param a_f     Final dimensionless spin
 */
void gw_qnm_220_mode(double *f_rd, double *tau, double M_f, double a_f);

/**
 * Final mass and spin from component masses (Barausse+ 2012 fitting).
 *
 * @param M_f     Output: final mass [kg]
 * @param a_f     Output: final spin
 * @param m1, m2  Component masses [kg]
 */
void gw_final_mass_spin(double *M_f, double *a_f, double m1, double m2);

/* ================================================================
 * L6 — Higher Modes
 * ================================================================ */

/**
 * Spherical harmonic (l,m) mode amplitude for a non-precessing binary.
 *
 * h(t; theta, phi) = sum_{l,m} h_{lm}(t) _{-2}Y_{lm}(theta, phi)
 *
 * This function computes the (l,m) = (2,2), (3,3), (4,4), (2,1) modes
 * relative to the dominant (2,2) mode, at leading PN order.
 *
 * @param amp_lm   Output: amplitude of (l,m) harmonic relative to (2,2)
 * @param l, m     Mode indices
 * @param eta      Symmetric mass ratio
 * @return 0 if mode is leading-order available, -1 for vanishing modes
 */
int gw_higher_mode_amplitude(double *amp_lm, int l, int m, double eta);

/**
 * Generate frequency-domain waveform including (2,2), (3,3), (4,4), (2,1)
 * modes.
 *
 * @param h_re, h_im  Output arrays [nf]
 * @param f_arr       Frequency array [Hz]
 * @param nf          Number of frequency points
 * @param params      Binary parameters
 * @param iota        Inclination angle [rad]
 * @param t_c, phi_c  Coalescence reference
 */
void gw_waveform_higher_modes(double *h_re, double *h_im,
    const double *f_arr, int nf,
    const GwBinaryParams *params, double iota,
    double t_c, double phi_c);

/* ================================================================
 * L7 — Strain Amplitude Scaling
 * ================================================================ */

/**
 * Frequency-domain strain amplitude (Newtonian, face-on):
 *
 *   |h(f)| = (1 / D_L) * sqrt( (5/24) pi^{-4/3} )
 *          * (G M_c / c^3)^{5/6} * f^{-7/6} * c
 */
double gw_strain_amplitude_fd(double Mc, double D_L, double f);

/**
 * Characteristic strain: h_c(f) = 2 f |h(f)| for stochastic sources.
 */
double gw_characteristic_strain(double Mc, double D_L, double f);

#endif /* GW_WAVEFORMS_H */
