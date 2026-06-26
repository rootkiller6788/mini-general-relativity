/**
 * @file    gw_pulsar.h
 * @brief   Continuous Gravitational Waves from Rotating Neutron Stars
 *
 * Reference: Maggiore Vol.1 Ch.7, Jaranowski, Krolak & Schutz (1998),
 *            LIGO-Virgo-KAGRA CW search papers (O1-O3)
 *
 * L2 — Continuous GW emission from non-axisymmetric neutron stars
 * L4 — Strain from ellipticity / r-mode oscillations
 * L6 — Pulsar timing arrays, spin-down limits
 * L7 — Narrow-band CW searches, F-statistic
 */

#ifndef GW_PULSAR_H
#define GW_PULSAR_H

#include "gw_core.h"

/* ================================================================
 * L2 — Ellipticity-based GW Strain
 * ================================================================ */

/**
 * Continuous GW strain from a triaxial neutron star.
 *
 *   h_0 = (4 pi^2 G / c^4) * (epsilon I_zz f_rot^2) / D
 *
 * where epsilon = (I_xx - I_yy) / I_zz is the equatorial ellipticity,
 * I_zz is the moment of inertia about rotation axis,
 * f_rot is the rotation frequency.
 *
 * The GW frequency is f_gw = 2 f_rot for the dominant l=m=2 mode.
 *
 * @param epsilon  Equatorial ellipticity
 * @param I_zz     Moment of inertia about rotation axis [kg m^2]
 * @param f_rot    Rotation frequency [Hz]
 * @param D        Distance to source [m]
 * @return Characteristic strain amplitude h_0
 */
double gw_cw_strain(double epsilon, double I_zz, double f_rot, double D);

/**
 * GW frequency for continuous wave source: f_gw = 2 * f_rot (m=2 mode).
 */
double gw_cw_frequency(double f_rot);

/* ================================================================
 * L4 — Spin-down Limit
 * ================================================================ */

/**
 * Spin-down limit on GW strain.
 *
 * Assuming all rotational kinetic energy loss goes into GWs:
 *
 *   h_0^{sd} = sqrt( (5/2) * (G/c^3) * I_zz * |f_dot| / (D^2 f_gw) )
 *
 * This is a UPPER LIMIT: the true GW strain cannot exceed this
 * without violating energy conservation.
 *
 * @param I_zz   Moment of inertia [kg m^2]
 * @param f_gw   GW frequency [Hz] (= 2 f_rot)
 * @param f_dot  Frequency derivative (spin-down) [Hz/s] (negative)
 * @param D      Distance [m]
 * @return Spin-down limit strain
 */
double gw_spindown_limit(double I_zz, double f_gw, double f_dot, double D);

/**
 * Age-based spin-down limit (characteristic age).
 *
 *   tau = -f_rot / (2 f_dot)    (characteristic age assuming n=3 braking)
 *
 *   h_0^{age} = sqrt( (5/2) * (G/c^3) * I_zz / (D^2 * tau) ) / (2 pi f_rot)
 *
 * @param I_zz   Moment of inertia [kg m^2]
 * @param f_rot  Rotation frequency [Hz]
 * @param f_dot  Frequency derivative [Hz/s]
 * @param D      Distance [m]
 * @return Age-based spin-down limit
 */
double gw_spindown_limit_age(double I_zz, double f_rot,
                             double f_dot, double D);

/* ================================================================
 * L4 — R-mode Instability
 * ================================================================ */

/**
 * GW strain from r-mode oscillations in a neutron star.
 *
 * R-modes are fluid perturbations restored by Coriolis force.
 * The GW frequency is f_gw = (4/3) f_rot for l=m=2 r-mode.
 *
 *   h_0 = sqrt( (8 pi G / c^3) * (alpha^2 * J_tilde * M * R^3 * Omega^3 / D^2) )
 *
 * @param alpha     Mode amplitude (saturation amplitude ~ 1e-6 to 1e-3)
 * @param M         NS mass [kg]
 * @param R         NS radius [m]
 * @param f_rot     Rotation frequency [Hz]
 * @param D         Distance [m]
 * @return R-mode GW strain
 */
double gw_rmode_strain(double alpha, double M, double R,
                       double f_rot, double D);

/**
 * R-mode GW frequency: f_gw = 4/3 * f_rot for dominant l=m=2 mode.
 */
double gw_rmode_frequency(double f_rot);

/* ================================================================
 * L6 — Pulsar Timing Array (PTA)
 * ================================================================ */

/**
 * GW-induced timing residual for a single pulsar.
 *
 *   R(t) = integral_0^t dt' [ nu(t') - nu_0 ] / nu_0
 *
 * For a monochromatic GW:
 *   R(t) = (h_0 / (2 pi f_gw)) * F_plus * sin(2 pi f_gw t + phi)
 *
 * @param t       Observation time [s]
 * @param h_0     GW strain amplitude
 * @param f_gw    GW frequency [Hz]
 * @param F_plus  Pulsar antenna pattern (depends on sky location)
 * @param phi0    Initial phase [rad]
 * @return Timing residual [s]
 */
double gw_pta_timing_residual(double t, double h_0, double f_gw,
                              double F_plus, double phi0);

/**
 * Hellings-Downs correlation: expected correlation between
 * two pulsars separated by angle xi on the sky, for an isotropic
 * stochastic GW background.
 *
 *   C(xi) = (1/2) - (1/4) * (1 - cos xi) / 2
 *          + (3/2) * (1 - cos xi) / 2 * ln((1-cos xi)/2)
 *
 * (Normalized to C(0) = 1/2.)
 *
 * @param xi  Angular separation between pulsars [rad]
 * @return Hellings-Downs correlation value
 */
double gw_hellings_downs(double xi);

/* ================================================================
 * L7 — F-Statistic for CW Searches
 * ================================================================ */

/**
 * F-statistic (Jaranowski, Krolak & Schutz 1998):
 *
 * Detection statistic for continuous GWs, maximizing over
 * four unknown amplitude parameters (h_0, cos iota, psi, phi0).
 *
 * F = 2 * (data | data) maximized over extrinsic parameters.
 *
 * In Gaussian noise, 2F ~ chi^2_4 (non-central if signal present).
 *
 * @param data_re, data_im  Data in frequency domain (narrow band)
 * @param nf                Number of frequency bins
 * @param df                Frequency resolution [Hz]
 * @param psd               Noise PSD
 * @param F_plus_a, F_cross_a  Antenna patterns at time a
 * @param F_plus_b, F_cross_b  Antenna patterns at time b
 * @return F-statistic value
 */
double gw_f_statistic(const double *data_re, const double *data_im,
                      int nf, double df, const double *psd,
                      double F_plus_a, double F_cross_a,
                      double F_plus_b, double F_cross_b);

/**
 * False-alarm probability for F-statistic:
 *
 *   P_F(2F > 2F_0) = exp(-F_0)   (single trial, 4 dof chi^2)
 */
double gw_fstat_fap(double F0);

/**
 * Pulsar antenna pattern (Earth-term only):
 *
 *   F_plus  = cos(2 psi) * a(t) - sin(2 psi) * b(t)
 *   F_cross = sin(2 psi) * a(t) + cos(2 psi) * b(t)
 *
 * where a,b depend on source sky location (ra, dec) and detector
 * position. Earth-term only, no pulsar term.
 *
 * @param F_plus, F_cross  Output antenna patterns
 * @param ra     Source right ascension [rad]
 * @param dec    Source declination [rad]
 * @param psi    GW polarization angle [rad]
 */
void gw_pulsar_antenna_pattern(double *F_plus, double *F_cross,
                               double ra, double dec, double psi);

#endif /* GW_PULSAR_H */
