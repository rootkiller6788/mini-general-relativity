/**
 * @file    gw_quadrupole.h
 * @brief   Quadrupole Radiation — Einstein's Quadrupole Formula
 *
 * Reference: Einstein (1918), Wald Ch.4.4, Maggiore Ch.3
 *
 * L1 — Quadrupole moment tensor I_{ij}, reduced quadrupole Q_{ij}
 * L2 — Quadrupole formula: h_{ij}^{TT} = (2G/c^4 r) d^2 I_{ij}/dt^2
 *      Energy loss to GWs via quadrupole luminosity
 * L4 — Einstein quadrupole formula (first-order GW generation)
 */

#ifndef GW_QUADRUPOLE_H
#define GW_QUADRUPOLE_H

#include "gw_core.h"

/* ================================================================
 * L1 — Mass Quadrupole Moment
 * ================================================================ */

/**
 * Compute mass quadrupole moment I_{ij} for N point masses.
 *
 * I_{ij} = sum_a m_a ( x_a^i x_a^j - (1/3) delta^{ij} r_a^2 )
 *
 * where x_a^i is the i-th coordinate of particle a.
 *
 * @param I_out  Output quadrupole tensor
 * @param m      Mass array [kg], length N
 * @param x      Position array, flat [x0,y0,z0, x1,y1,z1, ...], length 3N
 * @param N      Number of particles
 */
void gw_quadrupole_moment(GwTensor3 *I_out, const double *m,
                          const double *x, int N);

/**
 * Compute reduced quadrupole moment Q_{ij} for N point masses.
 *
 * Q_{ij} = I_{ij} - (1/3) delta_{ij} Tr[I]
 *
 * (Q_{ij} is traceless by construction — required for TT gauge.)
 *
 * @param Q_out  Output reduced quadrupole tensor
 * @param m      Mass array [kg]
 * @param x      Position array, flat [x0,y0,z0, ...]
 * @param N      Number of particles
 */
void gw_reduced_quadrupole(GwTensor3 *Q_out, const double *m,
                           const double *x, int N);

/* ================================================================
 * L1 — Mass Multipole Moments
 * ================================================================ */

/**
 * Current (mass-current) quadrupole moment J_{ij} for N point masses.
 *
 * J_{ij} = sum_a m_a ( x_a^i (x_a x v_a)^j + ... ) antisymmetrized
 *
 * @param J_out  Output current quadrupole
 * @param m      Mass array [kg]
 * @param x      Position array [m], flat [x0,y0,z0,...]
 * @param v      Velocity array [m/s], flat [vx0,vy0,vz0,...]
 * @param N      Number of particles
 */
void gw_current_quadrupole(GwTensor3 *J_out, const double *m,
                           const double *x, const double *v, int N);

/**
 * Mass octupole moment (L=3) for N point masses.
 *
 * M_{ijk} = sum_a m_a x_a^i x_a^j x_a^k  (raw, not STF)
 *
 * @param M_out  Output octupole as 27 doubles [i][j][k] flat
 * @param m      Mass array [kg]
 * @param x      Position array [m], flat
 * @param N      Number of particles
 */
void gw_mass_octupole(double M_out[27], const double *m,
                      const double *x, int N);

/* ================================================================
 * L2 — Quadrupole Formula (Einstein 1918)
 * ================================================================ */

/**
 * Compute GW strain h_{ij}^{TT} from second time derivative of
 * the quadrupole moment.
 *
 *   h_{ij}^{TT}(t,r) = (2G / c^4 r) * d^2 Q_{ij}(t - r/c) / dt^2
 *
 * where Q_{ij} is the reduced (traceless) quadrupole.
 *
 * @param h_tt     Output TT-gauge strain tensor (dimensionless)
 * @param Q_ddot   Second time derivative of reduced quadrupole [kg m^2 / s^2]
 * @param r        Distance to source [m]
 */
void gw_quadrupole_strain(GwTensor3 *h_tt, const GwTensor3 *Q_ddot, double r);

/**
 * Compute quadrupole strain directly from particle state at retarded time.
 * Combines reduced-quadrupole + strain computation.
 *
 * @param h_tt   Output strain tensor
 * @param m      Mass array [kg]
 * @param x      Position array [m]
 * @param a      Acceleration array [m/s^2]
 * @param N      Number of particles
 * @param r      Distance to observer [m]
 *
 * Uses: d^2 I_{ij}/dt^2 = 2 sum m_a (v_a^i v_a^j + x_a^i a_a^j)
 */
void gw_quadrupole_strain_from_state(GwTensor3 *h_tt,
    const double *m, const double *x, const double *v,
    const double *a, int N, double r);

/* ================================================================
 * L2 — Quadrupole Luminosity (GW Power)
 * ================================================================ */

/**
 * GW luminosity from third time derivative of quadrupole:
 *
 *   L_GW = (G / 5 c^5) * < d^3 Q_{ij}/dt^3 d^3 Q^{ij}/dt^3 >
 *
 * where <...> denotes time averaging over several GW cycles.
 *
 * @param Q_dddot_ij  Instantaneous third derivative (no averaging)
 * @return Instantaneous GW power [W]
 *
 * Note: For periodic systems, average over one period for physical L_GW.
 */
double gw_quadrupole_luminosity(const GwTensor3 *Q_dddot);

/**
 * Time-averaged GW luminosity for an array of Q_dddot samples.
 *
 * @param Q_dddot_arr  Array of Q_dddot tensors, one per time step
 * @param n_samples    Number of time samples
 * @return Averaged GW luminosity [W]
 */
double gw_quadrupole_luminosity_avg(const GwTensor3 *Q_dddot_arr, int n_samples);

/**
 * GW flux at distance r: F_GW = L_GW / (4 pi r^2) [W/m^2]
 */
double gw_flux(double L_GW, double r);

/* ================================================================
 * L2 — Energy-Momentum of GWs (Isaacson Tensor)
 * ================================================================ */

/**
 * Effective stress-energy tensor of GWs (Isaacson 1968):
 *
 *   T_{mu nu}^{GW} = (c^4 / 32 pi G) < partial_mu h_{alpha beta}
 *                                         partial_nu h^{alpha beta} >
 *
 * For spatial components in TT gauge:
 *   T_{00}^{GW} = (c^2 / 32 pi G) < h_ij^{TT} h^{ij,TT} >_dot
 *
 * @param h_plus  h_plus amplitude (peak)
 * @param h_cross h_cross amplitude (peak)
 * @param omega   Angular frequency [rad/s]
 * @return Energy density rho_GW [J/m^3]
 */
double gw_energy_density(double h_plus, double h_cross, double omega);

/* ================================================================
 * L4 — Quadrupole Formula for Simple Systems
 * ================================================================ */

/**
 * Quadratic (h_plus, h_cross) for a binary in circular orbit
 * (Newtonian order, quadrupole formula):
 *
 *   h_plus  = -(2 G^2 M mu) / (c^4 r a) * (1+cos^2 iota) * cos(2 Phi)
 *   h_cross = -(4 G^2 M mu) / (c^4 r a) * cos(iota) * sin(2 Phi)
 *
 * where M = m1+m2, mu = m1*m2/M, a = separation,
 * Phi = orbital phase = omega*t + phi0.
 *
 * @param h_plus   Output + strain
 * @param h_cross  Output x strain
 * @param M        Total mass [kg]
 * @param mu       Reduced mass [kg]
 * @param a        Orbital separation [m]
 * @param r        Distance to observer [m]
 * @param iota     Inclination [rad]
 * @param Phi      Orbital phase [rad]
 */
void gw_binary_quadrupole_strain(double *h_plus, double *h_cross,
    double M, double mu, double a, double r, double iota, double Phi);

/**
 * Circular binary GW luminosity (Newtonian quadrupole):
 *
 *   L_GW = (32/5) (G^4 / c^5) (mu^2 M^3) / a^5
 *
 * @return GW luminosity [W]
 */
double gw_binary_luminosity_circular(double M, double mu, double a);

#endif /* GW_QUADRUPOLE_H */
