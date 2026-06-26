/**
 * @file einstein.h
 * @brief Einstein field equations and Einstein-Hilbert action.
 *
 * Reference: Wald (1984), Ch.4, Ch.E (Appendix E)
 *            Carroll (2004), Ch.4
 *            MIT 8.962 General Relativity
 *            Weinberg "Gravitation and Cosmology" Ch.7
 *
 * Knowledge map:
 *   L4: Einstein field equations: G_{mu,nu} + Lambda * g_{mu,nu} = kappa * T_{mu,nu}
 *   L4: Einstein-Hilbert action: S = (c^4/16pi G) ∫ R sqrt(-g) d^4x
 *   L4: Trace-reversed form
 *   L7: Initial value formulation (Lichnerowicz equation)
 */

#ifndef MINI_EINSTEIN_EINSTEIN_H
#define MINI_EINSTEIN_EINSTEIN_H

#include "curvature.h"
#include "stress_energy.h"

/* ---------------------------------------------------------------------------
 * L4: Physical constants
 * -------------------------------------------------------------------------*/

/** Newton's gravitational constant in SI (m^3 kg^{-1} s^{-2}) */
#define G_NEWTON 6.67430e-11

/** Speed of light in vacuum (m/s) */
#define C_LIGHT 299792458.0

/** Einstein gravitational coupling constant:
 *   kappa = 8 * pi * G / c^4 (SI)
 *          = 8 * pi (geometric units, G=c=1) */
#define KAPPA_SI (8.0 * M_PI * G_NEWTON / (C_LIGHT * C_LIGHT * C_LIGHT * C_LIGHT))

/** Geodesic coupling constant (geometric units):
 *   kappa_geom = 8 * pi */
#define KAPPA_GEOM (8.0 * M_PI)

/* ---------------------------------------------------------------------------
 * L4: Einstein equation evaluation
 * -------------------------------------------------------------------------*/

/**
 * Einstein field equations:
 *   G_{mu,nu} + Lambda * g_{mu,nu} = kappa * T_{mu,nu}
 *
 * @param m        metric
 * @param G        Einstein tensor G_{mu,nu}
 * @param Lambda   cosmological constant
 * @param kappa    coupling constant (KAPPA_SI or KAPPA_GEOM)
 * @param T        stress-energy tensor
 * @param residual output: LHS - RHS (should be ≈ 0 for a solution)
 */
void einstein_equation_residual(const Metric *m, const Tensor2 G,
                                 double Lambda, double kappa,
                                 const Tensor2 T, Tensor2 residual);

/**
 * Compute LHS of Einstein equation:
 *   E_{mu,nu} = G_{mu,nu} + Lambda * g_{mu,nu}
 */
void einstein_lhs(const Metric *m, const Tensor2 G, double Lambda,
                  Tensor2 E);

/**
 * Compute RHS of Einstein equation:
 *   S_{mu,nu} = kappa * T_{mu,nu}
 */
void einstein_rhs(double kappa, const Tensor2 T, Tensor2 S);

/**
 * Compute trace of Einstein equation:
 *   -R + 4*Lambda = kappa * T
 * (Since G^{mu}_{mu} = -R and g^{mu}_{mu} = 4.)
 *
 * @param R_scalar  Ricci scalar
 * @param Lambda    cosmological constant
 * @param T_trace   trace of stress-energy tensor
 * @param kappa     coupling constant
 * @return residual of trace equation
 */
double einstein_trace_equation(double R_scalar, double Lambda,
                                double T_trace, double kappa);

/**
 * Compute trace-reversed Einstein equation:
 *   R_{mu,nu} = kappa * (T_{mu,nu} - (1/2) * g_{mu,nu} * T)
 *
 * Often more convenient for solving the field equations.
 *
 * @param m      metric
 * @param T      stress-energy tensor T_{mu,nu}
 * @param T_trace trace of T
 * @param kappa  coupling constant
 * @param R_ricci output: Ricci tensor = RHS of trace-reversed equation
 */
void einstein_trace_reversed(const Metric *m, const Tensor2 T,
                              double T_trace, double kappa, Tensor2 R_ricci);

/**
 * Compute Einstein tensor from stress-energy (inverse problem):
 *   G_{mu,nu} = kappa * T_{mu,nu} - Lambda * g_{mu,nu}
 */
void einstein_from_matter(const Metric *m, const Tensor2 T,
                           double kappa, double Lambda, Tensor2 G);

/* ---------------------------------------------------------------------------
 * L4: Einstein-Hilbert action
 * -------------------------------------------------------------------------*/

/**
 * Gravitational Lagrangian density:
 *   L_G = R * sqrt(-g)
 *
 * The Einstein-Hilbert action is S = ∫ L_G d^4x.
 */
double einstein_hilbert_lagrangian_density(double R_scalar, double det_g);

/**
 * Matter Lagrangian contribution to total action.
 * For a perfect fluid: L_M = -2 * rho * sqrt(-g)
 */
double matter_lagrangian_density_perfect_fluid(double rho, double det_g);

/**
 * Total action density: L = L_G + L_M
 */
double total_lagrangian_density(double R_scalar, double rho, double det_g);

/* ---------------------------------------------------------------------------
 * L7: Vacuum equation check
 * -------------------------------------------------------------------------*/

/**
 * Check if a given metric satisfies the vacuum Einstein equations:
 *   G_{mu,nu} = 0  (or  G_{mu,nu} + Lambda g_{mu,nu} = 0)
 *
 * @return max absolute component of G_{mu,nu} + Lambda g_{mu,nu}
 */
double einstein_vacuum_check(const Metric *m, const Tensor2 G, double Lambda);

/**
 * Compute the ADM mass (Komar mass) for stationary spacetimes:
 *   M_Komar = (1/8pi) * ∮_{S} ∇^{mu} xi^{nu} dS_{mu,nu}
 *
 * For Schwarzschild: equals M.
 *
 * @param m          metric
 * @param Gamma      Christoffel symbols
 * @param xi         timelike Killing vector
 * @param radius     radius of enclosing 2-sphere
 * @return Komar mass
 */
double komar_mass(const Metric *m, const Tensor3 Gamma,
                  const Vector4 xi, double radius);

#endif /* MINI_EINSTEIN_EINSTEIN_H */
