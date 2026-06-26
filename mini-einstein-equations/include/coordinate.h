/**
 * @file coordinate.h
 * @brief Coordinate systems and transformations in GR.
 *
 * Reference: Wald (1984), Ch.6; Carroll (2004), Ch.5
 *            MIT 8.962 General Relativity
 *
 * Knowledge map:
 *   L3: Coordinate transformations, Jacobian matrices
 *   L6: Schwarzschild coordinates, Kruskal-Szekeres, Eddington-Finkelstein
 *   L6: Kerr metric in Boyer-Lindquist coordinates
 *   L7: Tortoise coordinate, Regge-Wheeler coordinate
 */

#ifndef MINI_EINSTEIN_COORDINATE_H
#define MINI_EINSTEIN_COORDINATE_H

#include "metric.h"

/* ---------------------------------------------------------------------------
 * L3: General coordinate transformations
 * -------------------------------------------------------------------------*/

/**
 * Transform rank-2 covariant tensor under coordinate change x -> x'(x):
 *   T'_{mu,nu} = (∂x^{alpha}/∂x'^{mu}) * (∂x^{beta}/∂x'^{nu}) * T_{alpha,beta}
 *
 * @param J_inv   inverse Jacobian: J_inv[mu][alpha] = ∂x^{alpha}/∂x'^{mu}
 * @param T       tensor in original coordinates
 * @param T_prime output tensor in new coordinates
 */
void transform_tensor2_covariant(const Tensor2 J_inv, const Tensor2 T,
                                  Tensor2 T_prime);

/**
 * Transform Christoffel symbols:
 *   Gamma'^{lambda}_{mu,nu} = (∂x'^{lambda}/∂x^{rho})
 *     * [∂x^{sigma}/∂x'^{mu} * ∂x^{tau}/∂x'^{nu} * Gamma^{rho}_{sigma,tau}
 *        + ∂^2 x^{rho} / (∂x'^{mu} ∂x'^{nu})]
 *
 * @param J       Jacobian: J[lambda][rho] = ∂x'^{lambda}/∂x^{rho}
 * @param J_inv   inverse Jacobian
 * @param H       Hessian: H[rho][mu][nu] = ∂^2 x^{rho} / (∂x'^{mu} ∂x'^{nu})
 * @param Gamma   Christoffel symbols in original coordinates
 * @param Gamma_p output in new coordinates
 */
void transform_christoffel(const Tensor2 J, const Tensor2 J_inv,
                            const Tensor3 H, const Tensor3 Gamma,
                            Tensor3 Gamma_p);

/**
 * Compute Jacobian of coordinate transformation numerically via
 * symmetric finite differences.
 *
 * @param x        point at which to compute
 * @param transform function implementing x' = f(x)
 * @param h        step size
 * @param J        output Jacobian: J[mu][nu] = ∂x'^{mu}/∂x^{nu}
 */
typedef void (*CoordinateTransform)(const double x[4], double xp[4]);

void compute_jacobian(const double x[4], CoordinateTransform transform,
                       double h, Tensor2 J);

/* ---------------------------------------------------------------------------
 * L6: Schwarzschild coordinates
 * -------------------------------------------------------------------------*/

/**
 * Schwarzschild metric in standard (t, r, theta, phi) coordinates:
 *   ds^2 = -(1 - 2M/r) dt^2 + (1 - 2M/r)^{-1} dr^2
 *         + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
 *
 * @param M     mass parameter
 * @param r     radial coordinate (> 2M for exterior)
 * @param theta polar angle
 * @param m     output metric
 */
void schwarzschild_metric(double M, double r, double theta, Metric *m);

/**
 * Eddington-Finkelstein (ingoing) coordinates:
 * Uses v = t + r + 2M * ln|r/(2M) - 1| (advanced time)
 *
 *   ds^2 = -(1 - 2M/r) dv^2 + 2 dv dr
 *         + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
 *
 * Regular at r = 2M (event horizon).
 *
 * @param M     mass parameter
 * @param r     radial coordinate (> 0)
 * @param theta polar angle
 * @param m     output metric
 */
void eddington_finkelstein_ingoing(double M, double r, double theta, Metric *m);

/**
 * Eddington-Finkelstein (outgoing) coordinates:
 * Uses u = t - r - 2M * ln|r/(2M) - 1| (retarded time)
 *
 *   ds^2 = -(1 - 2M/r) du^2 - 2 du dr
 *         + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
 */
void eddington_finkelstein_outgoing(double M, double r, double theta, Metric *m);

/**
 * Kruskal-Szekeres maximal extension.
 * Kruskal coordinates (T, X, theta, phi) for which:
 *
 *   ds^2 = (32M^3 / r) * exp(-r/(2M)) * (-dT^2 + dX^2)
 *         + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
 *
 * where r is defined implicitly by T^2 - X^2 = (1 - r/(2M)) * exp(r/(2M)).
 */
void kruskal_szekeres_metric(double M, double T, double X, double theta,
                              Metric *m);

/**
 * Tortoise coordinate r_*:
 *   r_* = r + 2M * ln|r/(2M) - 1|
 *
 * Maps r ∈ (2M, ∞) to r_* ∈ (-∞, ∞).
 */
double tortoise_coordinate(double r, double M);

/**
 * Inverse of tortoise coordinate (numerical).
 */
double inverse_tortoise_coordinate(double rstar, double M);

/* ---------------------------------------------------------------------------
 * L6: Kerr metric
 * -------------------------------------------------------------------------*/

/**
 * Kerr metric in Boyer-Lindquist coordinates (t, r, theta, phi):
 *
 *   ds^2 = -(1 - 2Mr/rho^2) dt^2
 *         - (4Mar sin^2(theta)/rho^2) dt dphi
 *         + (rho^2/Delta) dr^2
 *         + rho^2 dtheta^2
 *         + (r^2 + a^2 + 2Ma^2 r sin^2(theta)/rho^2) sin^2(theta) dphi^2
 *
 * where:
 *   rho^2 = r^2 + a^2 cos^2(theta)
 *   Delta = r^2 - 2Mr + a^2
 *
 * @param M     mass parameter
 * @param a     spin parameter (J/M), |a| <= M
 * @param r     Boyer-Lindquist radial coordinate
 * @param theta Boyer-Lindquist polar angle
 * @param m     output metric
 */
void kerr_metric(double M, double a, double r, double theta, Metric *m);

/**
 * Kerr ergosphere surface condition: g_{tt} = 0
 *   r_ergo(theta) = M + sqrt(M^2 - a^2 cos^2(theta))
 *
 * Outside the ergosphere, timelike Killing vector exists everywhere.
 */
double kerr_ergosphere_radius(double M, double a, double theta);

/**
 * Kerr outer event horizon:
 *   r_+ = M + sqrt(M^2 - a^2)
 */
double kerr_horizon_radius(double M, double a);

/**
 * Kerr inner Cauchy horizon:
 *   r_- = M - sqrt(M^2 - a^2)
 */
double kerr_cauchy_horizon_radius(double M, double a);

/* ---------------------------------------------------------------------------
 * L7: Isotropic coordinates (for numerical relativity)
 * -------------------------------------------------------------------------*/

/**
 * Schwarzschild metric in isotropic coordinates:
 *   ds^2 = -[(1 - M/(2r_bar))/(1 + M/(2r_bar))]^2 dt^2
 *         + [1 + M/(2r_bar)]^4 * [dr_bar^2 + r_bar^2 dOmega^2]
 *
 * Useful for initial data in numerical relativity (conformally flat).
 *
 * @param M       mass parameter
 * @param r_bar   isotropic radial coordinate
 * @param theta   polar angle
 * @param m       output metric
 */
void schwarzschild_isotropic_metric(double M, double r_bar, double theta,
                                     Metric *m);

#endif /* MINI_EINSTEIN_COORDINATE_H */
