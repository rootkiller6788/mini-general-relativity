/**
 * schwarzschild_numerical.h — Numerical methods for Schwarzschild spacetime
 *
 * Reference: Thijssen (2007) "Computational Physics"; Press et al. (2007) NR3
 * MIT 8.962 — General Relativity
 *
 * Covers:
 *   L5 — Computational Methods: RK4, shooting method, root finding, ray tracing
 */

#ifndef SCHWARZSCHILD_NUMERICAL_H
#define SCHWARZSCHILD_NUMERICAL_H

#include "schwarzschild_defs.h"
#include "schwarzschild_geodesics.h"

/* ==========================================================================
 * ODE Integrators for Geodesics
 * ========================================================================== */

/**
 * Single RK4 step for the 8D geodesic system.
 *
 * Given state y[8] = (t, r, theta, phi, p_t, p_r, p_theta, p_phi)
 * and affine parameter lambda, advances to lambda + dlambda.
 *
 * Complexity: O(1), 4 evaluations of schwarzschild_geodesic_derivs
 */
void schwarzschild_rk4_step_geodesic(double lambda, const double y[8],
                                       double dlambda, double y_out[8], double rs);

/**
 * Integrate a geodesic using RK4 from lambda0 to lambda_end with fixed steps.
 *
 * Stores the trajectory at each step in arrays (caller allocates).
 * Returns the number of steps taken.
 *
 * @param y0        initial state (8 components)
 * @param lambda0   initial affine parameter
 * @param lambda_end final affine parameter
 * @param n_steps   number of steps (step size = (lambda_end - lambda0) / n_steps)
 * @param rs        Schwarzschild radius
 * @param trajectory 2D array [n_steps+1][8] for output
 * @return number of steps successfully integrated (may be fewer if geodesic
 *         encounters singularity at r=0 or crosses horizon)
 */
int schwarzschild_integrate_geodesic_rk4(const double y0[8],
                                          double lambda0, double lambda_end,
                                          int n_steps, double rs,
                                          double *trajectory);

/**
 * Adaptive RK45 (Runge-Kutta-Fehlberg) step for geodesic integration.
 *
 * Takes a step and estimates the error. If error > tolerance,
 * the step is rejected and retried with a smaller step.
 *
 * @param lambda       current affine parameter (updated)
 * @param y            state vector (updated in place)
 * @param dlambda      step size (adjusted based on error)
 * @param tolerance    desired relative accuracy
 * @param rs           Schwarzschild radius
 * @return 0 if step succeeded, -1 if minimum step size reached
 */
int schwarzschild_rk45_adaptive_step(double *lambda, double y[8],
                                       double *dlambda, double tolerance,
                                       double rs);

/**
 * Integrate a geodesic using adaptive RK45 from lambda0 to lambda_end.
 *
 * The trajectory is stored as a sequence of (lambda, y[8]) pairs.
 * Caller provides output arrays and max_steps.
 */
int schwarzschild_integrate_geodesic_adaptive(const double y0[8],
                                                double lambda0, double lambda_end,
                                                double tol, double rs,
                                                double *lambda_out,
                                                double *trajectory,
                                                int max_steps,
                                                int *actual_steps);

/* ==========================================================================
 * Root Finding for Orbital Parameters
 * ========================================================================== */

/**
 * Find the circular orbit radius for given angular momentum L using
 * the secant method: solve dV_eff/dr = 0.
 *
 * Returns ISCO radius if no stable circular orbit exists for this L.
 */
double schwarzschild_find_circular_orbit_secant(double L, double rs,
                                                  double r_guess1,
                                                  double r_guess2,
                                                  double tol, int max_iter);

/**
 * Find the innermost stable circular orbit by solving d^2V_eff/dr^2 = 0
 * simultaneously with dV_eff/dr = 0.
 *
 * This solves r^2 - 6*m*r + 9*m^2 = 0 => r = 6*m - 3*m = ... no wait:
 * The ISCO condition is: r^2 - 6*m*r + 9*m^2 = 0 actually the derivative
 * of the circular orbit condition. The analytic ISCO is r = 6*m.
 *
 * This function verifies numerically.
 */
double schwarzschild_find_isco_numerical(double rs, double tol);

/**
 * Find radial turning points of a bound geodesic (periastron and apastron)
 * by solving E^2 = V_eff(r) using bisection.
 *
 * @param E_sq     energy squared (per unit mass) E^2
 * @param L_sq     angular momentum squared L^2
 * @param rs       Schwarzschild radius
 * @param r_min    lower bound for search
 * @param r_max    upper bound for search
 * @param roots    output array, roots found in ascending order
 * @return number of roots found in [r_min, r_max]
 */
int schwarzschild_find_turning_points_bisection(double E_sq, double L_sq,
                                                  double rs, double r_min,
                                                  double r_max,
                                                  double tol, double roots[4]);

/**
 * Find the photon sphere radius by solving for the unstable null circular orbit.
 * For Schwarzschild: r_ph = 3*m. Verified by solving d/dr(V_eff_null) = 0.
 */
double schwarzschild_find_photon_sphere_numerical(double rs, double tol);

/* ==========================================================================
 * Null Geodesic Ray Tracing
 * ========================================================================== */

/**
 * Ray tracing result for a single photon trajectory.
 */
typedef struct {
    int captured;
    double deflection_angle;
    double closest_approach;
    double time_of_flight;
    int n_revolutions;
} SchwarzschildRayResult;

/**
 * Trace a null geodesic from infinity towards a black hole.
 *
 * The photon starts at large r with impact parameter b.
 * If b < b_crit (3*sqrt(3)*m), the photon is captured.
 * If b > b_crit, the photon is deflected and returns to infinity.
 *
 * Uses adaptive RK45 integration with step rejection near
 * the photon sphere for accurate capture determination.
 *
 * @param impact_param  impact parameter b [m]
 * @param rs            Schwarzschild radius [m]
 * @param r_start       starting radius (large)
 * @param tol           integration tolerance
 * @param max_lambda    maximum affine parameter
 * @param result        output ray tracing result
 */
void schwarzschild_trace_null_geodesic(double impact_param, double rs,
                                         double r_start, double tol,
                                         double max_lambda,
                                         SchwarzschildRayResult *result);

/**
 * Compute the deflection angle vs impact parameter for a range of b values.
 * Generates a deflection curve that shows the logarithmic divergence
 * near the critical impact parameter.
 *
 * @param b_values      array of impact parameters
 * @param n_values      number of impact parameters
 * @param rs            Schwarzschild radius
 * @param r_start       starting radius
 * @param tol           tolerance
 * @param deflection    output deflection angles [rad]
 */
void schwarzschild_deflection_curve(const double *b_values, int n_values,
                                      double rs, double r_start, double tol,
                                      double *deflection);

/* ==========================================================================
 * Shooting Method for Bound Orbits
 * ========================================================================== */

/**
 * Use the shooting method to find initial conditions for a bound geodesic
 * with specified periastron r_p and apastron r_a.
 *
 * The shooting method adjusts E and L to match the desired turning points.
 *
 * @param r_periastron  desired periastron radius
 * @param r_apastron    desired apastron radius
 * @param rs            Schwarzschild radius
 * @param E_out         output energy
 * @param L_out         output angular momentum
 * @param tol           tolerance for turning point matching
 * @param max_iter      maximum iterations
 * @return 0 on success, 1 if no bound orbit exists for these radii
 */
int schwarzschild_shooting_bound_orbit(double r_periastron, double r_apastron,
                                         double rs, double *E_out, double *L_out,
                                         double tol, int max_iter);

/**
 * Compute the periastron advance for a full radial libration:
 * integrate dphi/dr from r_p to r_a and back.
 *
 * Delta_phi = 2 * integral_{r_p}^{r_a} (dphi/dr) dr - 2*pi
 *
 * Uses numerical quadrature with the exact expression for dphi/dr
 * from the geodesic equation.
 */
double schwarzschild_periastron_advance_numerical(double r_p, double r_a,
                                                    double rs, int n_points);

/* ==========================================================================
 * Interpolation on the Geodesic
 * ========================================================================== */

/**
 * Cubic Hermite interpolation on a geodesic trajectory.
 *
 * Given 4 consecutive points (lambda_i, r_i, dr_dlambda_i),
 * interpolates r at lambda between lambda_1 and lambda_2.
 *
 * Useful for finding exact crossing points (e.g., when a photon
 * crosses the equatorial plane).
 */
double schwarzschild_cubic_hermite_interp(double lambda,
                                            const double lambda_vals[4],
                                            const double r_vals[4],
                                            const double dr_vals[4]);

/**
 * Locate the event (r crossing r_target) along a geodesic trajectory
 * using bisection with RK4 substeps.
 *
 * @param y0         initial state
 * @param lambda0    initial lambda
 * @param lambda_max max lambda to search
 * @param r_target   target radius
 * @param rs         Schwarzschild radius
 * @param tol        tolerance on r
 * @param lambda_cross output lambda at crossing
 * @return 0 if crossing found, 1 if not found before lambda_max or capture
 */
int schwarzschild_find_radius_crossing(const double y0[8],
                                         double lambda0, double lambda_max,
                                         double r_target, double rs,
                                         double tol, double *lambda_cross);

#endif /* SCHWARZSCHILD_NUMERICAL_H */