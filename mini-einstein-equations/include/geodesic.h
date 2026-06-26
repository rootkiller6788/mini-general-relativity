/**
 * @file geodesic.h
 * @brief Geodesic equation integration, geodesic deviation, effective potentials.
 *
 * Reference: Wald (1984), Ch.3; Carroll (2004), Ch.3
 *            MIT 8.962 General Relativity
 *            Hartle "Gravity" Ch.6-9
 *
 * Knowledge map:
 *   L2: Geodesic principle — free particles follow geodesics
 *   L4: Geodesic equation: d^2 x^mu / dtau^2 + Gamma^mu_{alpha,beta} U^alpha U^beta = 0
 *   L5: RK4 numerical integration
 *   L6: Schwarzschild geodesics, effective potential
 */

#ifndef MINI_EINSTEIN_GEODESIC_H
#define MINI_EINSTEIN_GEODESIC_H

#include "curvature.h"

/* ---------------------------------------------------------------------------
 * L1: Geodesic state
 * -------------------------------------------------------------------------*/

/** State for geodesic integration: position + 4-velocity */
typedef struct {
    double x[4];   /**< spacetime coordinates x^{mu} = (t, r, theta, phi) */
    double u[4];   /**< 4-velocity u^{mu} = dx^{mu}/dtau */
    double tau;    /**< proper time (affine parameter) */
} GeodesicState;

/* ---------------------------------------------------------------------------
 * L4: Geodesic equation RHS
 * -------------------------------------------------------------------------*/

/**
 * Right-hand side of geodesic equation:
 *   d^2 x^{mu} / dtau^2 = -Gamma^{mu}_{alpha,beta} * u^{alpha} * u^{beta}
 *   du^{mu} / dtau      = -Gamma^{mu}_{alpha,beta} * u^{alpha} * u^{beta}
 *
 * This is a system of 8 first-order ODEs (4 positions + 4 velocities).
 *
 * @param Gamma_fun  function pointer computing Gamma at given position
 * @param state      current GeodesicState
 * @param dxdt       output: dxdt[0..3] = u^{mu}, dxdt[4..7] = a^{mu}
 */
typedef void (*ChristoffelFun)(const double x[4], Tensor3 Gamma);

void geodesic_equation_rhs(ChristoffelFun Gamma_fun,
                           const GeodesicState *state,
                           double dxdt[8]);

/* ---------------------------------------------------------------------------
 * L5: RK4 numerical integration
 * -------------------------------------------------------------------------*/

/**
 * Single RK4 step for geodesic equation.
 * Advances state by proper time interval dtau.
 *
 * @param Gamma_fun  Christoffel computer at a point
 * @param state      current state (modified in place)
 * @param dtau       proper time step
 */
void geodesic_rk4_step(ChristoffelFun Gamma_fun,
                       GeodesicState *state, double dtau);

/**
 * Integrate geodesic over N steps.
 *
 * @param Gamma_fun     Christoffel computer
 * @param initial_state initial conditions
 * @param dtau          step size
 * @param N             number of steps
 * @param trajectory    output array of N+1 states
 */
void geodesic_integrate(ChristoffelFun Gamma_fun,
                        const GeodesicState *initial_state,
                        double dtau, int N,
                        GeodesicState *trajectory);

/* ---------------------------------------------------------------------------
 * L6: Schwarzschild geodesics — effective potential
 * -------------------------------------------------------------------------*/

/**
 * Schwarzschild effective potential for massive particles:
 *   V_eff(r) = (1 - 2M/r) * (1 + L^2 / r^2)
 *
 * In geometric units (G = c = 1).
 *
 * @param r       radial coordinate (r > 2M)
 * @param M       mass parameter
 * @param L       angular momentum per unit mass
 * @return V_eff(r)
 */
double schwarzschild_eff_potential(double r, double M, double L);

/**
 * Roots of effective potential derivative — circular orbit radii.
 * Returns the two roots (stable and unstable circular orbits).
 *
 * @param M    mass parameter
 * @param L    angular momentum
 * @param r_outer output: outer (stable) circular orbit radius
 * @param r_inner output: inner (unstable) circular orbit radius
 * @return number of roots found (0, 1, or 2)
 */
int schwarzschild_circular_orbits(double M, double L,
                                   double *r_outer, double *r_inner);

/**
 * Innermost Stable Circular Orbit (ISCO):
 *   r_ISCO = 6 * M
 *
 * @param M  mass parameter
 * @return r_ISCO
 */
double schwarzschild_isco(double M);

/**
 * Photon sphere radius:
 *   r_photon = 3 * M
 *
 * Radius at which light can orbit in a circle (unstable).
 */
double schwarzschild_photon_sphere(double M);

/* ---------------------------------------------------------------------------
 * L4: Geodesic deviation equation
 * -------------------------------------------------------------------------*/

/**
 * Geodesic deviation (Jacobi equation):
 *   ∇_{U} ∇_{U} xi^{mu} = R^{mu}_{nu,rho,sigma} * U^{nu} * U^{rho} * xi^{sigma}
 *
 * Describes how nearby geodesics converge/diverge.
 *
 * @param R_riemann  Riemann tensor at the point
 * @param U          tangent vector (4-velocity of reference geodesic)
 * @param xi         deviation vector
 * @param accel      output: ∇_U ∇_U xi^{mu}
 */
void geodesic_deviation(const Tensor4 R_riemann,
                        const Vector4 U, const Vector4 xi,
                        Vector4 accel);

/**
 * Compute 4-velocity normalization:
 *   g_{mu,nu} * u^{mu} * u^{nu} = -1 (timelike, massive)
 *                                 =  0 (null, massless)
 * @return g_{mu,nu} u^{mu} u^{nu}
 */
double four_velocity_norm(const Metric *m, const Vector4 u);

#endif /* MINI_EINSTEIN_GEODESIC_H */
