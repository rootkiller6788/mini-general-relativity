/**
 * @file geodesic.c
 * @brief Geodesic equation integration and related computations.
 *
 * Reference: Wald (1984), Ch.3; Carroll (2004), Ch.3
 *            Hartle "Gravity" Ch.6-9
 *            Chandrasekhar "The Mathematical Theory of Black Holes" (1983)
 *
 * The geodesic equation describes the motion of freely falling test
 * particles in curved spacetime.
 */

#include "geodesic.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ========================================================================
 * L4: Geodesic equation RHS
 * ========================================================================*/

/**
 * Right-hand side of geodesic equation as 8 first-order ODEs.
 *
 * System:
 *   dx^{mu}/dtau   = u^{mu}
 *   du^{mu}/dtau   = -Γ^{mu}_{α,β} * u^{α} * u^{β}
 *
 * The state vector layout is:
 *   y[0..3] = x^{mu}   (positions)
 *   y[4..7] = u^{mu}   (4-velocity)
 *
 * Knowledge point: The geodesic equation is the generalization of
 * Newton's first law to curved spacetime. Free particles follow
 * geodesics — paths that parallel transport their own tangent vector.
 */
void geodesic_equation_rhs(ChristoffelFun Gamma_fun,
                           const GeodesicState *state,
                           double dxdt[8])
{
    /* Get Christoffel symbols at current position */
    Tensor3 Gamma;
    Gamma_fun(state->x, Gamma);

    /* Position derivatives = velocity */
    for (int mu = 0; mu < 4; mu++) {
        dxdt[mu] = state->u[mu];
    }

    /* Velocity derivatives = -Γ^{mu}_{α,β} u^{α} u^{β} */
    Vector4 accel;
    contract_gamma_uu(Gamma, state->u, accel);

    for (int mu = 0; mu < 4; mu++) {
        dxdt[4 + mu] = -accel[mu];
    }
}

/* ========================================================================
 * L5: RK4 integration (classical fourth-order Runge-Kutta)
 * ========================================================================*/

/**
 * Single RK4 step for geodesic equation.
 *
 * Knowledge point: RK4 is the workhorse of numerical ODE integration,
 * providing O(h^4) local error with four evaluations per step.
 * It naturally handles the geodesic equation's nonlinearity.
 */
void geodesic_rk4_step(ChristoffelFun Gamma_fun,
                       GeodesicState *state, double dtau)
{
    double y[8], k1[8], k2[8], k3[8], k4[8], y_temp[8];
    GeodesicState s_temp;

    /* Pack state into y */
    for (int i = 0; i < 4; i++) { y[i] = state->x[i]; }
    for (int i = 0; i < 4; i++) { y[4 + i] = state->u[i]; }

    /* k1 */
    geodesic_equation_rhs(Gamma_fun, state, k1);

    /* k2: advance by half step */
    for (int i = 0; i < 8; i++) { y_temp[i] = y[i] + 0.5 * dtau * k1[i]; }
    for (int i = 0; i < 4; i++) { s_temp.x[i] = y_temp[i]; s_temp.u[i] = y_temp[4 + i]; }
    geodesic_equation_rhs(Gamma_fun, &s_temp, k2);

    /* k3 */
    for (int i = 0; i < 8; i++) { y_temp[i] = y[i] + 0.5 * dtau * k2[i]; }
    for (int i = 0; i < 4; i++) { s_temp.x[i] = y_temp[i]; s_temp.u[i] = y_temp[4 + i]; }
    geodesic_equation_rhs(Gamma_fun, &s_temp, k3);

    /* k4: full step */
    for (int i = 0; i < 8; i++) { y_temp[i] = y[i] + dtau * k3[i]; }
    for (int i = 0; i < 4; i++) { s_temp.x[i] = y_temp[i]; s_temp.u[i] = y_temp[4 + i]; }
    geodesic_equation_rhs(Gamma_fun, &s_temp, k4);

    /* Update y */
    for (int i = 0; i < 8; i++) {
        y[i] += (dtau / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    /* Unpack */
    for (int i = 0; i < 4; i++) { state->x[i] = y[i]; state->u[i] = y[4 + i]; }
    state->tau += dtau;
}

/**
 * Integrate geodesic over N steps.
 *
 * Knowledge point: Long-duration geodesic integration is the basis for
 * ray-tracing in curved spacetime (gravitational lensing, black hole
 * imaging). The trajectory array records the worldline.
 */
void geodesic_integrate(ChristoffelFun Gamma_fun,
                        const GeodesicState *initial_state,
                        double dtau, int N,
                        GeodesicState *trajectory)
{
    trajectory[0] = *initial_state;
    for (int i = 0; i < N; i++) {
        trajectory[i + 1] = trajectory[i];
        geodesic_rk4_step(Gamma_fun, &trajectory[i + 1], dtau);

        /* Safety: stop if we're inside the horizon (r < small) or hit r=0 */
        if (trajectory[i + 1].x[1] < 1e-8) {
            /* Fill rest with last state */
            for (int j = i + 2; j <= N; j++) {
                trajectory[j] = trajectory[i + 1];
            }
            break;
        }
    }
}

/* ========================================================================
 * L6: Schwarzschild effective potential
 * ========================================================================*/

/**
 * Schwarzschild effective potential for massive particles (timelike geodesics).
 *
 * Derivation: The geodesic Lagrangian in Schwarzschild spacetime yields
 * conserved quantities E (energy) and L (angular momentum) due to the
 * static and spherical symmetries (Killing vectors ∂_t and ∂_φ).
 *
 * From the normalization condition g_{μ,ν} u^μ u^ν = -1:
 *   (1/2) (dr/dτ)^2 + V_eff(r) = (1/2) (E^2 - 1)
 *
 * where:
 *   V_eff(r) = -M/r + L^2/(2r^2) - M L^2/r^3
 *
 *   -M/r           : Newtonian potential
 *   +L^2/(2r^2)   : centrifugal barrier
 *   -M L^2/r^3     : GR correction (stronger attraction at small r)
 *
 * Knowledge point: The effective potential approach reduces the
 * 4D geodesic equation to a 1D radial problem. The GR term -M L^2/r^3
 * is what causes the ISCO at r=6M, absent in Newtonian gravity.
 */
double schwarzschild_eff_potential(double r, double M, double L)
{
    double term1 = -M / r;
    double term2 = 0.5 * L * L / (r * r);
    double term3 = -M * L * L / (r * r * r);
    return term1 + term2 + term3;
}

/**
 * Circular orbit radii from dV_eff/dr = 0.
 *
 * The derivative condition gives a quadratic in r:
 *   r^2 - (L^2/M) * r + 3 * L^2 = 0
 *
 * Solutions:
 *   r_outer = (L^2/M + sqrt(L^4/M^2 - 12 L^2)) / 2   (stable)
 *   r_inner = (L^2/M - sqrt(L^4/M^2 - 12 L^2)) / 2   (unstable)
 *
 * Real solutions require L^2 >= 12 M^2, corresponding to ISCO at r=6M.
 *
 * Knowledge point: Circular orbits in GR differ fundamentally from
 * Newtonian mechanics. The GR term creates an inner unstable orbit
 * and an ISCO inside which no stable circular orbits exist.
 */
int schwarzschild_circular_orbits(double M, double L,
                                   double *r_outer, double *r_inner)
{
    double L2 = L * L;
    double discriminant = L2 * L2 / (M * M) - 12.0 * L2;
    if (discriminant < 0.0) return 0;

    double sqrt_disc = sqrt(discriminant);
    double r_plus  = (L2 / M + sqrt_disc) / 2.0;
    double r_minus = (L2 / M - sqrt_disc) / 2.0;

    *r_outer = r_plus;
    if (discriminant > 0.0) {
        *r_inner = r_minus;
        return 2;
    }
    /* discriminant == 0 → ISCO */
    *r_inner = r_plus;
    return 1;
}

double schwarzschild_isco(double M)
{
    return 6.0 * M;
}

double schwarzschild_photon_sphere(double M)
{
    return 3.0 * M;
}

/* ========================================================================
 * L4: Geodesic deviation
 * ========================================================================*/

/**
 * Geodesic deviation (Jacobi equation):
 *   ∇_U ∇_U ξ^μ = R^μ_{ν,ρ,σ} U^ν U^ρ ξ^σ
 *
 * where:
 *   U^μ = dx^μ/dτ is the tangent to the reference geodesic,
 *   ξ^μ is the separation vector to a nearby geodesic.
 *
 * Knowledge point: Geodesic deviation is the physical manifestation of
 * spacetime curvature. The relative acceleration of nearby geodesics
 * provides a direct operational definition of the Riemann tensor.
 * This equation underlies the tidal force explanation of GR.
 */
void geodesic_deviation(const Tensor4 R_riemann,
                        const Vector4 U, const Vector4 xi,
                        Vector4 accel)
{
    for (int mu = 0; mu < 4; mu++) {
        double sum = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    sum += R_riemann[mu][nu][rho][sigma]
                         * U[nu] * U[rho] * xi[sigma];
                }
            }
        }
        accel[mu] = sum;
    }
}

/**
 * Compute norm: g_{μ,ν} u^μ u^ν.
 *
 * For timelike (massive): result = -1 (with c=1 convention)
 * For null (massless):    result = 0
 * For spacelike:          result > 0
 */
double four_velocity_norm(const Metric *m, const Vector4 u)
{
    double norm = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            norm += m->g[mu][nu] * u[mu] * u[nu];
        }
    }
    return norm;
}
