/**
 * @file demo_schwarzschild.c
 * @brief Full Schwarzschild geodesic orbit integration demo.
 *
 * Demonstrates:
 *   - Schwarzschild metric construction
 *   - Christoffel symbols computation
 *   - RK4 geodesic integration
 *   - Effective potential analysis
 *   - ISCO and photon sphere verification
 *
 * Reference: MIT 8.962, Wald Ch.6, Hartle Ch.9
 *
 * Usage: ./demo_schwarzschild
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "tensor.h"
#include "metric.h"
#include "curvature.h"
#include "geodesic.h"
#include "coordinate.h"
#include "einstein.h"

/* --- Christoffel provider for Schwarzschild --- */
typedef struct {
    double M;
    double r;
    double theta;
} SchwParams;

static void schw_christoffel_fun(const double x[4], Tensor3 Gamma)
{
    SchwParams *p = (SchwParams*)(void*)x; /* simplified: reuse x array as params ptr */
    /* Actually, we need x array for coordinates. Let's pass data differently. */
    double M = x[4]; /* extended: x[4]=M, x[5]=unused */
    double r = x[1];
    double theta = x[2];

    double f = 1.0 - 2.0 * M / r;
    double df_dr = 2.0 * M / (r * r);
    double st = sin(theta);
    double ct = cos(theta);

    memset(Gamma, 0, sizeof(Tensor3));

    if (r < 2.0 * M + 1e-10) return; /* inside horizon */

    /* Γ^t_{t,r} = Γ^t_{r,t} = (1/2) g^{tt} ∂_r g_{tt} = f'/(2f) */
    double G_ttr = 0.5 * df_dr / f;
    Gamma[0][0][1] = G_ttr;
    Gamma[0][1][0] = G_ttr;

    /* Γ^r_{t,t} = -(1/2) g^{rr} ∂_r g_{tt} = (f/2) * f' */
    Gamma[1][0][0] = 0.5 * f * df_dr;

    /* Γ^r_{r,r} = -(1/2) g^{rr} ∂_r g_{rr} = -f'/(2f) */
    Gamma[1][1][1] = -0.5 * df_dr / f;

    /* Γ^r_{θ,θ} = -r f */
    Gamma[1][2][2] = -r * f;

    /* Γ^r_{φ,φ} = -r f sin^2(θ) */
    Gamma[1][3][3] = -r * f * st * st;

    /* Γ^θ_{r,θ} = Γ^θ_{θ,r} = 1/r */
    Gamma[2][1][2] = 1.0 / r;
    Gamma[2][2][1] = 1.0 / r;

    /* Γ^θ_{φ,φ} = -sin(θ) cos(θ) */
    Gamma[2][3][3] = -st * ct;

    /* Γ^φ_{r,φ} = Γ^φ_{φ,r} = 1/r */
    Gamma[3][1][3] = 1.0 / r;
    Gamma[3][3][1] = 1.0 / r;

    /* Γ^φ_{θ,φ} = Γ^φ_{φ,θ} = cot(θ) */
    Gamma[3][2][3] = ct / st;
    Gamma[3][3][2] = ct / st;
}

static void schw_christoffel(const double x[4], Tensor3 Gamma)
{
    double M = 1.0; /* Default M=1 */
    double r = x[1];
    double theta = x[2];

    double f = 1.0 - 2.0 * M / r;
    double df_dr = 2.0 * M / (r * r);
    double st = sin(theta);
    double ct = cos(theta);

    memset(Gamma, 0, sizeof(Tensor3));

    if (r < 2.0 * M + 1e-10) return;

    Gamma[0][0][1] = Gamma[0][1][0] = 0.5 * df_dr / f;
    Gamma[1][0][0] = 0.5 * f * df_dr;
    Gamma[1][1][1] = -0.5 * df_dr / f;
    Gamma[1][2][2] = -r * f;
    Gamma[1][3][3] = -r * f * st * st;
    Gamma[2][1][2] = Gamma[2][2][1] = 1.0 / r;
    Gamma[2][3][3] = -st * ct;
    Gamma[3][1][3] = Gamma[3][3][1] = 1.0 / r;
    Gamma[3][2][3] = Gamma[3][3][2] = ct / st;
}

int main(void)
{
    printf("=== Schwarzschild Geodesic Orbit Demo ===\n\n");

    double M = 1.0; /* geometric units */

    /* --- Effective Potential Analysis --- */
    printf("--- Effective Potential Analysis ---\n");
    printf("Mass M = %.2f (geometric units)\n", M);
    printf("ISCO radius r_ISCO = %.2f M\n", schwarzschild_isco(M));
    printf("Photon sphere r_ph = %.2f M\n", schwarzschild_photon_sphere(M));
    printf("\n");

    /* Effective potential for different L values */
    printf("Effective potential V_eff(r) for various angular momenta L:\n");
    printf("%8s %12s %12s %12s %12s\n", "r/M", "L=3.5", "L=3.464", "L=4.0", "L=5.0");
    for (int i = 0; i < 12; i++) {
        double r = 3.0 + i * 1.0;
        printf("%8.1f %12.6f %12.6f %12.6f %12.6f\n",
               r,
               schwarzschild_eff_potential(r, M, 3.5),
               schwarzschild_eff_potential(r, M, 3.464),
               schwarzschild_eff_potential(r, M, 4.0),
               schwarzschild_eff_potential(r, M, 5.0));
    }
    printf("\n");

    /* Circular orbit analysis */
    printf("--- Circular Orbit Analysis ---\n");
    double L_values[] = {3.6, 4.0, 5.0, 10.0};
    for (int i = 0; i < 4; i++) {
        double L = L_values[i];
        double r_outer, r_inner;
        int n = schwarzschild_circular_orbits(M, L, &r_outer, &r_inner);
        printf("L = %.1f: %d orbit(s) found", L, n);
        if (n >= 1) printf(", r_outer = %.2f M", r_outer);
        if (n >= 2) printf(", r_inner = %.2f M (unstable)", r_inner);
        printf("\n");
    }
    printf("\n");

    /* --- Geodesic Integration Demo --- */
    printf("--- Geodesic Integration Demo ---\n");

    /* Circular orbit at r=10M (stable) */
    double r0 = 10.0 * M;
    double L_circ = r0 * sqrt(M / (r0 - 3.0 * M)); /* angular momentum for circular orbit */
    double E_circ = (r0 - 2.0*M) / sqrt(r0 * (r0 - 3.0*M)); /* energy */

    printf("Initial circular orbit: r0 = %.1f M, L = %.4f, E = %.4f\n",
           r0 / M, L_circ, E_circ);

    /* Set up initial state */
    double Omega = sqrt(M / (r0 * r0 * r0)); /* Keplerian angular velocity (approx) */
    GeodesicState state;
    state.x[0] = 0.0;                    /* t */
    state.x[1] = r0;                     /* r */
    state.x[2] = M_PI / 2.0;             /* theta = pi/2 (equatorial) */
    state.x[3] = 0.0;                    /* phi */
    state.u[0] = E_circ / (1.0 - 2.0*M/r0); /* dt/dtau */
    state.u[1] = 0.0;                    /* dr/dtau = 0 (circular) */
    state.u[2] = 0.0;                    /* dtheta/dtau = 0 */
    state.u[3] = L_circ / (r0 * r0);     /* dphi/dtau */
    state.tau = 0.0;

    /* Integrate one orbit period */
    double T_orbit = 2.0 * M_PI / Omega;
    int N_steps = 200;
    double dtau = T_orbit / N_steps;

    printf("Orbital period: T = %.2f M, dtau = %.4f M, N_steps = %d\n",
           T_orbit, dtau, N_steps);
    printf("\nTrajectory (first 20 steps + last 5):\n");
    printf("%8s %12s %12s %12s\n", "tau", "r", "phi", "dr/dtau");

    for (int step = 0; step <= N_steps; step++) {
        if (step <= 20 || step >= N_steps - 5) {
            printf("%8.2f %12.6f %12.6f %12.2e\n",
                   state.tau, state.x[1], state.x[3], state.u[1]);
        }
        if (step < N_steps) {
            geodesic_rk4_step(schw_christoffel, &state, dtau);
        }
    }

    /* Check radial drift after one orbit */
    double r_final = state.x[1];
    double drift = fabs(r_final - r0);
    printf("\nRadial drift after one orbit: %.2e M\n", drift);
    printf("(Should be small for stable circular orbit)\n\n");

    /* --- Curvature invariants at the orbit --- */
    printf("--- Curvature at r = %.1f M ---\n", r0);
    Metric m;
    schwarzschild_metric(M, r0, M_PI/2.0, &m);
    printf("  g_{tt} = %.6f\n", m.g[0][0]);
    printf("  g_{rr} = %.6f\n", m.g[1][1]);
    printf("  det(g) = %.6f\n", m.det_g);
    printf("  Kretschmann (expected) = 48M^2/r^6 = %.6e\n",
           48.0 * M * M / (r0 * r0 * r0 * r0 * r0 * r0));
    printf("  (Kretschmann at r=2M horizon = 3/(4M^4) = %.6f — finite!)\n",
           3.0 / (4.0 * M * M * M * M));

    printf("\n=== Demo Complete ===\n");
    return 0;
}
