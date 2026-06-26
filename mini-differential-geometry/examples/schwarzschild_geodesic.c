/**
 * @example schwarzschild_geodesic.c
 * @brief Compute timelike geodesic orbit in Schwarzschild spacetime
 *
 * Demonstrates: L6 Canonical System - Schwarzschild metric geodesics
 * Knowledge: geodesic integration, Schwarzschild metric, Christoffel symbols
 * Reference: Wald Sec.6.3, Carroll Sec.5.4
 * Course: MIT 8.962
 *
 * The Schwarzschild metric in Schwarzschild coordinates (t, r, theta, phi):
 *   ds^2 = -(1-2M/r)dt^2 + (1-2M/r)^{-1}dr^2 + r^2(dtheta^2 + sin^2(theta) dphi^2)
 *
 * For equatorial geodesics (theta = pi/2), constants of motion:
 *   E = (1-2M/r) dt/dtau   (energy per unit mass)
 *   L = r^2 dphi/dtau      (angular momentum per unit mass)
 */

#include "dg_geodesic.h"
#include "dg_metric.h"
#include "dg_connection.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Non-zero Christoffel symbols for Schwarzschild metric */
static void schwarzschild_connection(double r, double M, Connection *conn) {
    int dim = 4;
    /* Initialize to zero */
    for (int i = 0; i < dim*dim*dim; i++) conn->Gamma[i] = 0.0;

    double f = 1.0 - 2.0*M/r;
    double fp = 2.0*M/(r*r);  /* df/dr */

    /* Gamma^t_{tr} = Gamma^t_{rt} = M/(r^2 f) = fp/(2f) */
    double G_ttr = fp / (2.0 * f);
    connection_set(conn, 0, 0, 1, G_ttr);
    connection_set(conn, 0, 1, 0, G_ttr);

    /* Gamma^r_{tt} = (M/r^2) * f = fp*f/2 */
    connection_set(conn, 1, 0, 0, 0.5 * fp * f);

    /* Gamma^r_{rr} = -M/(r^2 f) = -fp/(2f) */
    connection_set(conn, 1, 1, 1, -fp/(2.0*f));

    /* Gamma^r_{theta theta} = -r*f */
    connection_set(conn, 1, 2, 2, -r * f);

    /* Gamma^r_{phi phi} = -r*f*sin^2(theta) */
    connection_set(conn, 1, 3, 3, -r * f);

    /* Gamma^theta_{r theta} = 1/r */
    connection_set(conn, 2, 1, 2, 1.0/r);
    connection_set(conn, 2, 2, 1, 1.0/r);

    /* Gamma^theta_{phi phi} = -sin(theta)*cos(theta) */
    connection_set(conn, 2, 3, 3, -1.0);  /* for theta=pi/2 */

    /* Gamma^phi_{r phi} = 1/r */
    connection_set(conn, 3, 1, 3, 1.0/r);
    connection_set(conn, 3, 3, 1, 1.0/r);
}

int main(void) {
    printf("=== Schwarzschild Geodesic Example ===\n\n");

    double M = 1.0;   /* Black hole mass (geometric units) */
    double r0 = 20.0; /* Starting radius */

    printf("Black hole mass M = %.2f\n", M);
    printf("Schwarzschild radius Rs = 2M = %.2f\n", 2.0*M);
    printf("Initial radius r0 = %.2f (%.1f Rs)\n", r0, r0/(2.0*M));

    /* Set up the metric at initial position */
    double g_comp[16] = {0};
    double f = 1.0 - 2.0*M/r0;
    g_comp[0*4+0] = -f;
    g_comp[1*4+1] = 1.0/f;
    g_comp[2*4+2] = r0*r0;
    g_comp[3*4+3] = r0*r0;  /* sin^2(pi/2) = 1 */
    Metric *metric = metric_create(4, SIG_LORENTZIAN_MPPP, g_comp, 0);
    printf("Metric at r=%.2f:\n", r0);
    metric_print(metric);

    /* Set up connection at initial position */
    Connection *conn = connection_alloc(4, 0);
    schwarzschild_connection(r0, M, conn);

    printf("\nNon-zero Christoffel symbols at r=%.2f:\n", r0);
    for (int rho = 0; rho < 4; rho++)
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++) {
                double G = connection_get(conn, rho, mu, nu);
                if (fabs(G) > 1e-12)
                    printf("  Gamma^%d_{%d %d} = %g\n", rho, mu, nu, G);
            }

    /* Set up geodesic initial conditions */
    GeodesicPoint start;
    geodesic_schwarzschild_init(r0, 5.0, M, &start);
    printf("\nInitial geodesic state:\n");
    geodesic_point_print(&start);

    /* Verify velocity normalization */
    double norm = geodesic_velocity_norm(&start, metric);
    printf("Velocity norm g(u,u) = %.6f (should be -1 for timelike)\n", norm);

    /* Compute angular momentum L = r^2 * dphi/dtau */
    double L = r0 * r0 * start.vel[3];
    printf("Angular momentum L = %.6f\n", L);

    /* Compute energy E = (1-2M/r) * dt/dtau */
    double E = f * start.vel[0];
    printf("Energy E = %.6f\n", E);

    /* Effective potential for Schwarzschild:
     * V_eff(r) = (1-2M/r)(1 + L^2/r^2)
     * For timelike: (dr/dtau)^2 = E^2 - V_eff(r) */
    double Veff = f * (1.0 + L*L/(r0*r0));
    printf("Effective potential V_eff(r0) = %.6f\n", Veff);
    printf("E^2 - V_eff = %.6f ", E*E - Veff);
    if (E*E - Veff > 0) {
        printf("(allowed region, particle can move)\n");
    } else {
        printf("(forbidden region)\n");
    }

    /* Integrate a few steps */
    printf("\nGeodesic integration (5 steps, dtau=0.5):\n");
    GeodesicPoint p = start;
    for (int step = 1; step <= 5; step++) {
        schwarzschild_connection(p.pos[1], M, conn);
        geodesic_step_rk4(&p, conn, 0.5);
        printf("  step %d: ", step);
        geodesic_point_print(&p);
    }

    /* Compute photon sphere radius: r = 3M */
    printf("\nPhoton sphere radius: r_ph = 3M = %.2f\n", 3.0*M);
    printf("ISCO (innermost stable circular orbit): r_isco = 6M = %.2f\n", 6.0*M);

    /* Cleanup */
    connection_free(conn);
    metric_free(metric);

    printf("\n=== End Schwarzschild example ===\n");
    return 0;
}