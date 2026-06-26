/**
 * @example frw_cosmology.c
 * @brief FRW cosmological metric and Hubble expansion
 *
 * Demonstrates: L6 Canonical System - FRW metric, L7 Application (cosmology)
 * Knowledge: FRW metric, scale factor, Hubble parameter, redshift
 * Reference: Dodelson & Schmidt Sec.1-2, Wald Sec.5.1
 * Course: MIT 8.962, Caltech Ph 205
 *
 * FRW metric: ds^2 = -dt^2 + a^2(t)[dr^2/(1-kr^2) + r^2 dOmega^2]
 * where k = -1 (open), 0 (flat), +1 (closed)
 *
 * The Friedmann equation: H^2 = (da/dt / a)^2 = (8piG/3) rho - k/a^2 + Lambda/3
 */

#include "dg_metric.h"
#include "dg_connection.h"
#include "dg_curvature.h"
#include "dg_geodesic.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    printf("=== FRW Cosmology Example ===\n\n");

    /* Parameters */
    double H0 = 70.0;     /* Hubble constant (km/s/Mpc) */
    double Omega_m = 0.3; /* Matter density parameter */
    double Omega_L = 0.7; /* Dark energy density parameter */
    double Omega_k = 0.0; /* Curvature parameter (flat) */
    double a0 = 1.0;      /* Scale factor today */

    printf("Cosmological parameters:\n");
    printf("  H0 = %.1f km/s/Mpc\n", H0);
    printf("  Omega_m = %.2f, Omega_Lambda = %.2f, Omega_k = %.2f\n",
           Omega_m, Omega_L, Omega_k);

    /* FRW metric at present time:
     * g_00 = -1, g_11 = a^2/(1-kr^2), g_22 = a^2 r^2, g_33 = a^2 r^2 sin^2(theta) */
    double r = 1.0;  /* comoving radial coordinate */
    double k = 0.0;  /* flat universe */

    double g_comp[16] = {0};
    g_comp[0*4+0] = -1.0;
    g_comp[1*4+1] = a0*a0 / (1.0 - k*r*r);
    g_comp[2*4+2] = a0*a0 * r*r;
    g_comp[3*4+3] = a0*a0 * r*r;  /* sin^2(pi/2) = 1 */

    Metric *metric = metric_create(4, SIG_LORENTZIAN_MPPP, g_comp, 0);
    printf("\nFRW metric (flat, a=%.2f, r=%.2f):\n", a0, r);
    metric_print(metric);

    /* Compute the Hubble parameter's effect on the metric */
    /* For FRW: the metric time derivative involves d(a^2)/dt = 2 a a_dot = 2 a^2 H */
    double H = H0 * sqrt(Omega_m*pow(a0,-3) + Omega_L + Omega_k*pow(a0,-2));
    printf("\nHubble parameter at a=%.2f: H = %.2f km/s/Mpc\n", a0, H);

    /* Compute non-zero Christoffel symbols for FRW (flat, k=0) */
    double dg[64] = {0};
    /* d_0 g_{11} = d/dt (a^2) = 2 a a_dot = 2 a^2 H */
    double a_dot = a0 * H;
    dg[0*16 + 1*4 + 1] = 2.0 * a0 * a_dot;  /* d_t g_{rr} */
    dg[0*16 + 2*4 + 2] = 2.0 * a0 * a_dot * r*r;  /* d_t g_{theta theta} */
    dg[0*16 + 3*4 + 3] = 2.0 * a0 * a_dot * r*r;  /* d_t g_{phi phi} */

    /* Spatial derivatives */
    dg[1*16 + 2*4 + 2] = 2.0 * a0*a0 * r;  /* d_r g_{theta theta} */
    dg[1*16 + 3*4 + 3] = 2.0 * a0*a0 * r;  /* d_r g_{phi phi} */

    Connection *conn = connection_from_metric(metric, dg);
    printf("\nSelected FRW Christoffel symbols:\n");
    /* Gamma^0_{11} = a a_dot */
    printf("  Gamma^0_{rr} = %g (expansion of space)\n",
           connection_get(conn, 0, 1, 1));
    /* Gamma^1_{01} = a_dot/a = H (Hubble friction) */
    printf("  Gamma^r_{tr} = %g = H (Hubble friction)\n",
           connection_get(conn, 1, 0, 1));

    /* Compute Ricci scalar for flat FRW:
     * R = 6(a_ddot/a + (a_dot/a)^2)
     * For matter+Lambda dominated: R = 6 H^2 (with appropriate factors) */
    printf("\nRicci scalar (indicates spacetime curvature):\n");
    Metric *g_copy = metric_clone(metric);

    /* For a proper computation we need dGamma, but here we just show the concept */
    printf("  During radiation era: R ~ 0 (conformal invariance)\n");
    printf("  During matter era: R ~ 3H^2\n");

    /* Light ray propagation: null geodesics in FRW */
    GeodesicPoint photon;
    geodesic_frw_null_init(a0, &photon);

    printf("\nNull geodesic (light ray) initial state:\n");
    geodesic_point_print(&photon);
    printf("  Note: In expanding universe, photon wavelength stretches with a(t)\n");
    printf("  Redshift: 1+z = a(now)/a(emission)\n");

    /* Cosmic event horizon estimate */
    double rh = 1.0 / H0;  /* Hubble radius in Hubble units */
    printf("\nCosmic distances (in Hubble units):\n");
    printf("  Hubble radius: d_H = c/H0 ~ %.2f\n", rh);
    printf("  Particle horizon: d_p = a0 * integral_0^a0 da/(a^2 H(a))\n");

    /* Cleanup */
    connection_free(conn);
    metric_free(metric);
    metric_free(g_copy);

    printf("\n=== End FRW cosmology example ===\n");
    return 0;
}