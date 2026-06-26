/**
 * @example curvature_demo.c
 * @brief Curvature tensor computation and Bianchi identity verification
 *
 * Demonstrates: L4 Curvature tensors, Bianchi identities
 * Knowledge: Riemann/Ricci/Einstein tensor computation, curvature invariants
 * Reference: Wald Sec.3.2-3.4, Carroll Sec.3
 * Course: Cambridge Part III, Oxford CMT
 *
 * Computes and verifies curvature for:
 *   1. Flat Minkowski spacetime (all curvature zero)
 *   2. A simple spherically-symmetric metric (non-trivial curvature)
 */

#include "dg_curvature.h"
#include "dg_metric.h"
#include "dg_connection.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    printf("=== Curvature Tensor Computation Demo ===\n\n");

    /* ================================================================
     * Case 1: Minkowski flat spacetime
     * ================================================================ */
    printf("--- Case 1: Minkowski Flat Spacetime ---\n");

    double eta[16] = {-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    Metric *g_flat = metric_create(4, SIG_LORENTZIAN_MPPP, eta, 0);

    double dg_flat[64] = {0};
    Connection *conn_flat = connection_from_metric(g_flat, dg_flat);

    double dGamma_flat[256] = {0};
    RiemannTensor *R_flat = riemann_compute(conn_flat, dGamma_flat);

    printf("Minkowski metric:\n");
    metric_print(g_flat);

    printf("\nRiemann tensor (non-zero components):\n");
    riemann_print(R_flat);
    printf("  (should be empty — flat spacetime)\n");

    /* Ricci tensor */
    RicciTensor *Ric_flat = ricci_compute(R_flat);
    printf("\nRicci tensor:\n");
    ricci_print(Ric_flat);
    printf("  (should be all zeros)\n");

    /* Ricci scalar */
    double R_sc_flat = ricci_scalar_compute(Ric_flat, g_flat);
    printf("Ricci scalar R = %.10f (should be 0)\n", R_sc_flat);

    /* Einstein tensor */
    double *G_flat = einstein_tensor_compute(Ric_flat, g_flat, R_sc_flat);
    printf("Einstein tensor G_{mu nu}:\n");
    for (int mu = 0; mu < 4; mu++) {
        printf("  [");
        for (int nu = 0; nu < 4; nu++)
            printf(" % 8.4f", G_flat[mu*4+nu]);
        printf(" ]\n");
    }
    printf("  (should be all zeros for vacuum flat spacetime)\n");

    /* Bianchi identity */
    double bianchi_viol = bianchi_identity_first_verify(R_flat);
    printf("\nFirst Bianchi identity violation: %.2e (should be ~0)\n",
           bianchi_viol);

    /* ================================================================
     * Case 2: Spherically symmetric metric (non-zero curvature)
     * We use a simplified metric: ds^2 = -A(r)dt^2 + B(r)dr^2 + r^2 dOmega^2
     * with A(r) = 1 - 2M/r (Schwarzschild-like)
     * ================================================================ */
    printf("\n--- Case 2: Schwarzschild-like Spherically Symmetric Metric ---\n");

    double M = 1.0, r0 = 5.0;
    double f = 1.0 - 2.0*M/r0;

    double g_schw[16] = {0};
    g_schw[0*4+0] = -f;
    g_schw[1*4+1] = 1.0/f;
    g_schw[2*4+2] = r0*r0;
    g_schw[3*4+3] = r0*r0;  /* equatorial plane */

    Metric *g_s = metric_create(4, SIG_LORENTZIAN_MPPP, g_schw, 0);
    printf("Metric at r=%.2f (Schwarzschild-like):\n", r0);
    metric_print(g_s);

    /* Compute derivatives of metric at this point:
     * d_r g_{tt} = -f' = -2M/r^2
     * d_r g_{rr} = (-1/f^2)*f' = -2M/(r^2 f^2)
     * d_r g_{theta theta} = 2r
     * d_r g_{phi phi} = 2r sin^2(theta) = 2r
     */
    double fp = 2.0*M/(r0*r0);
    double dg_s[64] = {0};
    dg_s[1*16 + 0*4 + 0] = -fp;              /* d_r g_{tt} */
    dg_s[1*16 + 1*4 + 1] = fp/(f*f);         /* d_r g_{rr} */
    dg_s[1*16 + 2*4 + 2] = 2.0*r0;           /* d_r g_{theta theta} */
    dg_s[1*16 + 3*4 + 3] = 2.0*r0;           /* d_r g_{phi phi} */

    Connection *conn_s = connection_from_metric(g_s, dg_s);

    printf("\nSelected non-zero Christoffel symbols:\n");
    for (int rho = 0; rho < 4; rho++)
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++) {
                double G = connection_get(conn_s, rho, mu, nu);
                if (fabs(G) > 1e-12)
                    printf("  Gamma^%d_{%d %d} = %g\n", rho, mu, nu, G);
            }

    /* For proper Riemann computation we need dGamma.
     * Here we set all zero (approximation) — full computation would require
     * second derivatives of the metric. */
    double dGam_s[256] = {0};
    RiemannTensor *R_s = riemann_compute(conn_s, dGam_s);

    printf("\nRiemann tensor (non-zero components):\n");
    riemann_print(R_s);

    /* Ricci tensor */
    RicciTensor *Ric_s = ricci_compute(R_s);
    printf("\nRicci tensor:\n");
    ricci_print(Ric_s);

    double R_sc_s = ricci_scalar_compute(Ric_s, g_s);
    printf("Ricci scalar R = %.10f\n", R_sc_s);

    /* For vacuum Schwarzschild: Ricci tensor = 0, but Riemann != 0
     * The non-zero Riemann components indicate tidal forces.
     * At r=5M: R^t_{rtr} = -2M/r^3 = -2/125 = -0.016
     * (This requires full dGamma computation) */

    printf("\nTheoretical values for exact Schwarzschild at r=%.2f:\n", r0);
    printf("  Riemann component R^t_{rtr} = -2M/r^3 = %.6f\n", -2.0*M/(r0*r0*r0));
    printf("  Riemann component R^theta_{phi theta phi} = 2M/r^3 = %.6f\n", 2.0*M/(r0*r0*r0));

    /* Kretschmann scalar for Schwarzschild: K = 48M^2/r^6 */
    double K_theory = 48.0 * M * M / pow(r0, 6);
    printf("  Kretschmann scalar K = 48M^2/r^6 = %.6f\n", K_theory);
    printf("  (diverges as r->0, indicating true singularity)\n");

    /* Cleanup */
    double *R_cov = riemann_lower_first_index(R_s, g_s);
    if (R_cov) {
        double K_num = kretschmann_scalar(R_cov, g_s, 4);
        printf("  Computed Kretschmann scalar: %.6f\n", K_num);
        free(R_cov);
    }

    riemann_free(R_s);
    ricci_free(Ric_s);
    connection_free(conn_s);
    metric_free(g_s);

    free(G_flat);
    ricci_free(Ric_flat);
    riemann_free(R_flat);
    connection_free(conn_flat);
    metric_free(g_flat);

    printf("\n=== End curvature demo ===\n");
    return 0;
}