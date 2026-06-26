/**
 * @file example_schwarzschild.c
 * @brief Example: Schwarzschild black hole initial data and constraint checking.
 *
 * Sets up a Schwarzschild black hole in isotropic coordinates, verifies
 * the Hamiltonian and momentum constraints, and computes the ADM mass.
 *
 * Run: ./example_schwarzschild
 *
 * Knowledge: L4 (constraints), L6 (Schwarzschild solution)
 */

#include "../include/nr_initial.h"
#include "../include/nr_adm.h"
#include "../include/nr_grid.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Schwarzschild Black Hole Initial Data ===\n\n");

    int N = 12;
    nr_grid_t grid;
    nr_grid_init(&grid, N, N, N, -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, 3);

    printf("Grid: %d³ points, dx = %.4f\n", N, grid.dx);
    printf("Domain: [%.1f, %.1f]³\n\n", grid.xmin, grid.xmax);

    double M = 1.0;
    nr_adm_state_t *state = nr_adm_alloc(N, N, N, 3);
    nr_init_schwarzschild_isotropic(state, &grid, M);

    printf("Constructed Schwarzschild initial data with M = %.1f\n", M);
    printf("Isotropic radius of apparent horizon: r_AH = M/2 = %.2f\n", M/2.0);

    /* Sample the conformal factor along the x-axis */
    printf("\n--- Radial Profile (x-axis) ---\n");
    printf("%8s  %12s  %12s  %12s\n", "x", "psi", "alpha", "g_xx");
    int ng = 3;
    for (int i = ng; i < N + ng; i += 2) {
        double x = grid.x_coords[i];
        if (x >= 0) {
            int j = N/2 + ng, k = N/2 + ng;
            double gxx = state->gamma->comp[0]->data[
                i*state->gamma->comp[0]->stride_x
                + j*state->gamma->comp[0]->stride_y + k];
            double psi = pow(gxx, 0.25);
            double alpha = state->alpha->data[
                i*state->alpha->stride_x
                + j*state->alpha->stride_y + k];
            printf("%8.3f  %12.6f  %12.6f  %12.6f\n", x, psi, alpha, gxx);
        }
    }

    /* Check constraints */
    printf("\n--- Constraint Satisfaction ---\n");
    double norm_H, norm_M;
    nr_adm_constraint_norms(state, grid.dx, &norm_H, &norm_M);
    printf("||H||_2  = %e\n", norm_H);
    printf("||M||_2  = %e\n", norm_M);
    printf("(Both should be small; FD truncation error ~ O(dx^2))\n");

    /* ADM mass */
    nr_gf_t *psi_gf = nr_gf_alloc(N, N, N, 3);
    for (int i = ng; i < N + ng; i++) {
        for (int j = ng; j < N + ng; j++) {
            for (int k = ng; k < N + ng; k++) {
                double gxx = state->gamma->comp[0]->data[
                    i*state->gamma->comp[0]->stride_x
                    + j*state->gamma->comp[0]->stride_y + k];
                psi_gf->data[i*psi_gf->stride_x + j*psi_gf->stride_y + k]
                    = pow(gxx, 0.25);
            }
        }
    }
    double M_adm = nr_adm_mass(psi_gf, &grid, 4.0);
    printf("\nM_ADM = %.6f (expected %.1f)\n", M_adm, M);

    /* Analytic horizon: r = M/2 in isotropic coordinates */
    double r_AH = M / 2.0;
    /* ψ at horizon: ψ(r=M/2) = 1 + M/(2·M/2) = 2 */
    printf("\nAt horizon (r=%.2f): ψ = %.1f (analytic: 2.0)\n",
           r_AH, 1.0 + M/(2.0*r_AH));

    nr_gf_free(psi_gf);
    nr_adm_free(state);
    nr_grid_free(&grid);

    printf("\n=== Done ===\n");
    return 0;
}
