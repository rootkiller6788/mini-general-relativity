/**
 * @file example_puncture.c
 * @brief Example: Binary black hole puncture initial data (Brill-Lindquist).
 *
 * Constructs initial data for two black holes using the Brill-Lindquist
 * puncture method. Verifies the Hamiltonian constraint and computes
 * the total ADM mass.
 *
 * Run: ./example_puncture
 *
 * Knowledge: L6 (Brill-Lindquist binary BH), L4 (constraint equations)
 */

#include "../include/nr_initial.h"
#include "../include/nr_adm.h"
#include "../include/nr_grid.h"
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Binary Black Hole Puncture Data ===\n\n");

    int N = 16;
    nr_grid_t grid;
    nr_grid_init(&grid, N, N, N, -6.0, 6.0, -6.0, 6.0, -6.0, 6.0, 3);

    printf("Grid: %d³ points, dx = %.4f\n", N, grid.dx);
    printf("Domain: [%.1f, %.1f]³\n\n", grid.xmin, grid.xmax);

    /* Two equal-mass black holes on x-axis */
    int n_punc = 2;
    double masses[2] = {0.5, 0.5};
    double pos_x[2]  = {2.0, -2.0};
    double pos_y[2]  = {0.0, 0.0};
    double pos_z[2]  = {0.0, 0.0};

    /* Zero momentum and spin for this example (time-symmetric data) */
    double Px[2] = {0, 0}, Py[2] = {0, 0}, Pz[2] = {0, 0};
    double Sx[2] = {0, 0}, Sy[2] = {0, 0}, Sz[2] = {0, 0};

    nr_adm_state_t *state = nr_adm_alloc(N, N, N, 3);
    nr_init_bowen_york(state, &grid, n_punc, masses,
                        pos_x, pos_y, pos_z, Px, Py, Pz, Sx, Sy, Sz);

    printf("Constructed Brill-Lindquist data for %d black holes:\n", n_punc);
    for (int p = 0; p < n_punc; p++) {
        printf("  BH %d: M=%.1f at (%.1f, %.1f, %.1f)\n",
               p, masses[p], pos_x[p], pos_y[p], pos_z[p]);
    }

    /* Sample the conformal factor along the x-axis */
    printf("\n--- Conformal Factor along x-axis ---\n");
    printf("%8s  %12s  %12s\n", "x", "ψ", "g_xx");
    int ng = 3;
    for (int i = ng; i < N + ng; i++) {
        double x = grid.x_coords[i];
        int j = N/2 + ng, k = N/2 + ng;
        double gxx = state->gamma->comp[0]->data[
            i*state->gamma->comp[0]->stride_x
            + j*state->gamma->comp[0]->stride_y + k];
        double psi = pow(gxx, 0.25);
        printf("%8.3f  %12.6f  %12.6f\n", x, psi, gxx);
    }

    /* Check constraints */
    printf("\n--- Constraint Satisfaction ---\n");
    double norm_H, norm_M;
    nr_adm_constraint_norms(state, grid.dx, &norm_H, &norm_M);
    printf("||H||_2  = %e  (Hamiltonian constraint)\n", norm_H);
    printf("||M||_2  = %e  (Momentum constraint)\n", norm_M);

    /* Estimate total ADM mass from conformal factor at large r */
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
    double M_adm = nr_adm_mass(psi_gf, &grid, 5.0);
    printf("\nM_ADM = %.6f (expected sum of bare masses = %.1f)\n",
           M_adm, masses[0] + masses[1]);

    /* Distance between black holes */
    double separation = pos_x[0] - pos_x[1];
    printf("Separation: D = %.1f\n", separation);
    printf("Approximate orbital period (Newtonian): T ≈ 2π D^{3/2} / sqrt(M_total) ≈ %.1f\n",
           2.0 * M_PI * pow(separation, 1.5) / sqrt(masses[0] + masses[1]));

    nr_gf_free(psi_gf);
    nr_adm_free(state);
    nr_grid_free(&grid);

    printf("\n=== Done ===\n");
    return 0;
}
