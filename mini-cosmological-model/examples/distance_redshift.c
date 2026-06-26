/**
 * @file distance_redshift.c
 * @brief Example: Compute distance-redshift relations for LCDM cosmology.
 *
 * Demonstrates the different cosmological distance measures and how they
 * relate to each other. Shows that D_A peaks at z~1.6 (objects at high
 * redshift appear larger on the sky) while D_L grows monotonically.
 *
 * Usage: make examples && ./build/examples/distance_redshift
 */

#include <stdio.h>
#include <math.h>
#include "../include/cosmo_model.h"
#include "../include/friedmann.h"
#include "../include/distances.h"

int main(void)
{
    printf("================================================\n");
    printf("  Cosmological Distance-Redshift Relations\n");
    printf("================================================\n\n");

    CosmoModel model;
    cosmo_model_init_planck(&model);

    printf("Cosmology: Planck 2018 LCDM\n");
    printf("  H0 = %.1f km/s/Mpc\n", model.H0);
    printf("  Omega_m = %.4f\n", model.Omega_m);
    printf("  Omega_Lambda = %.4f\n", model.Omega_Lambda);
    printf("  Omega_k = %.6f\n\n", model.Omega_k);

    printf("%8s  %12s  %12s  %12s  %12s  %10s\n",
           "z", "D_C [Mpc]", "D_A [Mpc]", "D_L [Mpc]", "D_L [Glyr]", "mu [mag]");
    printf("%s\n",
           "--------  ------------  ------------  ------------  ------------  ----------");

    double z_vals[] = {0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0,
                       7.0, 10.0, 100.0, 500.0, 1000.0, 1089.0};
    int n_z = sizeof(z_vals) / sizeof(z_vals[0]);
    int i;

    for (i = 0; i < n_z; i++) {
        double z = z_vals[i];
        double D_C = distances_comoving_line_of_sight(&model, z) / COSMO_MPC;
        double D_A = distances_angular_diameter(&model, z) / COSMO_MPC;
        double D_L = distances_luminosity(&model, z) / COSMO_MPC;
        double D_L_Glyr = D_L * COSMO_MPC / (COSMO_C_LIGHT * COSMO_GYR_SEC);
        double mu = distances_distance_modulus(&model, z);

        printf("%8.2f  %12.1f  %12.1f  %12.1f  %12.2f  %10.2f\n",
               z, D_C, D_A, D_L, D_L_Glyr, mu);
    }

    printf("\n");
    printf("Key observations:\n");
    printf("  1. D_A peaks at z~1.6 and decreases at higher z - objects\n");
    printf("     at z>1.6 appear LARGER on the sky than at z=1.6!\n");
    printf("     (This is a unique prediction of expanding universe.)\n\n");
    printf("  2. D_L grows monotonically, roughly as cz/H0 at low z.\n\n");
    printf("  3. The distance modulus mu is ~43 at z=1, ~46 at z=3,\n");
    printf("     ~48 at recombination (z=1089).\n\n");
    printf("  4. SNIa at z~0.5 have mu~42 (apparent mag m~24 for\n");
    printf("     absolute mag M~-19), requiring large telescopes.\n");
    printf("================================================\n");

    return 0;
}
