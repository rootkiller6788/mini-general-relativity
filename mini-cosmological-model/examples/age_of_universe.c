/**
 * @file age_of_universe.c
 * @brief Example: Compute the age of the universe for various cosmological models.
 *
 * Demonstrates how the age of the universe depends on cosmological parameters
 * (H0, Omega_m, Omega_Lambda). Compares Einstein-de Sitter, open, and LCDM.
 *
 * Usage: make examples && ./build/examples/age_of_universe
 */

#include <stdio.h>
#include <math.h>
#include "../include/cosmo_model.h"
#include "../include/friedmann.h"
#include "../include/distances.h"

int main(void)
{
    printf("===========================================\n");
    printf("  Age of the Universe Calculator\n");
    printf("===========================================\n\n");

    CosmoModel model;

    /* Model 1: Einstein-de Sitter (Omega_m=1, flat, no Lambda) */
    cosmo_model_init(&model, 70.0, 1.0, 0.0, 0.0, -1.0);
    double age_EdS = cosmo_age_of_universe(&model, 0.0);
    double H0_si = model.H0_si;
    double age_EdS_analytic = 2.0 / (3.0 * H0_si) / COSMO_GYR_SEC;

    printf("Model 1: Einstein-de Sitter (Omega_m=1, Omega_Lambda=0)\n");
    printf("  H0 = %.1f km/s/Mpc\n", model.H0);
    printf("  Age (numerical) = %.4f Gyr\n", age_EdS);
    printf("  Age (analytic)  = %.4f Gyr  [t0 = 2/(3*H0)]\n", age_EdS_analytic);
    printf("  This model is ruled out: universe is too young.\n\n");

    /* Model 2: Open universe (Omega_m=0.3, Omega_Lambda=0) */
    cosmo_model_init(&model, 70.0, 0.3, 0.0, 0.0, -1.0);
    double age_open = cosmo_age_of_universe(&model, 0.0);

    printf("Model 2: Open universe (Omega_m=0.3, Omega_Lambda=0, Omega_k=0.7)\n");
    printf("  H0 = %.1f km/s/Mpc\n", model.H0);
    printf("  Age = %.4f Gyr\n", age_open);
    printf("  Older than EdS but still somewhat young.\n\n");

    /* Model 3: Planck 2018 LCDM */
    cosmo_model_init_planck(&model);
    double age_LCDM = cosmo_age_of_universe(&model, 0.0);
    (void)H_z; (void)t_z;  /* Reserved for H(z), t(z) lookup table demonstration */
    double z_vals[] = {0.0, 0.5, 1.0, 2.0, 5.0};
    int i;

    printf("Model 3: Planck 2018 LCDM (Omega_m=0.315, Omega_Lambda=0.685)\n");
    printf("  H0 = %.1f km/s/Mpc\n", model.H0);
    printf("  Age of universe today = %.4f Gyr\n", age_LCDM);
    printf("  (Planck 2018 result: 13.787 +/- 0.020 Gyr)\n\n");

    printf("  Cosmic timeline:\n");
    printf("  %8s  %12s  %18s\n", "z", "t [Gyr]", "Lookback [Gyr]");

    for (i = 0; i < 5; i++) {
        double z = z_vals[i];
        double age_z = cosmo_age_of_universe(&model, z);
        double lookback = cosmo_lookback_time(&model, z);
        printf("  %8.1f  %12.4f  %18.4f\n", z, age_z, lookback);
    }

    /* Model 4: Higher H0 (SH0ES-like, H0=73) */
    cosmo_model_init(&model, 73.0, 0.3, 0.0, 0.7, -1.0);
    double age_SH0ES = cosmo_age_of_universe(&model, 0.0);

    printf("\nModel 4: Higher H0 = 73 km/s/Mpc (SH0ES-like)\n");
    printf("  Age = %.4f Gyr\n", age_SH0ES);
    printf("  Lower H0 -> older universe (Hubble tension effect).\n");

    printf("\n===========================================\n");
    printf("  Key insight: Different cosmological\n");
    printf("  parameters predict ages from ~9 to ~14 Gyr.\n");
    printf("  The LCDM model resolves the age crisis by\n");
    printf("  including dark energy (Lambda > 0).\n");
    printf("===========================================\n");

    return 0;
}
