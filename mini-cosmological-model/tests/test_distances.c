/** @file test_distances.c
 *  @brief Test cosmological distance measures against known results.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../include/cosmo_model.h"
#include "../include/friedmann.h"
#include "../include/distances.h"

#define TOL 1e-4

/* Test 1: In a flat universe, transverse comoving = comoving line-of-sight */
static void test_flat_transverse_equals_comoving(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double z = 1.0;
    double D_C = distances_comoving_line_of_sight(&model, z);
    double D_M = distances_transverse_comoving(&model, z);

    assert(fabs(D_M - D_C) / D_C < TOL);

    printf("  PASS: Flat universe D_M = D_C (z=1)\n");
}

/* Test 2: Etherington reciprocity: D_L = (1+z)^2 * D_A */
static void test_etherington_reciprocity(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double z = 2.0;
    double D_L = distances_luminosity(&model, z);
    double D_A = distances_angular_diameter(&model, z);

    double ratio = D_L / (D_A * (1.0 + z) * (1.0 + z));
    assert(fabs(ratio - 1.0) < TOL);

    printf("  PASS: Etherington reciprocity D_L = (1+z)^2 D_A\n");
}

/* Test 3: At z=0, all distance measures approach the same value (0) */
static void test_distances_z0(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double D_C = distances_comoving_line_of_sight(&model, 0.0);
    double D_A = distances_angular_diameter(&model, 0.0);

    assert(fabs(D_C) < 1.0e-6);
    assert(fabs(D_A) < 1.0e-6);

    printf("  PASS: All distances = 0 at z=0\n");
}

/* Test 4: Distance modulus at low z should follow Hubble law.
 * For z << 1: D_L ~ cz/H0, mu ~ 5*log10(cz/H0) + 25 */
static void test_distance_modulus_low_z(void)
{
    CosmoModel model;
    cosmo_model_init(&model, 70.0, 0.3, 0.0, 0.7, -1.0);

    double z = 0.01;
    double mu = distances_distance_modulus(&model, z);
    double D_L_approx = COSMO_C_LIGHT * z / model.H0_si;
    double mu_approx = 5.0 * log10(D_L_approx / COSMO_MPC) + 25.0;

    /* Should agree within ~0.02 mag at z=0.01 */
    assert(fabs(mu - mu_approx) < 0.1);

    printf("  PASS: Low-z distance modulus (~ %.2f, approx ~ %.2f)\n", mu, mu_approx);
}

/* Test 5: Particle horizon for Planck LCDM should be finite and ~ 46 Glyr */
static void test_particle_horizon(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double d_p = distances_particle_horizon(&model, 0.0);
    double d_p_Glyr = d_p / (COSMO_C_LIGHT * COSMO_GYR_SEC);

    /* Particle horizon should be ~ 46 Glyr (comoving) */
    assert(d_p_Glyr > 40.0 && d_p_Glyr < 50.0);

    printf("  PASS: Particle horizon ~ %.1f Glyr (expect ~46 Glyr)\n", d_p_Glyr);
}

/* Test 6: Age of universe for Planck LCDM should be ~ 13.8 Gyr */
static void test_age_of_universe(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double age = cosmo_age_of_universe(&model, 0.0);

    /* Planck 2018: t_0 = 13.787 +/- 0.020 Gyr */
    assert(age > 13.0 && age < 15.0);

    printf("  PASS: Age of universe ~ %.3f Gyr (expect ~13.8 Gyr)\n", age);
}

/* Test 7: Lookback time to z=1 should be less than age */
static void test_lookback_time(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double age = cosmo_age_of_universe(&model, 0.0);
    double lookback_z1 = cosmo_lookback_time(&model, 1.0);

    assert(lookback_z1 > 0.0);
    assert(lookback_z1 < age);

    printf("  PASS: Lookback time to z=1 ~ %.3f Gyr (< age %.3f Gyr)\n",
           lookback_z1, age);
}

int main(void)
{
    printf("=== Distance Measures Tests ===\n");

    test_flat_transverse_equals_comoving();
    test_etherington_reciprocity();
    test_distances_z0();
    test_distance_modulus_low_z();
    test_particle_horizon();
    test_age_of_universe();
    test_lookback_time();

    printf("=== All %d tests passed ===\n", 7);
    return 0;
}
