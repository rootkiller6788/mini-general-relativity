/** @file test_friedmann.c
 *  @brief Test Friedmann equations against analytic solutions.
 *
 *  Verifies the correctness of H(z), E(z), scale factor evolution,
 *  and epoch transitions against known analytic results.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../include/cosmo_model.h"
#include "../include/friedmann.h"

#define TOL 1e-6

/* Test 1: Flat matter-dominated universe (Einstein-de Sitter).
 * Omega_m = 1, Omega_Lambda = 0, Omega_r = 0, Omega_k = 0.
 * E(z) = (1+z)^{3/2}
 * H(z) = H0 * (1+z)^{3/2}
 */
static void test_EdS_Ez(void)
{
    CosmoModel model;
    cosmo_model_init(&model, 70.0, 1.0, 0.0, 0.0, -1.0);

    /* E(z) should be exactly (1+z)^{3/2} */
    double z_test[] = {0.0, 0.5, 1.0, 2.0, 5.0};
    int i;
    for (i = 0; i < 5; i++) {
        double z = z_test[i];
        double E_calc = friedmann_E_z(&model, z);
        double E_expect = pow(1.0 + z, 1.5);
        assert(fabs(E_calc - E_expect) < TOL);
    }

    printf("  PASS: Einstein-de Sitter E(z) test\n");
}

/* Test 2: Flat radiation-dominated universe.
 * Omega_r = 1, all others = 0.
 * E(z) = (1+z)^2
 */
static void test_radiation_Ez(void)
{
    CosmoModel model;
    cosmo_model_init(&model, 70.0, 0.0, 1.0, 0.0, -1.0);

    double z = 3.0;
    double E_calc = friedmann_E_z(&model, z);
    double E_expect = (1.0 + z) * (1.0 + z);  /* = 16 */
    assert(fabs(E_calc - E_expect) < TOL);

    printf("  PASS: Radiation-dominated E(z) test\n");
}

/* Test 3: Flat Lambda-only universe (de Sitter).
 * Omega_Lambda = 1, all others = 0.
 * E(z) = 1 (constant)
 */
static void test_deSitter_Ez(void)
{
    CosmoModel model;
    cosmo_model_init(&model, 70.0, 0.0, 0.0, 1.0, -1.0);

    double z = 10.0;
    double E_calc = friedmann_E_z(&model, z);
    double E_expect = 1.0;
    assert(fabs(E_calc - E_expect) < TOL);

    printf("  PASS: de Sitter E(z) test\n");
}

/* Test 4: Planck 2018 LCDM - E(0) = 1.
 * At z=0, E(z) must equal 1 for any valid cosmology.
 */
static void test_planck_E0_is_1(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double E0 = friedmann_E_z(&model, 0.0);
    assert(fabs(E0 - 1.0) < TOL);

    printf("  PASS: Planck LCDM E(0)=1 test\n");
}

/* Test 5: Matter-radiation equality for EdS.
 * For Omega_m=1, Omega_r>0: z_eq = Omega_m/Omega_r - 1
 */
static void test_matter_radiation_equality(void)
{
    CosmoModel model;
    double Omega_r = 1.0e-4;
    cosmo_model_init(&model, 70.0, 1.0, Omega_r, 0.0, -1.0);

    double z_eq = friedmann_z_eq_matter_radiation(&model);
    double z_eq_expect = 1.0 / Omega_r - 1.0;  /* = 9999 */
    assert(fabs(z_eq - z_eq_expect) < 1.0);

    printf("  PASS: Matter-radiation equality test\n");
}

/* Test 6: Acceleration onset for flat LCDM.
 * For Omega_m=0.3, Omega_Lambda=0.7: z_acc ~ (2*0.7/0.3)^{1/3} - 1 ~ 0.63
 */
static void test_acceleration_onset(void)
{
    CosmoModel model;
    cosmo_model_init(&model, 70.0, 0.3, 0.0, 0.7, -1.0);

    double z_acc = friedmann_z_acceleration_onset(&model);
    double z_acc_expect = pow(2.0 * 0.7 / 0.3, 1.0/3.0) - 1.0;
    assert(fabs(z_acc - z_acc_expect) < 1.0e-6);

    printf("  PASS: Acceleration onset test (z_acc ~ %.4f, expect ~ %.4f)\n",
           z_acc, z_acc_expect);
}

/* Test 7: Density evolution rho(a) = rho0 * a^{-3(1+w)}.
 * For w=1/3 (radiation): rho ~ a^{-4}
 * For w=0 (matter): rho ~ a^{-3}
 * For w=-1 (Lambda): rho = constant
 */
static void test_density_evolution(void)
{
    double rho0 = 1.0e-26;  /* typical critical density [kg/m^3] */
    double a = 0.5;

    /* Radiation: rho(a=0.5) = rho0 * 0.5^{-4} = rho0 * 16 */
    double rho_r = friedmann_rho_a(rho0, a, 1.0/3.0);
    assert(fabs(rho_r - rho0 * 16.0) < TOL * rho0);

    /* Matter: rho(a=0.5) = rho0 * 0.5^{-3} = rho0 * 8 */
    double rho_m = friedmann_rho_a(rho0, a, 0.0);
    assert(fabs(rho_m - rho0 * 8.0) < TOL * rho0);

    /* Lambda: rho(a=0.5) = rho0 * 0.5^0 = rho0 */
    double rho_L = friedmann_rho_a(rho0, a, -1.0);
    assert(fabs(rho_L - rho0) < TOL * rho0);

    printf("  PASS: Density evolution test\n");
}

/* Test 8: Effective equation of state at z=0 for Planck LCDM.
 * w_eff(0) = Omega_r*(1/3) + Omega_m*0 + Omega_Lambda*(-1)
 *          = 9.146e-5/3 - 0.685 ~ -0.68497
 */
static void test_w_eff_planck(void)
{
    CosmoModel model;
    cosmo_model_init_planck(&model);

    double w_eff = friedmann_w_eff(&model, 0.0);
    double w_expect = model.Omega_r * (1.0/3.0) - model.Omega_Lambda;
    assert(fabs(w_eff - w_expect) < 1.0e-6);

    /* w_eff should be negative (accelerating) */
    assert(w_eff < 0.0);

    printf("  PASS: Effective EoS test (w_eff(0) = %.6f)\n", w_eff);
}

/* Test 9: Analytic scale factor for matter-dominated universe.
 * a(t) = (t/t_0)^{2/3} => a(t_0) = 1
 */
static void test_analytic_scale_matter(void)
{
    double H0_si = COSMO_H0_SI;
    double t0 = 2.0 / (3.0 * H0_si);

    double a_t0 = friedmann_scale_matter_only(t0, H0_si);
    assert(fabs(a_t0 - 1.0) < 0.02);  /* ~2% tolerance due to approximate t_0 */

    /* a(t) at t = t0/8 should be (1/8)^{2/3} = 0.25 */
    double a_half = friedmann_scale_matter_only(t0 / 8.0, H0_si);
    assert(fabs(a_half - 0.25) < 0.01);

    printf("  PASS: Analytic scale factor (matter) test\n");
}

/* Test 10: Model validation rejects bad parameters. */
static void test_validation(void)
{
    CosmoModel model_good, model_bad;

    /* Good model */
    cosmo_model_init_planck(&model_good);
    assert(cosmo_model_validate(&model_good) == 1);

    /* Bad model: negative H0 */
    cosmo_model_init(&model_bad, -70.0, 0.3, 0.0, 0.7, -1.0);
    assert(cosmo_model_validate(&model_bad) == 0);

    /* Bad model: sum rule violated */
    cosmo_model_init(&model_bad, 70.0, 2.0, 1.0, 0.5, -1.0);
    /* Sum = 3.5, should fail */
    assert(cosmo_model_validate(&model_bad) == 0);

    printf("  PASS: Model validation test\n");
}

int main(void)
{
    printf("=== Friedmann Equation Tests ===\n");

    test_EdS_Ez();
    test_radiation_Ez();
    test_deSitter_Ez();
    test_planck_E0_is_1();
    test_matter_radiation_equality();
    test_acceleration_onset();
    test_density_evolution();
    test_w_eff_planck();
    test_analytic_scale_matter();
    test_validation();

    printf("=== All %d tests passed ===\n", 10);
    return 0;
}
