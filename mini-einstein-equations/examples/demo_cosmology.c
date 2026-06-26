/**
 * @file demo_cosmology.c
 * @brief FLRW cosmology demo: scale factor evolution, observables.
 *
 * Demonstrates:
 *   - Friedmann equation evaluation
 *   - Scale factor numerical integration
 *   - Cosmological observables: age, distances
 *   - Deceleration parameter evolution
 *   - Comparison of ΛCDM with Einstein-de Sitter
 *
 * Reference: Dodelson & Schmidt (2020), Weinberg (2008)
 *
 * Usage: ./demo_cosmology
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "cosmology.h"
#include "einstein.h"

int main(void)
{
    printf("=== FLRW Cosmology Demo ===\n\n");

    /* --- Planck 2018 ΛCDM parameters --- */
    CosmologicalParams lcdm;
    cosmology_params_init_planck(&lcdm);

    printf("--- ΛCDM Cosmological Parameters (Planck 2018) ---\n");
    printf("  H_0      = %.2f km/s/Mpc\n", lcdm.H0);
    printf("  Ω_m      = %.4f\n", lcdm.Omega_m);
    printf("  Ω_Λ      = %.4f\n", lcdm.Omega_L);
    printf("  Ω_r      = %.3e\n", lcdm.Omega_r);
    printf("  Ω_k      = %.6f\n", lcdm.Omega_k);
    printf("\n");

    /* --- Hubble parameter vs scale factor --- */
    printf("--- Hubble Parameter H(a) ---\n");
    printf("%10s %15s %15s\n", "a", "H(a) [km/s/Mpc]", "H(a)/H0");
    double a_values[] = {0.001, 0.01, 0.1, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0};
    for (int i = 0; i < 9; i++) {
        double a = a_values[i];
        double H = friedmann_hubble(&lcdm, a);
        printf("%10.3f %15.2f %15.4f\n", a, H, H / lcdm.H0);
    }
    printf("\n");

    /* --- Deceleration parameter --- */
    printf("--- Deceleration Parameter q(a) ---\n");
    printf("%10s %15s\n", "a", "q(a)");
    printf("  q > 0: decelerating; q < 0: accelerating\n");
    for (int i = 0; i < 9; i++) {
        double a = a_values[i];
        double q = deceleration_parameter(&lcdm, a);
        printf("%10.3f %15.4f %s\n", a, q, q > 0 ? "(decel.)" : "(accel.)");
    }
    printf("\n");

    /* --- Effective equation of state --- */
    printf("--- Effective Equation of State w_eff(a) ---\n");
    printf("%10s %15s\n", "a", "w_eff");
    for (int i = 0; i < 9; i++) {
        double a = a_values[i];
        double w = effective_eos_parameter(&lcdm, a);
        printf("%10.3f %15.4f\n", a, w);
    }
    printf("\n");

    /* --- Scale factor integration --- */
    printf("--- Scale Factor Evolution a(t) ---\n");
    int n_steps = 200;
    double t_max = 2.0; /* Hubble times */
    double *t_arr = malloc((n_steps + 1) * sizeof(double));
    double *a_arr = malloc((n_steps + 1) * sizeof(double));

    scale_factor_rk4_integrate(&lcdm, 1e-4, t_max, n_steps, t_arr, a_arr);

    printf("t (Hubble times) → a(t):\n");
    for (int i = 0; i <= n_steps; i += n_steps / 20) {
        printf("  t=%.3f  a=%.4f\n", t_arr[i], a_arr[i]);
    }
    printf("\n");

    /* --- Age of universe --- */
    printf("--- Age of the Universe ---\n");
    double age = age_of_universe(&lcdm);
    printf("  ΛCDM:  t_0 = %.2f Gyr\n", age);

    /* Einstein-de Sitter (Ω_m=1, Ω_Λ=0) for comparison */
    CosmologicalParams eds;
    cosmology_params_init_lcdm(&eds, PLANCK_H0, 1.0, 0.0);
    double age_eds = age_of_universe(&eds);
    printf("  EdS (Ω_m=1): t_0 = %.2f Gyr (would be too young!)\n", age_eds);
    printf("\n");

    /* --- Cosmological distances --- */
    printf("--- Cosmological Distances ---\n");
    printf("%10s %15s %15s %15s\n", "z", "d_L [Mpc]", "d_A [Mpc]", "lookback [Gyr]");
    double z_values[] = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
    for (int i = 0; i < 7; i++) {
        double z = z_values[i];
        double a = 1.0 / (1.0 + z);
        double d_L = luminosity_distance(&lcdm, a) / (MPC_TO_M);
        double d_A = angular_diameter_distance(&lcdm, a) / (MPC_TO_M);
        double t_lb = lookback_time(&lcdm, a);
        printf("%10.2f %15.1f %15.1f %15.2f\n", z, d_L, d_A, t_lb);
    }
    printf("\n");

    /* --- Critical density --- */
    printf("--- Critical Density ---\n");
    double rhoc = critical_density(&lcdm);
    printf("  ρ_c = %.3e kg/m^3\n", rhoc);
    printf("  ≈ %.2f proton masses / m^3\n", rhoc / 1.67262192369e-27);
    printf("\n");

    /* --- Inflation: de Sitter expansion --- */
    printf("--- Inflation (de Sitter) ---\n");
    double H_inf = 1e36; /* 1/s (typical GUT-scale inflation) */
    double a_init = 1.0;
    double dt = 1e-35; /* seconds */
    double a_final = de_sitter_scale_factor(a_init, H_inf, dt);
    double N = inflation_e_folds(a_init, a_final, H_inf, dt);
    printf("  de Sitter: H = %.1e 1/s, dt = %.1e s\n", H_inf, dt);
    printf("  a_final / a_init = exp(H*dt) = %.3e\n", a_final / a_init);
    printf("  N = %.2f e-folds\n", N);
    printf("  (Standard inflation: N ≥ 50-60 required)\n");

    free(t_arr);
    free(a_arr);

    printf("\n=== Demo Complete ===\n");
    return 0;
}
