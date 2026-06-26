/**
 * example_schwarzschild_thermo.c
 *
 * Demonstrates Schwarzschild black hole thermodynamics:
 * - Computes and displays full thermodynamic state
 * - Shows temperature, entropy, luminosity, and lifetime
 * - Illustrates the inverse mass-temperature relationship
 *
 * L2 Core Concepts: BH thermodynamics
 * L6 Canonical Systems: Schwarzschild BH
 *
 * Build: gcc -I../include example_schwarzschild_thermo.c -L../build -lblackholes -lm -o example_schwarzschild_thermo
 */

#include "black_hole_metrics.h"
#include "black_hole_thermodynamics.h"
#include <stdio.h>

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   Schwarzschild Black Hole Thermodynamics            ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* Three example masses spanning the observable range */
    double masses[] = {
        1e12,           /* Primordial BH (evaporating today) */
        1e22,           /* Lunar-mass BH */
        BH_SOLAR_MASS,  /* Stellar-mass BH (10 M_sun) */
        10.0 * BH_SOLAR_MASS
    };
    const char *labels[] = {
        "Primordial (10^12 kg)",
        "Lunar mass (~10^22 kg)",
        "Stellar mass (~1 M_sun)",
        "Stellar mass (10 M_sun)"
    };

    for (int i = 0; i < 4; i++) {
        SchwarzschildBH bh;
        schwarzschild_init(&bh, masses[i]);
        SchwarzschildThermo thermo;
        schwarzschild_thermo_compute(&bh, &thermo);

        printf("──────────────────────────────────────────────────────\n");
        printf("  %s\n", labels[i]);
        printf("──────────────────────────────────────────────────────\n");
        printf("  Mass:              %12.4e kg\n", thermo.mass_kg);
        printf("  Schwarzschild r_s: %12.4e m  (= %.2f km)\n",
               bh.r_s, bh.r_s / 1000.0);
        printf("  Horizon area:      %12.4e m²\n", thermo.area);
        printf("  Entropy S_BH:      %12.4e J/K  (≈ %12.4e k_B)\n",
               thermo.entropy, thermo.entropy / BH_KB);
        printf("  Temperature T_H:   %12.4e K\n", thermo.temperature);
        printf("  Luminosity P_H:    %12.4e W\n", thermo.luminosity);
        printf("  Heat capacity C:   %12.4e J/K  (NEGATIVE!)\n",
               thermo.heat_capacity);

        double tau_yr = thermo.evaporation_time / (365.25 * 86400.0);
        if (tau_yr < 1e20) {
            printf("  Lifetime τ_evap:   %12.4e yr\n", tau_yr);
        } else {
            printf("  Lifetime τ_evap:   %12.4e yr  (>> age of universe!)\n",
                   tau_yr);
        }

        /* Physical interpretation */
        printf("\n  Physics:\n");
        if (thermo.temperature > 2.725) {
            printf("    T_H > T_CMB → BH is EVAPORATING (net energy loss)\n");
        } else {
            printf("    T_H < T_CMB → BH is GROWING (accreting CMB photons)\n");
        }
        if (thermo.temperature > 1e9) {
            printf("    T_H > 10^9 K → emitting electron-positron pairs\n");
        }
        printf("\n");
    }

    /* Critical mass for CMB */
    double M_crit = critical_mass_cmb(2.725);
    printf("──────────────────────────────────────────────────────\n");
    printf("  Critical mass (T_H = T_CMB): %12.4e kg\n", M_crit);
    printf("  ≈ %.2f%% of lunar mass\n", M_crit / 7.342e22 * 100.0);
    printf("──────────────────────────────────────────────────────\n\n");

    return 0;
}
