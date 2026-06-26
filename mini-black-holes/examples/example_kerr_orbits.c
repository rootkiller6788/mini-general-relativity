/**
 * example_kerr_orbits.c
 *
 * Demonstrates Kerr black hole orbital dynamics:
 * - ISCO radii for varying spin parameters
 * - Photon orbit and IBCO radii
 * - Penrose process efficiency
 * - Superradiance conditions
 *
 * L6 Canonical Systems: Kerr BH orbital dynamics
 * L7 Applications: Astrophysical BH spin evolution
 *
 * Build: gcc -I../include example_kerr_orbits.c -L../build -lblackholes -lm -o example_kerr_orbits
 */

#include "black_hole_metrics.h"
#include "black_hole_dynamics.h"
#include <stdio.h>

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   Kerr Black Hole Orbital Dynamics                   ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    double M = mass_to_geometric(10.0 * BH_SOLAR_MASS);
    printf("BH mass: 10 M_sun  (M_geom = %.1f km)\n\n", M / 1000.0);

    printf("──────────────────────────────────────────────────────\n");
    printf("  Characteristic Radii as Function of Spin (a/M)\n");
    printf("──────────────────────────────────────────────────────\n");
    printf("  a/M    ISCO(pro) ISCO(retro)  r_ph(pro)  r_IBCO(pro)  r_+\n");
    printf("  ────────────────────────────────────────────────────\n");

    double spins[] = {0.0, 0.3, 0.5, 0.7, 0.9, 0.99, 0.998};
    for (int i = 0; i < 7; i++) {
        double a = spins[i] * M;
        double isco_pro = kerr_isco_radius(M, a, 1);
        double isco_retro = kerr_isco_radius(M, a, 0);
        double rph_pro = kerr_photon_orbit_radius(M, a, 1);
        double ribco_pro = kerr_ibco_radius(M, a, 1);
        double rp = kerr_horizon_radius(M, a);

        printf("  %.3f  %8.2f M  %8.2f M   %8.2f M  %8.2f M   %8.2f M\n",
               spins[i],
               isco_pro / M, isco_retro / M,
               rph_pro / M, ribco_pro / M, rp / M);
    }

    printf("\n──────────────────────────────────────────────────────\n");
    printf("  Penrose Process Efficiency\n");
    printf("──────────────────────────────────────────────────────\n");
    for (int i = 0; i < 7; i++) {
        double a = spins[i] * M;
        double eff = penrose_process_efficiency(M, a);
        printf("  a/M = %.3f  →  η_max = %.1f%%\n",
               spins[i], eff * 100.0);
    }

    printf("\n──────────────────────────────────────────────────────\n");
    printf("  Superradiance Check (ω/Ω_H for m=1)\n");
    printf("──────────────────────────────────────────────────────\n");
    for (int i = 0; i < 7; i++) {
        double a = spins[i] * M;
        if (a > 0) {
            double rp = kerr_horizon_radius(M, a);
            double Omega_H = a / (2.0 * M * rp);
            printf("  a/M = %.3f  →  Ω_H = %.4f / M  (superradiant for ω < Ω_H)\n",
                   spins[i], Omega_H * M);
        }
    }

    /* Bardeen spin evolution */
    printf("\n──────────────────────────────────────────────────────\n");
    printf("  Bardeen Spin Evolution (accretion of 1 M_sun)\n");
    printf("──────────────────────────────────────────────────────\n");
    double M_kg = 10.0 * BH_SOLAR_MASS;
    double a_init = 0.0;
    double a_final;
    double M_final = bardeen_spin_evolution(M_kg, a_init,
                                            1.0 * BH_SOLAR_MASS, &a_final);
    double chi_final = a_final / mass_to_geometric(M_final);
    printf("  Initial: M=10.0 M_sun, a/M=0.00\n");
    printf("  After accreting 1 M_sun: M=%.1f M_sun, a/M=%.3f\n",
           M_final / BH_SOLAR_MASS, chi_final);

    printf("\n  Note: Thorne (1974) limit a/M ≤ 0.998\n\n");

    return 0;
}
