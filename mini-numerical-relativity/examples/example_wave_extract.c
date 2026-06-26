/**
 * @file example_wave_extract.c
 * @brief Example: Gravitational wave extraction and analysis.
 *
 * Demonstrates the full wave extraction pipeline:
 *  1. Set up Schwarzschild initial data
 *  2. Compute Ψ_4 on an extraction sphere
 *  3. Decompose into spin-weighted spherical harmonics
 *  4. Compute gravitational wave strain
 *  5. Analyze for LIGO-relevant quantities
 *
 * Run: ./example_wave_extract
 *
 * Knowledge: L2 (Newman-Penrose formalism), L5 (wave extraction),
 *            L7 (LIGO/Virgo data analysis)
 */

#include "../include/nr_wave.h"
#include "../include/nr_adm.h"
#include "../include/nr_initial.h"
#include "../include/nr_grid.h"
#include "../include/nr_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Gravitational Wave Extraction Pipeline ===\n\n");

    int N = 16;
    nr_grid_t grid;
    nr_grid_init(&grid, N, N, N, -6.0, 6.0, -6.0, 6.0, -6.0, 6.0, 3);

    /* Set up Schwarzschild */
    nr_adm_state_t *state = nr_adm_alloc(N, N, N, 3);
    nr_init_schwarzschild_isotropic(state, &grid, 1.0);

    printf("Initial data: Schwarzschild with M=1.0\n");
    printf("Grid: %d³ points\n\n", N);

    /* Extract Ψ_4 on a sphere at r = 4.0 */
    double r_ext = 4.0;
    int n_theta = 10, n_phi = 20;
    int n_total = n_theta * n_phi;

    nr_complex_t *psi4_grid = (nr_complex_t*)malloc(n_total * sizeof(nr_complex_t));
    if (!psi4_grid) {
        printf("Memory allocation failed\n");
        nr_adm_free(state);
        nr_grid_free(&grid);
        return 1;
    }

    printf("Extracting Ψ_4 on sphere r = %.1f (%d × %d angular grid)...\n",
           r_ext, n_theta, n_phi);

    nr_wave_psi4_sphere(state->gamma, state->K, &grid, grid.dx,
                         r_ext, 0.0, 0.0, 0.0,
                         n_theta, n_phi, psi4_grid);

    /* Print angular distribution of |Ψ_4| */
    printf("\n|Ψ_4(θ, φ)| on extraction sphere:\n");
    printf("%8s  %8s  %14s\n", "θ", "φ", "|Ψ_4|");
    for (int it = 0; it < n_theta; it++) {
        double theta = M_PI * (it + 0.5) / n_theta;
        for (int ip = 0; ip < n_phi; ip++) {
            double phi = 2.0 * M_PI * ip / n_phi;
            double amp = sqrt(nr_complex_norm2(psi4_grid[it * n_phi + ip]));
            if (it == n_theta/2 && ip < 5) {  /* equatorial sample */
                printf("%8.3f  %8.3f  %14.6e\n", theta, phi, amp);
            }
        }
    }

    /* Decompose into modes */
    printf("\nΨ_4 mode decomposition (l ≤ 4):\n");
    nr_psi4_modes_t modes;
    modes.lmax = 4;
    nr_wave_decompose_psi4(psi4_grid, n_theta, n_phi, 4, &modes);

    printf("%4s %4s  %14s  %14s  %14s\n", "l", "m", "Re(Ψ_4^{lm})", "Im(Ψ_4^{lm})", "|Ψ_4^{lm}|");
    for (int l = 2; l <= 4; l++) {
        for (int m = -l; m <= l; m++) {
            int idx = l * (2*4 + 1) + (m + 4);
            nr_complex_t z = modes.psi4_mode[idx];
            double amp = sqrt(nr_complex_norm2(z));
            printf("%4d %4d  %14.6e  %14.6e  %14.6e\n", l, m, z.re, z.im, amp);
        }
    }

    /* Spin-weighted harmonics verification */
    printf("\nSpin-weight -2 spherical harmonics (l=2):\n");
    printf("%4s %4s  %8s  %14s  %14s\n", "l", "m", "θ", "Re(_{-2}Y)", "Im(_{-2}Y)");
    for (int m = -2; m <= 2; m++) {
        double theta = M_PI / 4.0;  /* 45 degrees */
        nr_complex_t y = nr_wave_spin_weighted_Ylm(2, m, theta, 0.0);
        printf("%4d %4d  %8.3f  %14.6f  %14.6f\n", 2, m, theta, y.re, y.im);
    }

    /* Build a strain series */
    printf("\n--- Gravitational Wave Strain ---\n");
    nr_strain_series_t *series = nr_strain_series_alloc(4, 100);
    series->dt = 0.1;
    series->t0 = 0.0;

    /* Simulate a toy GW signal (ringdown) */
    for (int t = 0; t < 100; t++) {
        series->ntimes = t + 1;
        double time = t * series->dt;

        /* Toy Ψ_4^{22}: damped sinusoid
         * Ψ_4(t) = A e^{−t/τ} e^{i ω t}
         */
        double amp = 0.01 * exp(-time / 10.0);
        double omega = 0.37;  /* Schwarzschild QNM ω_{22} ≈ 0.3737/M */
        nr_complex_t psi4_lm = nr_complex_make(
            amp * cos(omega * time),
            amp * sin(omega * time));

        nr_wave_integrate_psi4(psi4_lm, 2, 2, time, series->dt, 0.01, series);
    }

    /* Analyze the strain */
    double hp, hc;
    nr_wave_strain_at_angle(series, 90, M_PI/4.0, 0.0, &hp, &hc);
    printf("At t=%.1f: h_+ = %.6e, h_× = %.6e\n", 9.0, hp, hc);

    double E_rad = nr_wave_radiated_energy(series);
    printf("Radiated energy: E = %.6e (G=c=1 units)\n", E_rad);

    double t_peak, omega22, tau22;
    nr_wave_final_properties(series, &t_peak, &omega22, &tau22);
    printf("Peak time: t_peak = %.2f\n", t_peak);
    printf("Ringdown frequency: ω_{22} = %.4f\n", omega22);
    printf("Damping time: τ = %.4f\n", tau22);

    /* LIGO analysis */
    double M_chirp, M_total, SNR;
    nr_wave_ligo_analysis(series, 400.0, &M_chirp, &M_total, &SNR);
    printf("\n--- LIGO/Virgo Analysis ---\n");
    printf("Distance: 400 Mpc (GW150914-like)\n");
    printf("Chirp mass: M_chirp = %.1f M_sun\n", M_chirp);
    printf("Total mass: M_total = %.1f M_sun\n", M_total);
    printf("SNR estimate: %.1f\n", SNR);

    nr_strain_series_free(series);
    free(psi4_grid);
    nr_adm_free(state);
    nr_grid_free(&grid);

    printf("\n=== Done ===\n");
    return 0;
}
