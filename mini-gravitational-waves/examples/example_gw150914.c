/**
 * @file    example_gw150914.c
 * @brief   GW150914-like binary black hole inspiral simulation
 *
 * L6 — Canonical system: GW150914
 * Simulates the inspiral of a 36+29 Msun binary from 20 Hz to ISCO
 * and reports key numbers: chirp mass, cycles, time-to-merger, peak strain.
 * Compares with actual LIGO detection parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gw_core.h"
#include "gw_binary.h"
#include "gw_waveforms.h"
#include "gw_quadrupole.h"

int main(void) {
    printf("============================================================\n");
    printf("  GW150914-like Binary Black Hole Inspiral\n");
    printf("  Reference: LIGO-Virgo, PRL 116, 061102 (2016)\n");
    printf("============================================================\n\n");

    /* GW150914 source parameters (detector-frame approx) */
    double m1  = 36.0 * GW_MSUN;
    double m2  = 29.0 * GW_MSUN;
    double D_L = 410.0 * GW_MPC;

    /* Derived parameters */
    double M_total = m1 + m2;
    double Mc  = gw_chirp_mass(m1, m2);
    double eta = gw_symmetric_mass_ratio(m1, m2);
    double mu  = gw_reduced_mass(m1, m2);

    printf("--- Binary Parameters ---\n");
    printf("  Primary mass:       %.1f Msun\n", m1 / GW_MSUN);
    printf("  Secondary mass:     %.1f Msun\n", m2 / GW_MSUN);
    printf("  Total mass:         %.1f Msun\n", M_total / GW_MSUN);
    printf("  Chirp mass Mc:      %.1f Msun\n", Mc / GW_MSUN);
    printf("  Mass ratio eta:     %.6f\n", eta);
    printf("  Reduced mass:       %.1f Msun\n", mu / GW_MSUN);
    printf("  Luminosity distance: %.0f Mpc\n", D_L / GW_MPC);

    /* Frequency evolution */
    double f_start = 20.0;
    double f_isco  = gw_isco_frequency(M_total);
    printf("\n--- Inspiral Evolution ---\n");
    printf("  Starting frequency: %.1f Hz\n", f_start);
    printf("  ISCO frequency:     %.1f Hz\n", f_isco);

    double tau = gw_time_to_coalescence(Mc, f_start);
    printf("  Time to merger from 20 Hz: %.3f s\n", tau);
    printf("  Actual GW150914 signal: ~0.2 s in LIGO band\n");

    double N_cyc = gw_number_of_cycles(Mc, f_start, f_isco);
    printf("  GW cycles from 20 Hz to ISCO: %.0f\n", N_cyc);
    printf("  Actual GW150914: ~10 cycles visible\n");

    /* GW strain at ISCO (peak inspiral amplitude, face-on) */
    double a_isco = gw_orbital_separation(M_total, f_isco);
    double Phi = 0.0;
    double hp_isco, hx_isco;
    gw_binary_quadrupole_strain(&hp_isco, &hx_isco, M_total, mu,
                                a_isco, D_L, 0.0, Phi);
    printf("\n--- GW Strain at ISCO (face-on, peak) ---\n");
    printf("  h_plus peak:  %.2e\n", fabs(hp_isco));
    printf("  h_cross peak: %.2e\n", fabs(hx_isco));
    printf("  Actual GW150914: h_peak ~ 1e-21\n");

    /* GW luminosity at ISCO */
    double L_GW = gw_binary_luminosity_circular(M_total, mu, a_isco);
    printf("\n--- GW Luminosity at ISCO ---\n");
    printf("  L_GW: %.2e W\n", L_GW);
    printf("  ~ %.0f solar luminosities\n", L_GW / 3.828e26);
    printf("  Reference: peak GW luminosity ~ 3.6e49 W for GW150914\n");

    /* PN parameter */
    double x_isco = gw_pn_parameter(M_total, f_isco);
    printf("\n--- Post-Newtonian Regime ---\n");
    printf("  PN parameter x at ISCO: %.6f\n", x_isco);
    printf("  x = (v/c)^2 ~ 1/6 for Schwarzschild ISCO\n");

    /* Strain amplitude in frequency domain */
    double h_fd = gw_strain_amplitude_fd(Mc, D_L, f_isco);
    printf("  |h(f_isco)| ~ %.2e Hz^{-1/2}\n", h_fd);

    /* Summary */
    printf("\n============================================================\n");
    printf("  Summary: This simulation reproduces the inspiral properties\n");
    printf("  of GW150914 — the first directly detected GW signal.\n");
    printf("  Key physics: chirp mass sets frequency evolution,\n");
    printf("  quadrupole formula gives strain amplitude,\n");
    printf("  GW luminosity drives rapid inspiral.\n");
    printf("============================================================\n");

    return 0;
}
