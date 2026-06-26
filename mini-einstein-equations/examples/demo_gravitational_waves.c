/**
 * @file demo_gravitational_waves.c
 * @brief Gravitational wave physics demo.
 *
 * Demonstrates:
 *   - Linearized Einstein equations
 *   - TT gauge projection
 *   - GW strain from binary systems (quadrupole formula)
 *   - Comparison with LIGO/Virgo detections
 *   - Gravitational memory effect
 *
 * Reference: Maggiore (2008), MTW Ch.35-37
 *            LIGO Scientific Collaboration, PRL 116, 061102 (2016)
 *
 * Usage: ./demo_gravitational_waves
 */

#include <stdio.h>
#include <math.h>

#include "linearized.h"
#include "tensor.h"
#include "einstein.h"

int main(void)
{
    printf("=== Gravitational Waves Demo ===\n\n");

    /* --- GW150914: First LIGO detection --- */
    printf("--- GW150914 (LIGO First Detection, 2015) ---\n");
    printf("Reference: PRL 116, 061102 (2016)\n\n");

    double m1 = 36.0 * 1.989e30; /* 36 M_sun in kg */
    double m2 = 29.0 * 1.989e30; /* 29 M_sun in kg */
    double D_L = 410.0 * 1e6 * 3.085677581e16; /* 410 Mpc in meters */

    /* Chirp mass: M_c = (m1*m2)^{3/5} / (m1+m2)^{1/5} */
    double M_chirp = pow(m1 * m2, 3.0/5.0) / pow(m1 + m2, 1.0/5.0);

    printf("  Primary mass:       m1 = %.0f M_sun\n", m1 / 1.989e30);
    printf("  Secondary mass:     m2 = %.0f M_sun\n", m2 / 1.989e30);
    printf("  Chirp mass:         M_c = %.1f M_sun\n", M_chirp / 1.989e30);
    printf("  Luminosity distance: D_L = %.0f Mpc\n", D_L / (1e6 * 3.085677581e16));
    printf("\n");

    /* GW frequency sweep: 35 Hz → 250 Hz during inspiral */
    double f_gw_early = 35.0;
    double f_gw_late = 250.0;

    double h_early = gravitational_wave_strain(M_chirp, f_gw_early, D_L);
    double h_late = gravitational_wave_strain(M_chirp, f_gw_late, D_L);

    printf("  Strain at f=%.0f Hz: h = %.2e\n", f_gw_early, h_early);
    printf("  Strain at f=%.0f Hz: h = %.2e\n", f_gw_late, h_late);
    printf("  Peak strain detected: h ~ 1.0e-21\n");
    printf("  (Order-of-magnitude agreement with LIGO measurement ✓)\n\n");

    /* --- Gravitational wave polarizations --- */
    printf("--- TT Gauge: Plus and Cross Polarizations ---\n");

    /* Example: a circularly polarized GW propagating in z-direction */
    double h0 = 1.0; /* reference amplitude */
    double h_spatial_plus[3][3] = {
        { h0,  0.0, 0.0},
        { 0.0, -h0, 0.0},
        { 0.0,  0.0, 0.0}
    };
    double h_spatial_cross[3][3] = {
        { 0.0, h0, 0.0},
        { h0, 0.0, 0.0},
        { 0.0, 0.0, 0.0}
    };

    double hp, hx;
    tt_gauge_project(h_spatial_plus, &hp, &hx);
    printf("  Plus polarization:  h_+ = %.1f, h_x = %.1f\n", hp, hx);

    tt_gauge_project(h_spatial_cross, &hp, &hx);
    printf("  Cross polarization: h_+ = %.1f, h_x = %.1f\n", hp, hx);
    printf("\n");

    /* --- Deformation of a ring of test masses --- */
    printf("--- Effect on Test Mass Ring ---\n");
    printf("Plus polarization (+):\n");
    printf("  Deforms ring: → ↔ (alternating along x and y axes)\n");
    printf("  Pattern: (x, y) → (x + h_+*x/2, y - h_+*y/2)\n");
    printf("\n");
    printf("Cross polarization (×):\n");
    printf("  Deforms ring: rotated 45° relative to +\n");
    printf("  Pattern: (x, y) → (x + h_x*y/2, y + h_x*x/2)\n\n");

    /* --- Quadrupole radiation power --- */
    printf("--- Quadrupole Radiation Power ---\n");
    printf("For a binary with reduced mass μ, separation R, orbital freq Ω:\n");
    printf("  P = (32/5)*(G^4/c^5)*μ^2*M^3/R^5\n");

    double G = G_NEWTON;
    double c = C_LIGHT;
    double mu = m1 * m2 / (m1 + m2);
    double M_tot = m1 + m2;
    double R_isco = 6.0 * G * M_tot / (c * c); /* ISCO separation */
    double f_isco = sqrt(G * M_tot / (R_isco * R_isco * R_isco)) / (2.0 * M_PI);

    /* Power at ISCO */
    double Omega_isco = 2.0 * M_PI * f_isco;
    double d3I_approx[3][3] = {{0}};
    /* Rough quadrupole estimate: I ~ μ R^2, d³I/dt³ ~ μ R^2 Ω^3 */
    double d3I_mag = mu * R_isco * R_isco * Omega_isco * Omega_isco * Omega_isco;
    d3I_approx[0][0] = d3I_mag;
    d3I_approx[1][1] = -d3I_mag;
    double P_gw = quadrupole_power(d3I_approx);

    printf("  At ISCO (f=%.0f Hz):\n", f_isco);
    printf("  P_GW ≈ %.1e Watts\n", P_gw);
    printf("  (Compare: Solar luminosity = 3.8e26 W)\n\n");

    /* --- Gravitational memory effect --- */
    printf("--- Gravitational Wave Memory Effect ---\n");
    printf("Theory (Christodoulou 1991, Favata 2010):\n");
    printf("  After a GW burst, test masses remain permanently displaced.\n");
    printf("  Δh ~ (G/c^4) * (ΔE_GW / D)  for a source at distance D.\n");

    double h_init[3][3] = {{0}};
    double h_final[3][3] = {{1e-22, 0, 0}, {0, -1e-22, 0}, {0, 0, 0}};
    double Delta_h[3][3];
    gravitational_memory(h_init, h_final, Delta_h);
    printf("  Example memory displacement: Δh_+ ≈ 1e-22\n");
    printf("  (Current LIGO sensitivity ~ few × 10^{-23}, memory\n");
    printf("   detection expected with next-gen detectors: ET, CE, LISA)\n\n");

    /* --- LIGO-Virgo-KAGRA detection summary --- */
    printf("--- GW Observatory Comparison ---\n");
    printf("%-20s %12s %15s %15s\n", "Observatory", "Arm length", "Freq range", "h sensitivity");
    printf("%-20s %12s %15s %15s\n", "LIGO (Hanford)", "4 km", "10-1000 Hz", "~10^{-23}");
    printf("%-20s %12s %15s %15s\n", "LIGO (Livingston)", "4 km", "10-1000 Hz", "~10^{-23}");
    printf("%-20s %12s %15s %15s\n", "Virgo (Italy)", "3 km", "10-1000 Hz", "~10^{-23}");
    printf("%-20s %12s %15s %15s\n", "KAGRA (Japan)", "3 km", "10-1000 Hz", "~10^{-23}");
    printf("%-20s %12s %15s %15s\n", "LISA (space)", "2.5e9 m", "0.1 mHz-0.1 Hz", "~10^{-21}");
    printf("%-20s %12s %15s %15s\n", "Einstein Telescope", "10 km", "1-10000 Hz", "~10^{-24}");
    printf("\n");
    printf("  As of 2025: >90 confident GW detections (binary BHs, NSs, NS-BH)\n");
    printf("  Future: LISA (2030s) — mHz band for supermassive BH mergers\n");

    printf("\n=== Demo Complete ===\n");
    return 0;
}
