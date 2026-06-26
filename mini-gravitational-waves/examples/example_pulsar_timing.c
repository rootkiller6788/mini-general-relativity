/**
 * @file    example_pulsar_timing.c
 * @brief   Pulsar Timing Array — GW detection with millisecond pulsars
 *
 * L6/L7 — Canonical system: PTA, Hellings-Downs correlation
 * Demonstrates: timing residuals, spin-down limits,
 * continuous GW search parameters.
 */

#include <stdio.h>
#include <math.h>
#include "gw_core.h"
#include "gw_pulsar.h"

int main(void) {
    printf("============================================================\n");
    printf("  Pulsar Timing Array — Gravitational Wave Detection\n");
    printf("  NANOGrav / EPTA / PPTA / CPTA\n");
    printf("============================================================\n\n");

    /* Reference pulsar: PSR J0437-4715 (closest, brightest MSP) */
    double I_zz  = 1.0e38;     /* moment of inertia [kg m^2] */
    double f_rot = 173.7;      /* rotation frequency [Hz] */
    double f_dot = -1.73e-15;  /* spin-down rate [Hz/s] */
    double D     = 156.0 * 3.086e16;  /* distance ~156 pc */
    double f_gw  = gw_cw_frequency(f_rot);

    printf("--- Reference Pulsar: PSR J0437-4715 ---\n");
    printf("  Rotation frequency f_rot: %.2f Hz\n", f_rot);
    printf("  GW frequency f_gw:       %.2f Hz\n", f_gw);
    printf("  Spin-down rate f_dot:    %.2e Hz/s\n", f_dot);
    printf("  Distance:               %.0f pc\n", D / 3.086e16);

    /* Spin-down limit */
    double h_sd = gw_spindown_limit(I_zz, f_gw, f_dot, D);
    printf("\n--- Spin-Down Limit ---\n");
    printf("  h_0^{sd} = %.2e\n", h_sd);
    printf("  This is the absolute upper limit on GW strain\n");
    printf("  for this pulsar if all spin-down goes into GWs.\n");

    double h_age = gw_spindown_limit_age(I_zz, f_rot, f_dot, D);
    printf("  Age-based limit: %.2e\n", h_age);

    /* CW strain for realistic ellipticity */
    double eps_crust = 1.0e-7;  /* realistic crustal strain */
    double h_cw = gw_cw_strain(eps_crust, I_zz, f_rot, D);
    printf("\n--- Realistic CW Strain ---\n");
    printf("  Ellipticity epsilon = %.1e\n", eps_crust);
    printf("  h_0(epsilon) = %.2e\n", h_cw);

    double eps_needed = eps_crust * h_sd / h_cw;
    printf("  Ellipticity to saturate spin-down: %.2e\n", eps_needed);
    printf("  (Requires exotic NS structure)\n");

    /* R-mode strain */
    double h_rmode = gw_rmode_strain(1.0e-4, 1.4*GW_MSUN, 1.2e4,
                                     f_rot, D);
    printf("\n--- R-mode GW Strain ---\n");
    printf("  R-mode (alpha=1e-4): h_0 = %.2e\n", h_rmode);
    printf("  f_gw (r-mode) = %.2f Hz\n", gw_rmode_frequency(f_rot));

    /* PTA timing residual */
    double t_obs = 15.0 * 365.25 * 86400.0;  /* 15 years */
    double Fp = 0.3;  /* typical antenna pattern */
    double residual = gw_pta_timing_residual(t_obs, 1e-15, 1.0e-9, Fp, 0.0);
    printf("\n--- PTA Timing Residual ---\n");
    printf("  GW with h=1e-15 at f=1 nHz over 15 years:\n");
    printf("  Residual ~ %.2e s\n", fabs(residual));
    printf("  NANOGrav 15yr: residuals ~ 100 ns from GWB\n");

    /* Hellings-Downs correlation */
    printf("\n--- Hellings-Downs Correlation (NANOGrav 15yr Detection!) ---\n");
    printf("  Angular separation xi [deg] -> C(xi)\n");
    printf("  %8s  %12s\n", "xi [deg]", "C(xi)");
    printf("  %8s  %12s\n", "--------", "------------");
    double xi_vals[7] = {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0};
    for (int i = 0; i < 7; i++) {
        double xi_rad = xi_vals[i] * GW_PI / 180.0;
        double hd = gw_hellings_downs(xi_rad);
        printf("  %8.1f  %12.6f\n", xi_vals[i], hd);
    }
    printf("\n  Hellings-Downs curve is the 'smoking gun' of GW origin:\n");
    printf("  - Positive correlation at small angles\n");
    printf("  - Slightly negative at large angles\n");
    printf("  - This unique angular pattern cannot be mimicked\n");
    printf("    by non-GW systematic effects\n");
    printf("  - NANOGrav 15yr: HD correlation detected at ~4 sigma!\n");

    printf("\n============================================================\n");
    printf("  Multi-messenger future: PTAs + LIGO + LISA will cover\n");
    printf("  nanoHz to kHz GW spectrum — opening the full GW window\n");
    printf("  on the Universe.\n");
    printf("============================================================\n");

    return 0;
}
