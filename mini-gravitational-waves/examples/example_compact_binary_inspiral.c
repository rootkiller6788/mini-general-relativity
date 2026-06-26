/**
 * @file    example_compact_binary_inspiral.c
 * @brief   Simulate the full inspiral of a compact binary
 *
 * L6 — Canonical system: compact binary coalescence
 * Evolves orbital parameters under GW emission using the Peters ODE
 * and prints the inspiral track.
 */

#include <stdio.h>
#include <math.h>
#include "gw_core.h"
#include "gw_binary.h"
#include "gw_quadrupole.h"

/**
 * RK4 integrator for orbital evolution.
 * Evolves state y = [a, e, Phi] from t=0 until e < 0.001 or small a.
 */
static int integrate_inspiral(double *t_out, double *a_out, double *e_out,
                              double *f_out, int max_steps,
                              double a0, double e0,
                              double M_total, double mu) {
    double y[3] = {a0, e0, 0.0};
    double t = 0.0;
    double dt = 1.0;  /* initial time step [s] */

    for (int step = 0; step < max_steps; step++) {
        t_out[step] = t;
        a_out[step] = y[0];
        e_out[step] = y[1];

        double omega = sqrt(GW_G * M_total / (y[0]*y[0]*y[0]));
        f_out[step] = omega / GW_PI;  /* GW frequency ~ 2*f_orb, but f_orb here is per omega/pi */

        /* Adaptive time step: smaller as orbit shrinks */
        double T_orb = 2.0 * GW_PI / omega;
        dt = T_orb / 100.0;
        if (dt > 10.0) dt = 10.0;
        if (dt < 1e-6) dt = 1e-6;

        /* Stopping condition */
        if (y[0] < gw_schwarzschild_radius(M_total) * 3.0) break;
        if (y[1] < 0.001) break;

        /* RK4 step */
        double k1[3], k2[3], k3[3], k4[3], y_temp[3];

        gw_orbital_evolution_rhs(k1, y, M_total, mu);

        for (int i=0;i<3;i++) y_temp[i] = y[i] + 0.5*dt*k1[i];
        gw_orbital_evolution_rhs(k2, y_temp, M_total, mu);

        for (int i=0;i<3;i++) y_temp[i] = y[i] + 0.5*dt*k2[i];
        gw_orbital_evolution_rhs(k3, y_temp, M_total, mu);

        for (int i=0;i<3;i++) y_temp[i] = y[i] + dt*k3[i];
        gw_orbital_evolution_rhs(k4, y_temp, M_total, mu);

        for (int i=0;i<3;i++)
            y[i] += dt/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

        t += dt;
    }
    return 0;
}

int main(void) {
    printf("============================================================\n");
    printf("  Compact Binary Inspiral under GW Emission\n");
    printf("  Peters (1964) orbital evolution equations\n");
    printf("============================================================\n\n");

    /* System: 1.4+1.4 Msun NS-NS binary at 1000 km separation */
    double m1 = 1.4 * GW_MSUN;
    double m2 = 1.4 * GW_MSUN;
    double M = m1 + m2;
    double mu = m1 * m2 / M;
    double a0 = 1.0e8;  /* 100,000 km initial separation */
    double e0 = 0.1;    /* moderate eccentricity */

    printf("--- System Parameters ---\n");
    printf("  m1 = m2 = %.1f Msun (BNS)\n", m1/GW_MSUN);
    printf("  Initial separation: %.1e m\n", a0);
    printf("  Initial eccentricity: %.3f\n", e0);
    printf("  Initial GW frequency: %.3f Hz\n",
           sqrt(GW_G * M / (a0*a0*a0)) / GW_PI);

    /* Inspiral characteristic */
    printf("\n--- Inspiral Physics ---\n");
    printf("  GW luminosity scales as: L_GW ~ a^{-5}\n");
    printf("  Eccentricity decays faster than semi-major axis\n");
    printf("  (circularization timescale < inspiral timescale)\n");

    /* Show Peters enhancement factors */
    printf("\n--- Eccentricity Enhancement (Peters & Mathews 1963) ---\n");
    double e_vals[5] = {0.0, 0.1, 0.3, 0.5, 0.7};
    printf("  %8s  %12s  %12s\n", "e", "f(e) da/dt", "g(e) de/dt");
    printf("  %8s  %12s  %12s\n", "--------", "------------", "------------");
    for (int i = 0; i < 5; i++) {
        double fe = gw_eccentricity_enhancement_da(e_vals[i]);
        double ge = gw_eccentricity_enhancement_de(e_vals[i]);
        printf("  %8.3f  %12.2f  %12.2f\n", e_vals[i], fe, ge);
    }
    printf("  Note: at e=0.7, GW power is ~20x circular value!\n");

    /* Run short numerical integration */
    printf("\n--- Numerical Integration (first 100 steps) ---\n");
    #define NSTEPS 100
    double t_arr[NSTEPS], a_arr[NSTEPS], e_arr[NSTEPS], f_arr[NSTEPS];
    integrate_inspiral(t_arr, a_arr, e_arr, f_arr, NSTEPS, a0, e0, M, mu);

    printf("  %8s  %12s  %12s  %12s\n", "t [s]", "a [m]", "e", "f_gw [Hz]");
    printf("  %8s  %12s  %12s  %12s\n", "--------", "------------", "------------", "------------");
    for (int i = 0; i < 10; i++) {
        printf("  %8.1f  %12.4e  %12.6f  %12.4f\n",
               t_arr[i], a_arr[i], e_arr[i], f_arr[i]);
    }
    printf("  ...\n");

    /* Show that a decreases and e decreases */
    printf("\n  Initial a: %.4e m  ->  Final a (100 steps): %.4e m\n",
           a_arr[0], a_arr[99]);
    printf("  Initial e: %.6f    ->  Final e (100 steps): %.6f\n",
           e_arr[0], e_arr[99]);

    printf("\n============================================================\n");
    printf("  Key insight: GW emission circularizes and shrinks the\n");
    printf("  orbit. The frequency chirps upward — the 'chirp' signal\n");
    printf("  that LIGO/Virgo detects.\n");
    printf("============================================================\n");

    return 0;
}
