/**
 * @file    test_binary.c
 * @brief   Unit tests for binary system dynamics (chirp, inspiral, PN)
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "gw_core.h"
#include "gw_binary.h"

static int passed = 0, failed = 0;
#define T(n)  printf("  TEST: %s ... ", n)
#define P()   do { printf("PASS\n"); passed++; } while(0)
#define F(m)  do { printf("FAIL: %s\n", m); failed++; } while(0)
#define C(c,m) do { if(!(c)){F(m);return;} } while(0)

static void test_chirp_mass(void) {
    T("chirp mass computation");
    double mc = gw_chirp_mass(1.4*GW_MSUN, 1.4*GW_MSUN);
    /* For equal masses m: Mc = (m^2)^{3/5} / (2m)^{1/5} = m * 2^{-1/5}
     * 2^{-0.2} ~ 0.87055, so Mc ~ 1.2188 Msun */
    double expected = 1.4 * GW_MSUN * pow(2.0, -0.2);
    C(fabs(mc/expected - 1.0) < 1e-10, "chirp mass");
    P();
}

static void test_symmetric_mass_ratio(void) {
    T("symmetric mass ratio");
    double eta = gw_symmetric_mass_ratio(30.0, 30.0);
    C(fabs(eta - 0.25) < 1e-10, "equal masses -> 0.25");

    double eta2 = gw_symmetric_mass_ratio(30.0, 10.0);
    /* eta = 300 / 1600 = 0.1875 */
    C(fabs(eta2 - 0.1875) < 1e-10, "3:1 mass ratio");
    P();
}

static void test_frequency_evolution(void) {
    T("GW frequency evolution");
    double Mc = 28.0 * GW_MSUN;
    double f  = 30.0;

    double dfdt = gw_frequency_derivative(Mc, f);
    C(dfdt > 0, "df/dt positive");
    C(dfdt < 100.0, "df/dt reasonable for BBH at 30 Hz");

    double tau = gw_time_to_coalescence(Mc, f);
    C(tau > 0, "time to coalescence positive");
    C(tau < 1000.0, "realistic time for BBH at 30 Hz");

    double N = gw_number_of_cycles(Mc, 30.0, 100.0);
    C(N > 10, "significant cycles in band");
    C(N < 10000, "realistic cycle count");
    P();
}

static void test_pn_phase(void) {
    T("post-Newtonian phase");
    double M = 60.0 * GW_MSUN;
    double eta = 0.25;

    double phase = gw_taylor_t2_phase(50.0, M, eta, 0.0, 0.0);
    /* Phase should be many radians (PN terms dominate at t_c=0) */
    C(fabs(phase) > 1.0, "significant phase");

    double f2_phase = gw_taylor_f2_phase(50.0, M, eta, 0.0, 0.0);
    C(fabs(phase - f2_phase) < 1e-10, "T2 and F2 same phase");
    P();
}

static void test_taylor_t4(void) {
    T("TaylorT4 RHS");
    double M = 60.0 * GW_MSUN;
    double eta = 0.25;
    double f = 30.0;
    double df_dt, dphi_dt;

    gw_taylor_t4_rhs(f, M, eta, &df_dt, &dphi_dt);
    C(df_dt > 0, "df/dt positive");
    C(dphi_dt > 0, "dphi/dt positive");
    C(dphi_dt > df_dt, "phase evolves faster than frequency");
    P();
}

static void test_orbital_evolution(void) {
    T("eccentric orbital evolution ODE");
    double M = 2.0 * GW_MSUN;
    double mu = 0.5 * GW_MSUN;
    double y[3] = {1.0e9, 0.6, 0.0};  /* a, e, Phi */
    double dydt[3];

    gw_orbital_evolution_rhs(dydt, y, M, mu);
    C(dydt[0] < 0, "da/dt negative (inspiral)");
    C(dydt[1] < 0, "de/dt negative (circularization)");
    C(dydt[2] > 0, "dPhi/dt positive");
    P();
}

static void test_eccentricity_functions(void) {
    T("eccentricity enhancement functions");
    double f0 = gw_eccentricity_enhancement_da(0.0);
    C(fabs(f0 - 1.0) < 1e-10, "f(0)=1");

    double f_high = gw_eccentricity_enhancement_da(0.9);
    C(f_high > 100, "strong enhancement at high e");

    double g0 = gw_eccentricity_enhancement_de(0.0);
    C(fabs(g0 - 1.0) < 1e-10, "g(0)=1");
    P();
}

static void test_reference_systems(void) {
    T("GW150914 parameter init");
    GwBinaryParams p;
    gw_params_gw150914(&p);
    C(p.m1 > 30.0*GW_MSUN, "m1 ~ 36 Msun");
    C(p.m2 > 25.0*GW_MSUN, "m2 ~ 29 Msun");
    C(p.D_L > 400.0*GW_MPC, "D_L ~ 410 Mpc");
    C(p.eta < 0.25 && p.eta > 0.2, "realistic eta");
    P();
}

int main(void) {
    printf("=== gw_binary Tests ===\n\n");
    test_chirp_mass();
    test_symmetric_mass_ratio();
    test_frequency_evolution();
    test_pn_phase();
    test_taylor_t4();
    test_orbital_evolution();
    test_eccentricity_functions();
    test_reference_systems();
    printf("\nResults: %d passed, %d failed\n", passed, failed);
    return (failed > 0) ? 1 : 0;
}
