/**
 * @file    test_quadrupole.c
 * @brief   Unit tests for quadrupole radiation formulas
 *
 * Tests L1-L4: quadrupole moment, strain, luminosity
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "gw_core.h"
#include "gw_quadrupole.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)

/* ---------------------------------------------------------------
 * L1: Quadrupole moment for a binary system
 * --------------------------------------------------------------- */
static void test_quadrupole_moment_binary(void) {
    TEST("quadrupole moment for equal-mass binary");

    double m[2] = {10.0, 10.0};
    double x[6] = {5.0, 0.0, 0.0,  -5.0, 0.0, 0.0};

    GwTensor3 I;
    gw_quadrupole_moment(&I, m, x, 2);

    /* For two equal masses on x-axis: I_xx = 2*10*(5^2 - r^2/3) = 20*(25 - 25/3)
     * Actually I_xx = sum m_a (x^2 - r^2/3)
     * For a=1: x=5, y=0, z=0, r^2=25 -> x^2 - 25/3 = 25 - 8.333 = 16.667
     * I_xx_1 = 10 * 16.667 = 166.67
     * For a=2: x=-5, y=0, z=0, r^2=25 -> same: 166.67
     * Total I_xx = 333.33
     */
    double expected_Ixx = 2.0 * 10.0 * (25.0 - 25.0/3.0);
    CHECK(fabs(I.xx - expected_Ixx) < 1e-6, "I_xx mismatch");

    /* I_yy = sum m_a (y^2 - r^2/3) = 2*10*(0 - 25/3) = -166.667 */
    double expected_Iyy = 2.0 * 10.0 * (0.0 - 25.0/3.0);
    CHECK(fabs(I.yy - expected_Iyy) < 1e-6, "I_yy mismatch");

    /* Trace should be zero */
    double tr = gw_tensor_trace(&I);
    CHECK(fabs(tr) < 1e-10, "Trace not zero");

    PASS();
}

/* ---------------------------------------------------------------
 * L2: Quadrupole strain computation
 * --------------------------------------------------------------- */
static void test_quadrupole_strain(void) {
    TEST("quadrupole strain computation");

    /* Create a simple Q_ddot tensor */
    GwTensor3 Q_ddot = {1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    double r = 1.0e20;  /* 10 kpc */

    GwTensor3 h;
    gw_quadrupole_strain(&h, &Q_ddot, r);

    /* h = (2G/c^4 r) * Q_ddot
     * 2G/c^4 = 2*6.6743e-11 / (3e8)^4 ~ 1.65e-44
     * h_xx = 1.65e-44 / 1e20 * 1.0 = 1.65e-64 (extremely small)
     * Actually G/c^4 = 8.26e-45 [s^2/(kg*m)]
     * 2G/c^4 = 1.65e-44
     * h_xx = 1.65e-44 * 1.0 / 1e20 = 1.65e-64
     */
    double fac = 2.0 * GW_G / (GW_C*GW_C*GW_C*GW_C * r);
    CHECK(fabs(h.xx - fac*1.0) < 1e-70, "h_xx mismatch");
    CHECK(fabs(h.yy - fac*(-1.0)) < 1e-70, "h_yy mismatch");
    CHECK(h.xx > 0 && h.yy < 0, "polarization sign");

    PASS();
}

/* ---------------------------------------------------------------
 * L2: Binary quadrupole strain
 * --------------------------------------------------------------- */
static void test_binary_quadrupole_strain(void) {
    TEST("binary quadrupole strain (face-on)");

    double M   = 60.0 * GW_MSUN;
    double mu  = 15.0 * GW_MSUN;  /* m1=30, m2=30 -> mu = 15 */
    double a   = 1.0e6;   /* orbital separation ~ 10^6 m */
    double r   = 100.0 * GW_MPC;
    double iota = 0.0;    /* face-on */
    double Phi = 0.0;

    double hp, hx;
    gw_binary_quadrupole_strain(&hp, &hx, M, mu, a, r, iota, Phi);

    /* hp should be negative, hx should be zero (sin(0) = 0) */
    CHECK(hp < 0, "hp sign");
    CHECK(fabs(hx) < 1e-50, "hx ~ 0 for Phi=0");
    CHECK(hp > -1.0, "hp not unreasonably huge");

    PASS();
}

/* ---------------------------------------------------------------
 * L4: Circular binary luminosity
 * --------------------------------------------------------------- */
static void test_binary_luminosity(void) {
    TEST("circular binary GW luminosity");

    double M   = 2.0 * GW_MSUN;
    double mu  = 0.5 * GW_MSUN;  /* equal mass 1 Msun each */
    double a   = 1.0e9;  /* wide separation */

    double L = gw_binary_luminosity_circular(M, mu, a);

    /* L_GW should be positive */
    CHECK(L > 0, "Luminosity positive");
    /* For solar-mass objects at ~1e9 m, L_GW ~ 1e23 W (not small!) */
    CHECK(L < 1e30, "Luminosity < 1e30 W");

    /* For tighter orbit, should be larger */
    double L2 = gw_binary_luminosity_circular(M, mu, a/10.0);
    CHECK(L2 > L, "Luminosity increases as orbit shrinks");
    CHECK(L2 / L > 1e4, "Roughly ~ a^{-5} scaling");

    PASS();
}

/* ---------------------------------------------------------------
 * L2: Energy density of GWs
 * --------------------------------------------------------------- */
static void test_energy_density(void) {
    TEST("Isaacson energy density");

    double hp = 1.0e-21;
    double hx = 5.0e-22;
    double f  = 100.0;
    double omega = 2.0 * GW_PI * f;

    double rho = gw_energy_density(hp, hx, omega);

    /* rho should be positive and reasonable */
    CHECK(rho > 0, "Energy density positive");
    CHECK(rho < 1.0, "Energy density < 1 J/m^3");

    PASS();
}

/* ---------------------------------------------------------------
 * Main test runner
 * --------------------------------------------------------------- */
int main(void) {
    printf("=== gw_quadrupole Tests ===\n\n");

    test_quadrupole_moment_binary();
    test_quadrupole_strain();
    test_binary_quadrupole_strain();
    test_binary_luminosity();
    test_energy_density();

    printf("\nResults: %d passed, %d failed\n", tests_passed, tests_failed);
    return (tests_failed > 0) ? 1 : 0;
}
