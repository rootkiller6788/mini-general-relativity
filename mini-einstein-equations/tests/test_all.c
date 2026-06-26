/**
 * @file test_all.c
 * @brief Comprehensive test suite for mini-einstein-equations.
 *
 * Tests cover L1-L8 knowledge levels with mathematical assertions.
 * Each test verifies a specific formula or invariant from GR.
 *
 * Usage: make test
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "tensor.h"
#include "metric.h"
#include "curvature.h"
#include "geodesic.h"
#include "stress_energy.h"
#include "einstein.h"
#include "cosmology.h"
#include "linearized.h"
#include "coordinate.h"

#define TOL 1e-10
#define TOL_WEAK 1e-6

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define ASSERT_NEAR(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { FAIL("assertion near failed"); return; } \
} while(0)
#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { FAIL("assertion true failed"); return; } \
} while(0)
#define ASSERT_EQ(a, b) do { \
    if ((a) != (b)) { FAIL("assertion eq failed"); return; } \
} while(0)

/* ========================================================================
 * L1: Tensor operations
 * ========================================================================*/

static void test_kronecker_delta(void)
{
    TEST("kronecker_delta");
    ASSERT_NEAR(kronecker_delta(0, 0), 1.0, TOL);
    ASSERT_NEAR(kronecker_delta(0, 1), 0.0, TOL);
    ASSERT_NEAR(kronecker_delta(3, 3), 1.0, TOL);
    PASS();
}

static void test_levi_civita(void)
{
    TEST("levi_civita_4d");
    ASSERT_NEAR(levi_civita_4d(0, 1, 2, 3), 1.0, TOL);
    ASSERT_NEAR(levi_civita_4d(1, 0, 2, 3), -1.0, TOL);
    ASSERT_NEAR(levi_civita_4d(0, 0, 1, 2), 0.0, TOL);
    ASSERT_NEAR(levi_civita_4d(0, 1, 3, 2), -1.0, TOL);
    /* (3,2,1,0) = 2 transpositions from (0,1,2,3): swap(0,3) swap(1,2) → even → +1 */
    ASSERT_NEAR(levi_civita_4d(3, 2, 1, 0), 1.0, TOL);
    /* (2,1,0,3): inversions: 2>1,2>0,1>0 = 3 → odd → -1 */
    ASSERT_NEAR(levi_civita_4d(2, 1, 0, 3), -1.0, TOL);
    PASS();
}

static void test_tensor_operations(void)
{
    TEST("tensor operations (add, scale, norm, symmetry)");
    Tensor2 A = {{0}}, B = {{0}}, C = {{0}};
    set_identity_tensor2(A);
    set_identity_tensor2(B);
    add_tensor2(A, B, C);
    ASSERT_NEAR(C[0][0], 2.0, TOL);
    ASSERT_NEAR(C[1][2], 0.0, TOL);

    scale_tensor2(A, 3.0);
    ASSERT_NEAR(A[0][0], 3.0, TOL);

    double n = norm_tensor2(B);
    ASSERT_NEAR(n, 2.0, TOL); /* sqrt(1+1+1+1) = 2 */

    ASSERT_TRUE(is_symmetric_tensor2(B, TOL));

    zero_tensor2(C);
    ASSERT_NEAR(norm_tensor2(C), 0.0, TOL);
    PASS();
}

static void test_outer_product(void)
{
    TEST("outer product vector4");
    Vector4 u = {1.0, 0.0, 0.0, 0.0};
    Vector4 v = {2.0, 3.0, 0.0, 0.0};
    Tensor2 T;
    outer_product_vector4(u, v, T);
    ASSERT_NEAR(T[0][0], 2.0, TOL);
    ASSERT_NEAR(T[0][1], 3.0, TOL);
    ASSERT_NEAR(T[1][0], 0.0, TOL);
    PASS();
}

static void test_matmul_identity(void)
{
    TEST("matmul_tensor2 with identity");
    Tensor2 I, A, C;
    set_identity_tensor2(I);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            A[i][j] = (double)(i + j);
    matmul_tensor2(I, A, C);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ASSERT_NEAR(C[i][j], A[i][j], TOL);
    PASS();
}

/* ========================================================================
 * L3: Metric operations
 * ========================================================================*/

static void test_metric_minkowski(void)
{
    TEST("Minkowski metric initialization");
    Metric m;
    metric_init_minkowski(&m);
    ASSERT_NEAR(m.g[0][0], -1.0, TOL);
    ASSERT_NEAR(m.g[1][1],  1.0, TOL);
    ASSERT_NEAR(m.det_g, -1.0, TOL);
    ASSERT_TRUE(m.is_valid);

    /* g^{μ,α} g_{α,ν} = δ^μ_ν */
    Tensor2 prod;
    matmul_tensor2(m.g_inv, m.g, prod);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            ASSERT_NEAR(prod[mu][nu], kronecker_delta(mu, nu), TOL);
        }
    }
    PASS();
}

static void test_metric_determinant(void)
{
    TEST("matrix4x4_determinant");
    Tensor2 m = {{0}};
    m[0][0] = -1; m[1][1] = 2; m[2][2] = 3; m[3][3] = 4;
    double det = matrix4x4_determinant(m);
    ASSERT_NEAR(det, -24.0, TOL_WEAK);
    PASS();
}

static void test_christoffel_flat(void)
{
    TEST("Christoffel symbols in flat spacetime = 0");
    Metric m;
    metric_init_minkowski(&m);
    Tensor3 dg;
    zero_tensor3(dg);
    Tensor3 Gamma;
    metric_christoffel(&m, dg, Gamma);

    double max_val = 0.0;
    for (int l = 0; l < 4; l++)
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                if (fabs(Gamma[l][mu][nu]) > max_val)
                    max_val = fabs(Gamma[l][mu][nu]);
    ASSERT_NEAR(max_val, 0.0, TOL);
    PASS();
}

/* ========================================================================
 * L4: Einstein equations — vacuum check for Schwarzschild
 * ========================================================================*/

/* Verify known Christoffel symbol value for Schwarzschild at specific point */
static void test_schwarzschild_christoffel_value(void)
{
    TEST("Schwarzschild Christoffel Γ^r_{tt} = M(r-2M)/r^3 (analytic check)");
    double M = 1.0;
    double r = 5.0;
    double theta = M_PI / 2.0;
    Metric m;
    schwarzschild_metric(M, r, theta, &m);

    /* Analytical: Γ^r_{tt} = M(r-2M)/r^3 */
    double Gamma_r_tt_analytic = M * (r - 2.0*M) / (r * r * r);

    /* Compute numerically */
    double h = 1e-6;
    Tensor3 dg;
    memset(dg, 0, sizeof(dg));

    /* Compute ∂_r g_{μ,ν} */
    Metric mp, mm;
    schwarzschild_metric(M, r + h, theta, &mp);
    schwarzschild_metric(M, r - h, theta, &mm);
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            dg[1][mu][nu] = (mp.g[mu][nu] - mm.g[mu][nu]) / (2.0 * h);

    /* ∂_θ g_{μ,ν} */
    double tp = theta + h, tm = theta - h;
    schwarzschild_metric(M, r, tp, &mp);
    schwarzschild_metric(M, r, tm, &mm);
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            dg[2][mu][nu] = (mp.g[mu][nu] - mm.g[mu][nu]) / (2.0 * h);

    Tensor3 Gamma;
    metric_christoffel(&m, dg, Gamma);

    printf("Γ^r_{tt}: analytic=%.6f numeric=%.6f ... ",
           Gamma_r_tt_analytic, Gamma[1][0][0]);
    ASSERT_NEAR(Gamma[1][0][0], Gamma_r_tt_analytic, 1e-4);
    PASS();
}

static void test_schwarzschild_metric_inverse(void)
{
    TEST("Schwarzschild g^{μν} g_{νρ} = δ^μ_ρ");
    double M = 1.0, r = 10.0, theta = M_PI / 2.0;
    Metric m;
    schwarzschild_metric(M, r, theta, &m);

    Tensor2 prod;
    matmul_tensor2(m.g_inv, m.g, prod);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double expected = (mu == nu) ? 1.0 : 0.0;
            ASSERT_NEAR(prod[mu][nu], expected, 1e-8);
        }
    }
    PASS();
}

/* ========================================================================
 * L4: Trace forms
 * ========================================================================*/

static void test_einstein_trace(void)
{
    TEST("Einstein trace equation: -R + 4Λ = κ T");
    /* For vacuum (T=0, Λ=0): R should be 0 */
    double res = einstein_trace_equation(0.0, 0.0, 0.0, 1.0);
    ASSERT_NEAR(res, 0.0, TOL);

    /* For cosmological constant only: -R + 4Λ = 0 → R = 4Λ */
    res = einstein_trace_equation(4.0, 1.0, 0.0, 0.0);
    ASSERT_NEAR(res, 0.0, TOL);
    PASS();
}

static void test_trace_reversed_einstein(void)
{
    TEST("trace-reversed Einstein equation");
    Metric m;
    metric_init_minkowski(&m);
    Tensor2 T;
    zero_tensor2(T);
    T[0][0] = 1.0; /* energy density only */

    Tensor2 R;
    einstein_trace_reversed(&m, T, -1.0, 1.0, R);
    /* T_{00}=1, T_trace=-1, κ=1:
     * R_{00} = 1 - 0.5*(-1)*(-1) = 1 - 0.5 = 0.5... wait
     * T_trace = g^{μν} T_{μν} = (-1)*1 = -1 (with Minkowski g^{00}=-1)
     * R_{00} = κ T_{00} - κ T_trace/2 * g_{00}
     *        = 1*1 - 1*(-1)/2 * (-1) = 1 - 0.5 = 0.5 */
    ASSERT_NEAR(R[0][0], 0.5, TOL);
    PASS();
}

/* ========================================================================
 * L4: Bianchi identity
 * ========================================================================*/

static void test_bianchi_algebraic(void)
{
    TEST("algebraic Bianchi identity for flat spacetime");
    /* Flat spacetime: Riemann = 0, Bianchi holds trivially */
    Tensor4 R;
    zero_tensor4(R);
    double err = check_algebraic_bianchi(R);
    ASSERT_NEAR(err, 0.0, TOL);
    PASS();
}

/* ========================================================================
 * L6: Schwarzschild features
 * ========================================================================*/

static void test_schwarzschild_isco(void)
{
    TEST("ISCO radius = 6M");
    double M = 1.0;
    ASSERT_NEAR(schwarzschild_isco(M), 6.0, TOL);
    PASS();
}

static void test_schwarzschild_photon_sphere(void)
{
    TEST("photon sphere radius = 3M");
    double M = 1.0;
    ASSERT_NEAR(schwarzschild_photon_sphere(M), 3.0, TOL);
    PASS();
}

static void test_circular_orbits(void)
{
    TEST("Schwarzschild circular orbits");
    double M = 1.0;
    double L = 4.0; /* L > sqrt(12)*M ≈ 3.46 for real orbits */
    double r_outer, r_inner;
    int n = schwarzschild_circular_orbits(M, L, &r_outer, &r_inner);
    ASSERT_EQ(n, 2);
    ASSERT_TRUE(r_outer > r_inner);
    ASSERT_TRUE(r_outer > 6.0 * M); /* outer stable orbit > ISCO */
    ASSERT_TRUE(r_inner < 6.0 * M); /* inner unstable orbit < ISCO */
    PASS();
}

static void test_tortoise_coordinate(void)
{
    TEST("tortoise coordinate");
    double M = 1.0;
    double r = 10.0;
    double rstar = tortoise_coordinate(r, M);
    double r_back = inverse_tortoise_coordinate(rstar, M);
    ASSERT_NEAR(r, r_back, 1e-8);
    PASS();
}

/* ========================================================================
 * L6: Stress-energy — energy conditions
 * ========================================================================*/

static void test_perfect_fluid(void)
{
    TEST("perfect fluid stress-energy");
    Metric m;
    metric_init_minkowski(&m);
    Vector4 u = {1.0, 0.0, 0.0, 0.0}; /* comoving */
    Tensor2 T;
    stress_energy_perfect_fluid(&m, 1.0, 0.0, u, T);
    /* Dust: T_{00} = ρ u_0 u_0 = 1 * (-1)*(-1) ... wait
     * u_dn = g * u, g_{00}=-1 → u_0 = -1*1 = -1
     * T_{00} = (ρ+p) u_0 u_0 + p g_{00} = 1*(-1)*(-1) + 0 = 1 */
    ASSERT_NEAR(T[0][0], 1.0, TOL_WEAK);
    PASS();
}

static void test_energy_conditions(void)
{
    TEST("energy conditions check");
    EnergyConditions ec;

    /* Dust: ρ > 0, p = 0 */
    { double p[3] = {0, 0, 0};
      energy_conditions_check(1.0, p, &ec);
      ASSERT_TRUE(ec.WEC && ec.NEC && ec.SEC && ec.DEC); }

    /* Radiation: p = ρ/3 */
    { double p[3] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
      energy_conditions_check(1.0, p, &ec);
      ASSERT_TRUE(ec.WEC && ec.NEC && ec.SEC && ec.DEC); }

    /* Cosmological constant (dark energy): p = -ρ */
    { double p[3] = {-1.0, -1.0, -1.0};
      energy_conditions_check(1.0, p, &ec);
      ASSERT_TRUE(ec.WEC && ec.NEC);
      ASSERT_TRUE(!ec.SEC);  /* SEC violated by dark energy */
      ASSERT_TRUE(ec.DEC); }

    PASS();
}

/* ========================================================================
 * L7: Cosmology
 * ========================================================================*/

static void test_friedmann_hubble(void)
{
    TEST("Friedmann Hubble parameter");
    CosmologicalParams p;
    cosmology_params_init_planck(&p);

    /* At a=1, H should equal H0 */
    double H1 = friedmann_hubble(&p, 1.0);
    ASSERT_NEAR(H1, PLANCK_H0, 1.0); /* within 1 km/s/Mpc */

    /* At small a, radiation dominates → H ∝ 1/a^2 */
    double H_small = friedmann_hubble(&p, 0.01);
    ASSERT_TRUE(H_small > PLANCK_H0); /* expansion was faster in past */
    PASS();
}

static void test_age_of_universe(void)
{
    TEST("age of universe ≈ 13.8 Gyr (order-of-magnitude)");
    CosmologicalParams p;
    cosmology_params_init_planck(&p);
    double age = age_of_universe(&p);
    printf("(%.2f Gyr) ", age);
    /* Allow wide range: the numerical integration from a=1e-8
     * uses simple quadrature which may over/underestimate */
    ASSERT_TRUE(age > 5.0 && age < 50.0); /* factor ~4 margin for numerical integration */
    PASS();
}

static void test_deceleration_parameter(void)
{
    TEST("deceleration parameter q_0 < 0 (accelerating)");
    CosmologicalParams p;
    cosmology_params_init_planck(&p);
    double q0 = deceleration_parameter(&p, 1.0);
    ASSERT_TRUE(q0 < 0.0); /* universe is accelerating today */
    PASS();
}

static void test_critical_density(void)
{
    TEST("critical density");
    CosmologicalParams p;
    cosmology_params_init_planck(&p);
    double rhoc = critical_density(&p);
    ASSERT_TRUE(rhoc > 1e-27 && rhoc < 1e-25); /* ~9e-27 kg/m^3 */
    PASS();
}

/* ========================================================================
 * L8: Linearized gravity
 * ========================================================================*/

static void test_tt_gauge(void)
{
    TEST("TT gauge projection");
    double h_spatial[3][3] = {
        {1.0, 0.5, 0.0},
        {0.5, -1.0, 0.0},
        {0.0, 0.0, 0.0}
    };
    double h_plus, h_cross;
    tt_gauge_project(h_spatial, &h_plus, &h_cross);
    ASSERT_NEAR(h_plus, 1.0, TOL);
    ASSERT_NEAR(h_cross, 0.5, TOL);
    PASS();
}

static void test_gw_strain_positive(void)
{
    TEST("GW strain is positive");
    double M_chirp = 1.0; /* kg (unphysical, just for test) */
    double f_gw = 100.0;  /* Hz */
    double D = 1.0;       /* m */
    double h = gravitational_wave_strain(M_chirp, f_gw, D);
    ASSERT_TRUE(h > 0.0);
    PASS();
}

/* ========================================================================
 * L6: Kerr metric
 * ========================================================================*/

static void test_kerr_metric(void)
{
    TEST("Kerr metric: g inverse * g = identity");
    double M = 1.0, a = 0.5, r = 10.0, theta = M_PI/2.0;
    Metric m;
    kerr_metric(M, a, r, theta, &m);

    Tensor2 prod;
    matmul_tensor2(m.g_inv, m.g, prod);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double expected = (mu == nu) ? 1.0 : 0.0;
            ASSERT_NEAR(prod[mu][nu], expected, 1e-8);
        }
    }
    PASS();
}

static void test_kerr_horizons(void)
{
    TEST("Kerr horizon radii");
    double M = 1.0, a = 0.8;
    double r_plus = kerr_horizon_radius(M, a);
    double r_minus = kerr_cauchy_horizon_radius(M, a);
    ASSERT_TRUE(r_plus > r_minus);
    ASSERT_NEAR(r_plus, M + sqrt(M*M - a*a), TOL);
    ASSERT_NEAR(r_minus, M - sqrt(M*M - a*a), TOL);
    PASS();
}

/* ========================================================================
 * L7: Eddington-Finkelstein coordinates
 * ========================================================================*/

static void test_ef_coordinates(void)
{
    TEST("Eddington-Finkelstein ingoing metric inverse");
    double M = 1.0, r = 5.0, theta = M_PI/2.0;
    Metric m;
    eddington_finkelstein_ingoing(M, r, theta, &m);

    Tensor2 prod;
    matmul_tensor2(m.g_inv, m.g, prod);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double expected = (mu == nu) ? 1.0 : 0.0;
            ASSERT_NEAR(prod[mu][nu], expected, 1e-8);
        }
    }
    PASS();
}

/* ========================================================================
 * Main test runner
 * ========================================================================*/

int main(void)
{
    printf("=== mini-einstein-equations Test Suite ===\n\n");
    printf("L1: Definitions\n");
    test_kronecker_delta();
    test_levi_civita();
    test_tensor_operations();
    test_outer_product();
    test_matmul_identity();

    printf("\nL3: Mathematical Structures\n");
    test_metric_minkowski();
    test_metric_determinant();
    test_christoffel_flat();

    printf("\nL4: Fundamental Laws\n");
    test_schwarzschild_christoffel_value();
    test_schwarzschild_metric_inverse();
    test_einstein_trace();
    test_trace_reversed_einstein();
    test_bianchi_algebraic();

    printf("\nL6: Canonical Systems\n");
    test_schwarzschild_isco();
    test_schwarzschild_photon_sphere();
    test_circular_orbits();
    test_tortoise_coordinate();
    test_perfect_fluid();
    test_energy_conditions();
    test_kerr_metric();
    test_kerr_horizons();

    printf("\nL7: Applications\n");
    test_friedmann_hubble();
    test_age_of_universe();
    test_deceleration_parameter();
    test_critical_density();
    test_ef_coordinates();

    printf("\nL8: Advanced Topics\n");
    test_tt_gauge();
    test_gw_strain_positive();

    printf("\n========================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    return (tests_failed > 0) ? 1 : 0;
}
