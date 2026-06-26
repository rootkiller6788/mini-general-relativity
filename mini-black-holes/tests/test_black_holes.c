/**
 * @file test_black_holes.c
 * @brief Comprehensive test suite for mini-black-holes module.
 *
 * Tests cover all core APIs: metrics, thermodynamics, dynamics,
 * waves, and advanced topics. Uses standard assert().
 *
 * L1-L6: Complete coverage with mathematical assertion checks.
 * L7: Application-level tests (GW150914, EHT shadow).
 * L8-L9: Advanced topic assertions.
 */

#include "black_hole_metrics.h"
#include "black_hole_thermodynamics.h"
#include "black_hole_dynamics.h"
#include "black_hole_waves.h"
#include "black_hole_advanced.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define TOL 1e-8
#define NEAR(a,b,tol) (fabs((a)-(b)) < (tol) * (1.0 + fabs(b)))

static int tests_run = 0;
static int tests_failed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s ... ", name); } while(0)
#define CHECK(cond, msg) do { \
    if (!(cond)) { tests_failed++; printf("FAIL: %s\n", msg); } \
    else { printf("."); } \
} while(0)
#define ENDTEST() printf(" DONE\n")

/* ================================================================
 *  L1: METRIC TESTS
 * ================================================================ */

static void test_schwarzschild_init(void) {
    TEST("schwarzschild_init");
    SchwarzschildBH bh;
    schwarzschild_init(&bh, BH_SOLAR_MASS);
    CHECK(bh.M > 0, "M must be positive");
    CHECK(bh.r_s > 0, "r_s must be positive");
    CHECK(NEAR(bh.r_s, 2.0 * bh.M, TOL), "r_s = 2M in geometric units");
}

static void test_schwarzschild_metric_components(void) {
    TEST("schwarzschild_metric_components");
    SchwarzschildBH bh;
    schwarzschild_init(&bh, BH_SOLAR_MASS);

    BLPosition pos = {0, 10.0 * bh.r_s, M_PI/2, 0}; /* r = 10 r_s, equator */

    Metric4D metric;
    schwarzschild_metric(&bh, &pos, &metric);

    /* Check key components */
    double rs = bh.r_s;
    double r = pos.r;
    double expected_gtt = -(1.0 - rs / r);
    double expected_grr = 1.0 / (1.0 - rs / r);
    double expected_gthth = r * r;
    double expected_gphph = r * r; /* sin²(π/2) = 1 */

    CHECK(NEAR(metric.g[0][0], expected_gtt, TOL), "g_tt");
    CHECK(NEAR(metric.g[1][1], expected_grr, TOL), "g_rr");
    CHECK(NEAR(metric.g[2][2], expected_gthth, TOL), "g_θθ");
    CHECK(NEAR(metric.g[3][3], expected_gphph, TOL), "g_φφ");

    /* Off-diagonal must be zero */
    CHECK(NEAR(metric.g[0][1], 0.0, TOL), "g_tr = 0");
    CHECK(NEAR(metric.g[0][3], 0.0, TOL), "g_tφ = 0");
}

static void test_metric_inverse(void) {
    TEST("metric_inverse_identity");
    SchwarzschildBH bh;
    schwarzschild_init(&bh, BH_SOLAR_MASS);
    BLPosition pos = {0, 5.0 * bh.r_s, M_PI/3, 0.5};
    Metric4D m;
    schwarzschild_metric(&bh, &pos, &m);

    double max_dev = metric_inverse_check(&m);
    CHECK(max_dev < 1e-10, "g·g^{-1} = I");
}

static void test_kerr_metric(void) {
    TEST("kerr_metric");
    KerrBH bh;
    kerr_init(&bh, 10.0 * BH_SOLAR_MASS, 0.5 * BH_SOLAR_MASS
              * BH_SOLAR_MASS * BH_SOLAR_MASS * BH_G / BH_C);

    BLPosition pos = {0, 6.0 * bh.M, M_PI/2, 0};
    Metric4D m;
    kerr_metric(&bh, &pos, &m);

    /* g_tφ should be non-zero (frame dragging) */
    CHECK(fabs(m.g[0][3]) > 1e-30, "g_tφ non-zero (frame dragging)");

    /* g_tφ = g_φt (symmetry) */
    CHECK(NEAR(m.g[0][3], m.g[3][0], TOL), "g_tφ = g_φt");
}

static void test_rn_metric(void) {
    TEST("reissner_nordstrom_metric");
    ReissnerNordstromBH bh;
    reissner_nordstrom_init(&bh, 10.0 * BH_SOLAR_MASS, 1e20);

    /* For solar-scale mass, Q_geom << M_geom, so near-Schwarzschild */
    CHECK(bh.r_plus > 0, "Outer horizon exists");
    CHECK(bh.r_plus >= bh.r_minus, "Outer ≥ inner horizon");
}

static void test_kn_metric(void) {
    TEST("kerr_newman_metric");
    KerrNewmanBH bh;
    kerr_newman_init(&bh, 10.0 * BH_SOLAR_MASS,
                     0.3 * BH_SOLAR_MASS * BH_SOLAR_MASS
                     * BH_SOLAR_MASS * BH_G / BH_C,
                     1e20);

    BLPosition pos = {0, 8.0 * bh.M, M_PI/3, 0};
    Metric4D m;
    kerr_newman_metric(&bh, &pos, &m);
    CHECK(m.g[0][0] < 0, "g_tt negative (timelike at infinity)");
}

/* ================================================================
 *  L2-L4: THERMODYNAMICS TESTS
 * ================================================================ */

static void test_hawking_temperature(void) {
    TEST("hawking_temperature");
    /* Solar mass BH: T ~ 6e-8 K */
    double T = hawking_temperature_schwarzschild(BH_SOLAR_MASS);
    CHECK(T > 0, "Temperature positive");
    CHECK(T < 1e-6, "T_sun < 1 μK (very cold!)");

    /* T ∝ 1/M: compare 2M vs M */
    double T2 = hawking_temperature_schwarzschild(2.0 * BH_SOLAR_MASS);
    CHECK(NEAR(T2, T / 2.0, 0.01), "T(2M) ≈ T(M)/2");
}

static void test_bekenstein_hawking_entropy(void) {
    TEST("bekenstein_hawking_entropy");
    double A = schwarzschild_horizon_area(BH_SOLAR_MASS);
    double S = bekenstein_hawking_entropy(A);
    CHECK(S > 0, "Entropy positive");
    /* S ∝ M² */
    double A2 = schwarzschild_horizon_area(2.0 * BH_SOLAR_MASS);
    double S2 = bekenstein_hawking_entropy(A2);
    CHECK(NEAR(S2, 4.0 * S, 0.01), "S(2M) = 4·S(M)");
}

static void test_first_law(void) {
    TEST("first_law_schwarzschild");
    double M_geom = mass_to_geometric(BH_SOLAR_MASS);
    double residual = first_law_schwarzschild_check(M_geom, 100.0);
    CHECK(residual < 1e-10, "dM = T_H dS");
}

static void test_area_theorem(void) {
    TEST("area_increase_theorem");
    double A1 = 100.0, A2 = 150.0;
    double deltaA = area_increase_theorem(A1, A2);
    CHECK(deltaA >= 0, "δA ≥ 0 (area increase)");

    /* A_final < A_initial should give negative delta */
    deltaA = area_increase_theorem(A2, A1);
    CHECK(deltaA < 0, "δA < 0 for decreasing area (physically forbidden)");
}

static void test_heat_capacity(void) {
    TEST("schwarzschild_heat_capacity");
    double C = schwarzschild_heat_capacity(BH_SOLAR_MASS);
    CHECK(C < 0, "Heat capacity negative (Schwarzschild unstable)");
}

static void test_evaporation(void) {
    TEST("evaporation_time");
    double tau_sun = evaporation_time(BH_SOLAR_MASS);
    CHECK(tau_sun > 1e60, "Solar mass BH lifetime >> age of universe");

    /* Lighter BH evaporates faster: τ ∝ M³ */
    double tau_half = evaporation_time(0.5 * BH_SOLAR_MASS);
    CHECK(NEAR(tau_half, tau_sun / 8.0, 0.01), "τ ∝ M³");
}

/* ================================================================
 *  L5-L6: DYNAMICS TESTS
 * ================================================================ */

static void test_orbital_quantities(void) {
    TEST("orbital_quantities");
    double M_schw = mass_to_geometric(BH_SOLAR_MASS);
    double r_isco_schw = schwarzschild_isco(M_schw);
    CHECK(NEAR(r_isco_schw, 6.0 * M_schw, TOL), "r_ISCO = 6M (Schwarzschild)");

    double r_ph = photon_sphere_radius(M_schw);
    CHECK(NEAR(r_ph, 3.0 * M_schw, TOL), "r_ph = 3M");

    /* Kerr ISCO: prograde should be < 6M for a > 0 */
    double r_isco_kerr = kerr_isco_radius(M_schw, 0.5 * M_schw, 1);
    CHECK(r_isco_kerr < r_isco_schw, "Kerr prograde ISCO < Schwarzschild ISCO");

    /* Kerr retrograde ISCO > 6M */
    double r_isco_retro = kerr_isco_radius(M_schw, 0.5 * M_schw, 0);
    CHECK(r_isco_retro > r_isco_schw, "Kerr retrograde ISCO > Schwarzschild ISCO");
}

static void test_effective_potential(void) {
    TEST("effective_potential");
    double M_schw = mass_to_geometric(BH_SOLAR_MASS);
    double rs = 2.0 * M_schw;

    /* At ISCO (r=6M), L²=12M² gives V_eff extremum */
    double L_sq_isco = 12.0 * M_schw * M_schw;
    double V_eff = schwarzschild_effective_potential(rs, 6.0 * M_schw, L_sq_isco);

    /* Derivative should vanish at ISCO */
    double dV = schwarzschild_eff_potential_deriv(rs, 6.0 * M_schw, L_sq_isco);
    CHECK(fabs(dV) < 0.1, "dV_eff/dr ≈ 0 at ISCO");

    /* V_eff² < 1 for bound orbits */
    CHECK(V_eff < 1.0, "Bound orbit: V_eff < 1");
}

static void test_kretschmann_scalar(void) {
    TEST("kretschmann_scalar");
    double M = mass_to_geometric(BH_SOLAR_MASS);

    /* At horizon r = 2M: K = 3/(4M⁴) — finite */
    double K_horizon = schwarzschild_kretschmann(M, 2.0 * M);
    double K_expected = 3.0 / (4.0 * M * M * M * M);
    CHECK(NEAR(K_horizon, K_expected, TOL), "K(r=2M) = 3/(4M⁴)");

    /* K ∝ 1/r⁶: check ratio */
    double K_r2 = schwarzschild_kretschmann(M, 4.0 * M);
    /* K(4M)/K(2M) = (2M/4M)⁶ = (1/2)⁶ = 1/64 */
    CHECK(NEAR(K_r2 / K_horizon, 1.0 / 64.0, 0.01), "K ∝ 1/r⁶");
}

static void test_penrose_process(void) {
    TEST("penrose_process");
    double M = mass_to_geometric(BH_SOLAR_MASS);
    double eff = penrose_process_efficiency(M, 0.9 * M);
    CHECK(eff > 0.0, "Penrose efficiency > 0 for rotating BH");
    CHECK(eff < 0.3, "Penrose efficiency < 30%");

    /* Schwarzschild: a=0 → no Penrose process */
    double eff_schw = penrose_process_efficiency(M, 0.0);
    CHECK(NEAR(eff_schw, 0.0, TOL), "No Penrose process for a=0");
}

static void test_superradiance(void) {
    TEST("superradiance_condition");
    double M = mass_to_geometric(BH_SOLAR_MASS);
    double a = 0.9 * M;
    double Omega_H = a / (2.0 * M * (M + sqrt(M*M - a*a)));

    /* Superradiant: ω < m Ω_H */
    int sr = superradiance_condition(0.5 * Omega_H, 1, M, a);
    CHECK(sr == 1, "ω < Ω_H → superradiant");

    /* Non-superradiant: ω > m Ω_H */
    sr = superradiance_condition(2.0 * Omega_H, 1, M, a);
    CHECK(sr == 0, "ω > Ω_H → not superradiant");
}

/* ================================================================
 *  L7: APPLICATION TESTS
 * ================================================================ */

static void test_gw150914_parameters(void) {
    TEST("GW150914_parameters");
    double m1 = 36.0 * BH_SOLAR_MASS;
    double m2 = 29.0 * BH_SOLAR_MASS;
    double M_chirp = chirp_mass(m1, m2);
    CHECK(M_chirp > 0, "Chirp mass > 0");
    CHECK(NEAR(M_chirp / BH_SOLAR_MASS, 28.1, 0.5), "M_chirp ≈ 28 M_sun");

    (void)final_mass_after_merger(m1, m2, 0.3, 0.3);
    double chi_final = final_spin_after_merger(m1, m2, 0.3, 0.3);
    CHECK(chi_final > 0.5, "Final spin > 0.5");
    CHECK(chi_final < 1.0, "Final spin ≤ 1.0");
}

static void test_qnm_frequencies(void) {
    TEST("qnm_frequencies");
    QuasinormalMode qnm = schwarzschild_qnm_fundamental(2, 2);

    /* M ω_220 = 0.74734 - 0.17792 i */
    CHECK(NEAR(qnm.omega_R, 0.74734, 1e-4), "ω_R check");
    CHECK(NEAR(qnm.omega_I, -0.17792, 1e-4), "ω_I check");

    double f_hz = qnm_frequency_hz(&qnm, 65.0 * BH_SOLAR_MASS);
    CHECK(f_hz > 100.0, "f > 100 Hz");
    CHECK(f_hz < 500.0, "f < 500 Hz");
}

static void test_area_theorem_gw150914(void) {
    TEST("area_theorem_gw150914");
    /* GW150914 parameters (geometric units) */
    double M1 = mass_to_geometric(36.0 * BH_SOLAR_MASS);
    double M2 = mass_to_geometric(29.0 * BH_SOLAR_MASS);
    double Mf = mass_to_geometric(62.0 * BH_SOLAR_MASS);

    int holds = area_theorem_merger_check(M1, 0.3*M1, 0,
                                            M2, 0.3*M2, 0,
                                            Mf, 0.68*Mf, 0);
    CHECK(holds == 1, "Area theorem holds for GW150914");

    double energy_bound = hawking_area_bound_energy(M1, M2);
    double radiated = (M1 + M2 - Mf);
    CHECK(radiated <= energy_bound + 1e-10, "Radiated energy ≤ Hawking bound");
}

static void test_ringdown_waveform(void) {
    TEST("ringdown_waveform");
    QuasinormalMode qnm = schwarzschild_qnm_fundamental(2, 2);
    double f = qnm_frequency_hz(&qnm, 65.0 * BH_SOLAR_MASS);
    double tau = qnm_damping_time(&qnm, 65.0 * BH_SOLAR_MASS);

    /* At t=0: h = A cos(φ₀). Choose φ₀=0 → h(0)=A */
    double h0 = ringdown_waveform(0.0, 1.0, f, tau, 0.0);
    CHECK(NEAR(h0, 1.0, TOL), "h(0) = A cos(0) = A");

    /* At t=τ with φ₀ = -2πfτ: cos(2πfτ + φ₀) = cos(0) = 1
     * → h(τ) = A · e^{-1} · 1 */
    double phi_comp = -2.0 * M_PI * f * tau;
    double h_tau = ringdown_waveform(tau, 1.0, f, tau, phi_comp);
    CHECK(NEAR(fabs(h_tau), exp(-1.0), 0.01), "h(τ) = A/e (phase compensated)");
}

/* ================================================================
 *  L8-L9: ADVANCED TOPICS TESTS
 * ================================================================ */

static void test_no_hair(void) {
    TEST("no_hair_theorem");
    double M = mass_to_geometric(BH_SOLAR_MASS);
    double a = 0.5 * M;
    KerrMultipoles mp;
    kerr_multipoles(M, a, &mp);

    /* M₂/(J²/M) = -1 for Kerr */
    double ratio = no_hair_test_quadrupole(mp.M2, M, mp.J);
    CHECK(NEAR(ratio, -1.0, TOL), "M₂/(J²/M) = -1 (Kerr no-hair)");
}

static void test_holographic_principle(void) {
    TEST("holographic_bound");
    double A_schw = schwarzschild_horizon_area(BH_SOLAR_MASS);
    double N_bits = holographic_bound_bits(A_schw);
    CHECK(N_bits > 1e70, "BH contains enormous information in bits");

    /* BH saturates the holographic bound: S = k_B A/(4 l_P²) */
    double S = bekenstein_hawking_entropy(A_schw);
    double margin = holographic_bound_margin(S, A_schw);
    /* Because S and A involve huge numbers (~10^54), the absolute
     * margin may suffer floating-point cancellation. Check ratio. */
    double rel_margin = (S > 0) ? fabs(margin / S) : fabs(margin);
    CHECK(rel_margin < 1e-10, "BH saturates holographic bound (rel margin)");
}

static void test_information_paradox(void) {
    TEST("page_curve");
    double M0 = 1e12; /* 10^12 kg primordial BH */
    double t_page = page_time(M0);

    /* Before Page time: S_entanglement grows */
    double S_early = page_entanglement_entropy(M0, 0.1 * t_page);
    double S_mid = page_entanglement_entropy(M0, t_page);
    CHECK(S_mid > S_early, "Entanglement entropy grows to peak at t_Page");

    /* After Page time: S_entanglement decreases */
    double S_late = page_entanglement_entropy(M0, 1.5 * t_page);
    if (S_late > 0) {
        CHECK(S_late <= S_mid, "Entanglement entropy decreases after t_Page");
    }
}

static void test_cosmic_censorship(void) {
    TEST("cosmic_censorship");
    double M = mass_to_geometric(BH_SOLAR_MASS);
    CHECK(cosmic_censorship_check(M, 0.5*M, 0.0) == 1, "Kerr a<M → BH");
    CHECK(cosmic_censorship_check(M, M, 0.0) == 1, "Extremal Kerr → BH");
    CHECK(cosmic_censorship_check(M, 1.5*M, 0.0) == 0, "a>M → naked singularity");
}

static void test_ads_cft(void) {
    TEST("ads_cft_entropy");
    /* N=4 SYM at strong coupling: S ∝ N² V T³ */
    double S_ads = ads_cft_entropy(3.0, 1e-60, 1e30);
    CHECK(S_ads > 0, "AdS/CFT entropy positive");
}

/* ================================================================
 *  UTILITY AND EDGE CASE TESTS
 * ================================================================ */

static void test_unit_conversions(void) {
    TEST("unit_conversions");
    double M_kg = BH_SOLAR_MASS;
    double M_geom = mass_to_geometric(M_kg);
    double M_back = geometric_to_mass(M_geom);
    CHECK(NEAR(M_back, M_kg, 1e-6), "Mass conversion round-trip");

    double J_si = BH_SOLAR_MASS * 1e4;
    double a_geom = spin_to_geometric(J_si, BH_SOLAR_MASS);
    CHECK(a_geom >= 0, "Spin parameter non-negative");
}

static void test_proper_distance(void) {
    TEST("proper_distance");
    double rs = schwarzschild_radius(BH_SOLAR_MASS);
    double d_proper = schwarzschild_proper_distance(rs, 10.0*rs, 20.0*rs);
    CHECK(d_proper > (20.0 * rs - 10.0 * rs), "Proper distance > coordinate distance");
}

static void test_redshift(void) {
    TEST("gravitational_redshift");
    double rs = schwarzschild_radius(BH_SOLAR_MASS);
    double r_emit = 10.0 * rs;
    double z = schwarzschild_redshift(rs, r_emit);
    CHECK(z > 0, "Redshift positive");
    CHECK(z < 1.0, "z < 1 at r=10 r_s");

    /* As r → r_s, z → ∞ */
    double z_near = schwarzschild_redshift(rs, 1.01 * rs);
    CHECK(z_near > 5.0, "z >> 1 near horizon");
}

static void test_planck_units(void) {
    TEST("planck_units");
    double mP = planck_mass();
    double lP = planck_length();
    double tP = planck_time();

    /* tP = lP / c */
    CHECK(NEAR(tP, lP / BH_C, 1e-10), "t_P = l_P / c");
    /* mP lP = ħ / c */
    CHECK(NEAR(mP * lP, BH_HBAR / BH_C, 1e-10), "m_P · l_P = ħ/c");
}

/* ================================================================
 *  TEST RUNNER
 * ================================================================ */

int main(void) {
    printf("\n");
    printf("========================================================\n");
    printf("  MINI-BLACK-HOLES — Comprehensive Test Suite\n");
    printf("========================================================\n\n");

    printf("--- L1: Metric Definitions ---\n");
    test_schwarzschild_init();
    test_schwarzschild_metric_components();
    test_metric_inverse();
    test_kerr_metric();
    test_rn_metric();
    test_kn_metric();
    test_unit_conversions();

    printf("\n--- L2-L4: Thermodynamics ---\n");
    test_hawking_temperature();
    test_bekenstein_hawking_entropy();
    test_first_law();
    test_area_theorem();
    test_heat_capacity();
    test_evaporation();

    printf("\n--- L5-L6: Dynamics & Curvature ---\n");
    test_orbital_quantities();
    test_effective_potential();
    test_kretschmann_scalar();
    test_penrose_process();
    test_superradiance();
    test_proper_distance();
    test_redshift();

    printf("\n--- L7: Applications ---\n");
    test_gw150914_parameters();
    test_qnm_frequencies();
    test_area_theorem_gw150914();
    test_ringdown_waveform();

    printf("\n--- L8-L9: Advanced Topics ---\n");
    test_no_hair();
    test_holographic_principle();
    test_information_paradox();
    test_cosmic_censorship();
    test_ads_cft();
    test_planck_units();

    printf("\n========================================================\n");
    printf("  RESULTS: %d tests, %d failed\n", tests_run, tests_failed);
    printf("========================================================\n");

    if (tests_failed == 0) {
        printf("  ALL TESTS PASSED ✓\n\n");
        return 0;
    } else {
        printf("  %d TEST(S) FAILED ✗\n\n", tests_failed);
        return 1;
    }
}
