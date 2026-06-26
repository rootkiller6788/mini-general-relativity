/** tests/test_schwarzschild.c - Tests for mini-schwarzschild module */
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "schwarzschild_defs.h"
#include "schwarzschild_tensors.h"
#include "schwarzschild_geodesics.h"
#include "schwarzschild_observables.h"
#include "schwarzschild_numerical.h"
#include "schwarzschild_extensions.h"

#define EPS 1e-10
#define M_SUN 1.98847e30

static int tests_run = 0, tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)

int main(void) {
    printf("=== mini-schwarzschild Test Suite ===\n\n");

    printf("L1 - Definitions\n");
    TEST("schwarzschild_radius"); {
        double rs = schwarzschild_radius(M_SUN);
        double expected = 2.0 * SCHW_G_N * M_SUN / (SCHW_C * SCHW_C);
        CHECK(fabs(rs - expected) < EPS, "wrong rs");
    }
    TEST("geometric_mass"); {
        double m = geometric_mass(M_SUN);
        double expected = SCHW_G_N * M_SUN / (SCHW_C * SCHW_C);
        CHECK(fabs(m - expected) < EPS, "wrong geometric mass");
    }
    TEST("black_hole_init"); {
        SchwarzschildBlackHole bh;
        schwarzschild_black_hole_init(M_SUN, &bh);
        CHECK(bh.mass == M_SUN, "wrong mass");
        CHECK(bh.rs > 0.0, "rs not positive");
        CHECK(bh.horizon_area > 0.0, "horizon_area not positive");
    }
    TEST("schwarzschild_g_tt"); {
        double r = 2.0 * schwarzschild_radius(M_SUN);
        double gtt = schwarzschild_g_tt(r, schwarzschild_radius(M_SUN));
        CHECK(gtt < 0.0, "g_tt negative outside horizon");
        CHECK(fabs(gtt + 0.5) < EPS, "g_tt at r=2rs should be -0.5");
    }
    TEST("schwarzschild_g_rr"); {
        double rs = schwarzschild_radius(M_SUN);
        double r = 2.0 * rs;
        double grr = schwarzschild_g_rr(r, rs);
        CHECK(fabs(grr - 2.0) < EPS, "g_rr at r=2rs should be 2.0");
    }
    TEST("is_inside_horizon"); {
        double rs = 1000.0;
        CHECK(schwarzschild_is_inside_horizon(500.0, rs) == 1, "should be inside");
        CHECK(schwarzschild_is_inside_horizon(2000.0, rs) == 0, "should be outside");
    }

    printf("\nL2 - Core Concepts\n");
    TEST("time_dilation_factor"); {
        double rs = 1000.0;
        double f = schwarzschild_time_dilation_factor(10000.0, rs);
        CHECK(f > 0.0 && f < 1.0, "dilation in (0,1)");
    }
    TEST("gravitational_redshift"); {
        double rs = 1000.0;
        double z = schwarzschild_redshift_at_infinity(2000.0, rs);
        CHECK(z > 0.0, "redshift positive");
    }

    printf("\nL3 - Math Structures\n");
    TEST("christoffel_symbols"); {
        SchwarzschildPoint p = {0.0, 10000.0, M_PI/2.0, 0.0};
        double rs = 1000.0;
        SchwarzschildChristoffel chris;
        schwarzschild_christoffel(&p, rs, &chris);
        CHECK(fabs(chris.Gamma_thetartheta - 1.0/10000.0) < EPS, "Gamma theta_{r theta}=1/r");
    }
    TEST("riemann_tensor"); {
        double r = 10000.0, rs = 1000.0;
        double Rtrtr = schwarzschild_riemann_component(0, 1, 0, 1, r, M_PI/2.0, rs);
        CHECK(fabs(Rtrtr + rs/(r*r*r)) < EPS, "R_{trtr} = -rs/r^3");
    }
    TEST("ricci_vanishes"); {
        SchwarzschildPoint p = {0.0, 10000.0, M_PI/2.0, 0.0};
        SchwarzschildSymmetric4x4 ricci;
        schwarzschild_ricci(&p, 1000.0, &ricci);
        CHECK(fabs(ricci.g00)<EPS && fabs(ricci.g11)<EPS && fabs(ricci.g22)<EPS && fabs(ricci.g33)<EPS, "Ricci=0");
    }
    TEST("kretschmann_scalar"); {
        double K = schwarzschild_kretschmann(10000.0, 1000.0);
        CHECK(K > 0.0, "Kretschmann positive");
    }
    TEST("kretschmann_finite_at_horizon"); {
        double K = schwarzschild_kretschmann(1000.0, 1000.0);
        CHECK(isfinite(K), "Kretschmann finite at horizon");
    }

    printf("\nL4 - Fundamental Laws\n");
    TEST("einstein_vanishes"); {
        SchwarzschildPoint p = {0.0, 10000.0, M_PI/2.0, 0.0};
        SchwarzschildSymmetric4x4 einstein;
        schwarzschild_einstein(&p, 1000.0, &einstein);
        double maxc = fmax(fabs(einstein.g00), fmax(fabs(einstein.g11), fmax(fabs(einstein.g22), fabs(einstein.g33))));
        CHECK(maxc < EPS, "Einstein=0 in vacuum");
    }

    printf("\nL5 - Computational Methods\n");
    TEST("rk4_geodesic_step"); {
        double y[8] = {0.0, 100.0, M_PI/2.0, 0.0, -1.0, -0.01, 0.0, 3.0};
        double y_out[8];
        schwarzschild_rk4_step_geodesic(0.0, y, 1.0, y_out, 10.0);
        CHECK(isfinite(y_out[0]) && isfinite(y_out[1]), "RK4 valid");
    }
    TEST("adaptive_rk45_step"); {
        double y[8] = {0.0, 100.0, M_PI/2.0, 0.0, -1.0, -0.01, 0.0, 3.0};
        double lambda = 0.0, h = 1.0;
        int status = schwarzschild_rk45_adaptive_step(&lambda, y, &h, 1e-8, 10.0);
        CHECK(status >= -1 && status <= 1, "RK45 valid status");
    }

    printf("\nL6 - Canonical Systems\n");
    TEST("photon_sphere = 1.5*rs"); {
        CHECK(fabs(schwarzschild_photon_sphere_radius(1000.0) - 1500.0) < EPS, "r_ph=1.5rs");
    }
    TEST("isco = 3*rs"); {
        CHECK(fabs(schwarzschild_isco_radius(1000.0) - 3000.0) < EPS, "r_isco=3rs");
    }
    TEST("isco_energy"); {
        CHECK(fabs(schwarzschild_E_isco() - sqrt(8.0/9.0)) < EPS, "E_isco=sqrt(8/9)");
    }
    TEST("isco_angular_momentum"); {
        double m = 500.0;
        CHECK(fabs(schwarzschild_L_isco(2.0*m) - sqrt(12.0)*m) < EPS, "L_isco=sqrt(12)*m");
    }
    TEST("perihelion_precession"); {
        double dp = schwarzschild_perihelion_precession_per_orbit(5.79e10, 0.2056, M_SUN);
        CHECK(dp > 0.0, "precession positive");
        double arcsec_per_century = dp * (100.0/0.2409) * 180.0/M_PI * 3600.0;
        CHECK(fabs(arcsec_per_century - 42.98) < 0.5, "Mercury ~43 arcsec/century");
    }
    TEST("light_deflection"); {
        double alpha = schwarzschild_light_deflection_angle(6.957e8, M_SUN);
        double arcsec = alpha * 180.0/M_PI * 3600.0;
        CHECK(fabs(arcsec - 1.75) < 0.01, "Solar deflection ~1.75 arcsec");
    }
    TEST("shapiro_delay_example"); {
        double dt = schwarzschild_shapiro_earth_venus(1.5e11, 0.723*1.5e11, M_SUN);
        CHECK(dt > 0.0, "Shapiro delay positive");
    }

    printf("\nL7 - Applications\n");
    TEST("gps_time_dilation"); {
        double df = schwarzschild_gps_time_dilation(2.66e7, 6.371e6, 5.972e24, 3874.0);
        CHECK(df > 0.0, "GPS GR+SR net positive");
    }
    TEST("black_hole_shadow"); {
        double M = 6.5e9 * M_SUN;
        double D = 16.8e6 * SCHW_PC;
        double theta = schwarzschild_shadow_angular_radius(M, D);
        CHECK(theta > 0.0, "shadow radius positive");
    }
    TEST("accretion_efficiency"); {
        double eta = schwarzschild_accretion_efficiency();
        CHECK(fabs(eta - (1.0 - sqrt(8.0/9.0))) < EPS, "eta = 1-sqrt(8/9)");
    }
    TEST("accretion_disk_flux_outside_isco"); {
        double rs = schwarzschild_radius(10.0 * M_SUN);
        double F = schwarzschild_disk_flux(4.0 * rs, 10.0*M_SUN, 1e15);
        CHECK(F > 0.0, "flux positive outside ISCO");
    }
    TEST("accretion_disk_flux_inside_isco"); {
        double rs = schwarzschild_radius(10.0 * M_SUN);
        double F = schwarzschild_disk_flux(1.0 * rs, 10.0*M_SUN, 1e15);
        CHECK(F == 0.0, "flux zero inside ISCO");
    }

    printf("\nL8 - Advanced Topics\n");
    TEST("hawking_temperature"); {
        double T = schwarzschild_hawking_temperature(M_SUN);
        CHECK(T > 0.0 && T < 1e-6, "T_H ~ 6e-8 K for M_sun");
    }
    TEST("bekenstein_hawking_entropy"); {
        double S = schwarzschild_bekenstein_hawking_entropy(M_SUN);
        CHECK(S > 0.0, "BH entropy positive");
    }
    TEST("evaporation_time"); {
        double t = schwarzschild_evaporation_time(M_SUN);
        CHECK(t > 1e60, "t_evap >> age of universe for M_sun");
    }
    TEST("hawking_luminosity"); {
        double L = schwarzschild_hawking_luminosity(1e12);
        CHECK(L > 0.0, "Hawking luminosity positive for small BH");
    }
    TEST("rn_black_hole_init"); {
        ReissnerNordstromBH rn;
        rn_black_hole_init(M_SUN, 1.0e20, &rn);
        CHECK(!rn.is_naked_singularity, "small charge, not naked singularity");
    }
    TEST("rn_extremal_charge"); {
        double Q_max = rn_extremal_charge(M_SUN);
        CHECK(Q_max > 0.0, "extremal charge positive");
    }
    TEST("kerr_black_hole_init"); {
        KerrBH kerr;
        double a = 0.5 * SCHW_C * M_SUN * geometric_mass(M_SUN);
        kerr_black_hole_init(M_SUN, a, &kerr);
        CHECK(kerr.r_plus > 0.0, "Kerr outer horizon exists");
        CHECK(kerr.isco_prograde < kerr.isco_retrograde, "prograde ISCO < retrograde");
    }
    TEST("kerr_frame_dragging"); {
        KerrBH kerr;
        double a = 0.5 * SCHW_C * M_SUN * geometric_mass(M_SUN);
        kerr_black_hole_init(M_SUN, a, &kerr);
        double omega = kerr_frame_dragging_omega(2.0*kerr.r_plus, M_PI/2.0, &kerr);
        CHECK(omega > 0.0, "frame dragging positive");
    }
    TEST("kerr_ergosphere"); {
        KerrBH kerr;
        double a = 0.5 * SCHW_C * M_SUN * geometric_mass(M_SUN);
        kerr_black_hole_init(M_SUN, a, &kerr);
        double r_ergo = kerr_ergosphere_radius(M_PI/2.0, &kerr);
        CHECK(r_ergo > kerr.r_plus, "ergosphere outside horizon at equator");
    }
    TEST("penrose_efficiency_range"); {
        double a_spin = 0.5 * geometric_mass(M_SUN) * SCHW_C * M_SUN;
        double eff = kerr_penrose_max_efficiency(a_spin);
        CHECK(eff > 0.0 && eff < 0.3, "Penrose efficiency 0-30 percent");
    }
    TEST("qnm_fundamental"); {
        double freq, tau;
        schwarzschild_qnm_fundamental_l2_n0(30.0*M_SUN, &freq, &tau);
        CHECK(freq > 0.0, "QNM frequency positive");
        CHECK(tau > 0.0, "QNM damping time positive");
    }
    TEST("qnm_quality_factor"); {
        double Q = schwarzschild_qnm_quality_factor(M_SUN);
        CHECK(Q > 0.0, "Q factor positive");
    }

    printf("\nSafety - Edge Cases\n");
    TEST("zero_mass"); {
        CHECK(schwarzschild_radius(0.0) == 0.0, "zero mass -> zero rs");
        CHECK(schwarzschild_hawking_temperature(0.0) == 0.0, "zero mass -> zero T_H");
    }
    TEST("negative_mass"); {
        CHECK(schwarzschild_radius(-1.0) == 0.0, "negative mass -> zero rs");
    }
    TEST("null_pointer_defense"); {
        schwarzschild_metric_covariant(NULL, 1.0, NULL);
        schwarzschild_black_hole_init(1.0, NULL);
        schwarzschild_riemann(NULL, 1.0, NULL);
        schwarzschild_trace_null_geodesic(10.0, 1.0, 100.0, 1e-8, 1000.0, NULL);
        rn_black_hole_init(1.0, 1.0, NULL);
        kerr_black_hole_init(1.0, 1.0, NULL);
        PASS();
    }

    printf("\n=== Results: %d checks passed in %d tests ===\n", tests_passed, tests_run);
    return 0;
}
