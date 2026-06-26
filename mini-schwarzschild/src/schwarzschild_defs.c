/**
 * schwarzschild_defs.c — Core definitions and metric computation
 *
 * Reference: Wald (1984) Ch.6; Carroll (2004) Ch.5; MTW (1973) Ch.23
 * MIT 8.962 — General Relativity
 *
 * Knowledge coverage:
 *   L1 — Definitions: metric components, Schwarzschild radius, coordinate types
 *   L2 — Core Concepts: event horizon, coordinate singularities, tortoise coordinate
 */

#include "schwarzschild_defs.h"
#include <math.h>
#include <stdlib.h>

/* ==========================================================================
 * Schwarzschild radius and derived properties
 * ========================================================================== */

double schwarzschild_radius(double mass) {
    if (mass <= 0.0) return 0.0;
    return 2.0 * SCHW_G_N * mass / (SCHW_C * SCHW_C);
}

double geometric_mass(double mass) {
    if (mass <= 0.0) return 0.0;
    return SCHW_G_N * mass / (SCHW_C * SCHW_C);
}

double schwarzschild_surface_gravity(double mass) {
    if (mass <= 0.0) return 0.0;
    return SCHW_C * SCHW_C * SCHW_C * SCHW_C / (4.0 * SCHW_G_N * mass);
}

double schwarzschild_horizon_area(double mass) {
    double rs = schwarzschild_radius(mass);
    return 4.0 * M_PI * rs * rs;
}

/* ==========================================================================
 * Metric component computation
 * ========================================================================== */

double schwarzschild_g_tt(double r, double rs) {
    if (r <= 0.0 || rs <= 0.0) return 0.0;
    return -(1.0 - rs / r);
}

double schwarzschild_g_rr(double r, double rs) {
    if (r <= 0.0 || rs <= 0.0) return 0.0;
    if (r <= rs) return INFINITY;
    return 1.0 / (1.0 - rs / r);
}

double schwarzschild_g_theta_theta(double r) {
    return r * r;
}

double schwarzschild_g_phi_phi(double r, double theta) {
    return r * r * sin(theta) * sin(theta);
}

double schwarzschild_metric_determinant_root(double r, double theta) {
    return r * r * sin(theta);
}

void schwarzschild_metric_covariant(const SchwarzschildPoint *p, double rs,
                                     SchwarzschildSymmetric4x4 *metric) {
    if (!p || !metric) return;
    double r = p->r;
    double st = sin(p->theta);
    metric->g00 = -(1.0 - rs / r);
    metric->g01 = 0.0; metric->g02 = 0.0; metric->g03 = 0.0;
    metric->g11 = (r > rs) ? 1.0 / (1.0 - rs / r) : INFINITY;
    metric->g12 = 0.0; metric->g13 = 0.0;
    metric->g22 = r * r;
    metric->g23 = 0.0;
    metric->g33 = r * r * st * st;
}

void schwarzschild_metric_contravariant(const SchwarzschildPoint *p, double rs,
                                         SchwarzschildSymmetric4x4 *inv_metric) {
    if (!p || !inv_metric) return;
    double r = p->r;
    double st = sin(p->theta);
    double f = 1.0 - rs / r;
    inv_metric->g00 = (f != 0.0) ? -1.0 / f : INFINITY;
    inv_metric->g01 = 0.0; inv_metric->g02 = 0.0; inv_metric->g03 = 0.0;
    inv_metric->g11 = f;
    inv_metric->g12 = 0.0; inv_metric->g13 = 0.0;
    inv_metric->g22 = 1.0 / (r * r);
    inv_metric->g23 = 0.0;
    inv_metric->g33 = (st != 0.0) ? 1.0 / (r * r * st * st) : INFINITY;
}

/* ==========================================================================
 * Coordinate transformations
 * ========================================================================== */

double schwarzschild_to_tortoise(double r, double rs) {
    if (r <= 0.0 || rs <= 0.0) return 0.0;
    if (r <= rs) return -INFINITY;
    return r + rs * log((r / rs) - 1.0);
}

double tortoise_to_schwarzschild(double r_star, double rs) {
    /*
     * Solve r + rs*ln(r/rs - 1) = r* for r using Newton-Raphson.
     * f(r) = r + rs*ln((r/rs)-1) - r*
     * f'(r) = 1 + rs/(r - rs) = r/(r - rs)
     * Newton step: r_new = r - f(r)/f'(r) = r - (r + rs*ln((r/rs)-1) - r*) * (r-rs)/r
     */
    if (rs <= 0.0) return r_star;
    double r = r_star;
    if (r <= rs) r = rs * 1.001;
    for (int iter = 0; iter < 100; iter++) {
        double logarg = r / rs - 1.0;
        if (logarg <= 0.0) { r = rs * (1.0 + 1e-10); continue; }
        double f = r + rs * log(logarg) - r_star;
        double df = r / (r - rs);
        double dr = -f / df;
        r += dr;
        if (fabs(dr) < 1e-12 * fabs(r)) break;
    }
    return r;
}

double schwarzschild_to_isotropic(double r_schw, double rs) {
    /*
     * Inverse of: r_schw = r_iso * (1 + m/(2*r_iso))^2
     * where m = rs/2
     * Solve for r_iso.
     *
     * r_iso = (1/4) * ( sqrt(r_schw) + sqrt(r_schw - rs) )^2
     * This is the analytic inverse.
     */
    if (r_schw <= rs || rs <= 0.0) return 0.0;
    double s = sqrt(r_schw - rs);
    double t = sqrt(r_schw);
    return 0.25 * (t + s) * (t + s);
}

double isotropic_to_schwarzschild(double r_iso, double rs) {
    double m = rs / 2.0;
    double factor = 1.0 + m / (2.0 * r_iso);
    return r_iso * factor * factor;
}

void schwarzschild_to_kruskal(double t, double r, double rs,
                               double *T, double *X) {
    /*
     * For r > rs (Region I):
     *   T = sqrt(r/rs - 1) * exp(r/(2*rs)) * sinh(t/(2*rs))
     *   X = sqrt(r/rs - 1) * exp(r/(2*rs)) * cosh(t/(2*rs))
     *
     * For r < rs (Region II, black hole interior):
     *   T = sqrt(1 - r/rs) * exp(r/(2*rs)) * cosh(t/(2*rs))
     *   X = sqrt(1 - r/rs) * exp(r/(2*rs)) * sinh(t/(2*rs))
     *
     * Using geometric units (c = 1), t is in meters.
     * T^2 - X^2 = (1 - r/rs) * exp(r/rs)
     */
    if (!T || !X || rs <= 0.0) return;
    double exp_factor = exp(r / (2.0 * rs));
    if (r > rs) {
        double root = sqrt(r / rs - 1.0);
        double arg = t / (2.0 * rs);
        *T = root * exp_factor * sinh(arg);
        *X = root * exp_factor * cosh(arg);
    } else if (r < rs && r > 0.0) {
        double root = sqrt(1.0 - r / rs);
        double arg = t / (2.0 * rs);
        *T = root * exp_factor * cosh(arg);
        *X = root * exp_factor * sinh(arg);
    } else {
        *T = *X = 0.0;
    }
}

double schwarzschild_to_eddington_finkelstein_v(double t, double r, double rs) {
    double r_star = schwarzschild_to_tortoise(r, rs);
    return t + r_star;
}

double schwarzschild_to_eddington_finkelstein_u(double t, double r, double rs) {
    double r_star = schwarzschild_to_tortoise(r, rs);
    return t - r_star;
}

double schwarzschild_to_painleve_gullstrand_time(double t, double r, double rs) {
    if (r <= 0.0 || rs <= 0.0) return t;
    double root_term = sqrt(rs / r);
    double ln_arg;
    if (root_term < 1.0) {
        ln_arg = (1.0 + root_term) / (1.0 - root_term);
        return t + 2.0 * sqrt(r * rs) - rs * log(ln_arg);
    }
    return t;
}

double schwarzschild_proper_distance(double r1, double r2, double rs) {
    if (r1 <= rs || r2 <= rs || rs <= 0.0) return 0.0;
    if (r1 > r2) { double tmp = r1; r1 = r2; r2 = tmp; }
    /*
     * Analytic integral: integral dr/sqrt(1-rs/r) =
     *   sqrt(r*(r-rs)) + rs*ln(sqrt(r) + sqrt(r-rs)) + const
     *
     * In geometric terms:
     *   s(r) = sqrt(r^2 - rs*r) + (rs/2)*ln( (r-rs/2 + sqrt(r^2-rs*r)) / (rs/2) )
     *
     * Simpler: integrate numerically with Simpson's rule for robustness.
     */
    int n = 200;
    double dr = (r2 - r1) / n;
    double sum = 0.0;
    for (int i = 0; i <= n; i++) {
        double r = r1 + i * dr;
        double w = (i == 0 || i == n) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        sum += w / sqrt(1.0 - rs / r);
    }
    return sum * dr / 3.0;
}

double schwarzschild_volume_element(const SchwarzschildPoint *p, double rs) {
    (void)rs;
    if (!p) return 0.0;
    return schwarzschild_metric_determinant_root(p->r, p->theta);
}

/* ==========================================================================
 * Tetrad and local frame
 * ========================================================================== */

void schwarzschild_tetrad(const SchwarzschildPoint *p, double rs,
                           double tetrad[4][4]) {
    if (!p) return;
    double r = p->r;
    double f = sqrt(1.0 - rs / r);
    double st = sin(p->theta);
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            tetrad[a][b] = 0.0;
    if (r > rs) {
        tetrad[0][0] = f;
        tetrad[1][1] = 1.0 / f;
    }
    tetrad[2][2] = r;
    tetrad[3][3] = r * st;
}

void schwarzschild_project_to_tetrad(const SchwarzschildVector *v_coord,
                                      const double tetrad[4][4],
                                      SchwarzschildVector *v_tetrad) {
    if (!v_coord || !v_tetrad) return;
    double vc[4] = {v_coord->v0, v_coord->v1, v_coord->v2, v_coord->v3};
    double vt[4] = {0, 0, 0, 0};
    for (int a = 0; a < 4; a++)
        for (int mu = 0; mu < 4; mu++)
            vt[a] += tetrad[a][mu] * vc[mu];
    v_tetrad->v0 = vt[0]; v_tetrad->v1 = vt[1];
    v_tetrad->v2 = vt[2]; v_tetrad->v3 = vt[3];
}

/* ==========================================================================
 * Black hole initialization
 * ========================================================================== */

void schwarzschild_black_hole_init(double mass, SchwarzschildBlackHole *bh) {
    if (!bh) return;
    bh->mass = mass;
    bh->rs = schwarzschild_radius(mass);
    bh->geometric_mass = geometric_mass(mass);
    bh->horizon_area = schwarzschild_horizon_area(mass);
    bh->surface_gravity = (mass > 0.0) ? SCHW_C * SCHW_C * SCHW_C * SCHW_C / (4.0 * SCHW_G_N * mass) : 0.0;
    /*
     * Hawking temperature: T_H = hbar*c^3/(8*pi*G*M*k_B)
     * hbar = 1.054571817e-34, k_B = 1.380649e-23
     */
    double hbar = 1.054571817e-34;
    double kB  = 1.380649e-23;
    bh->hawking_temp = (mass > 0.0) ? (hbar * SCHW_C * SCHW_C * SCHW_C) / (8.0 * M_PI * SCHW_G_N * mass * kB) : 0.0;
    /* S = k_B*A*c^3/(4*G*hbar) */
    bh->bekenstein_entropy = (mass > 0.0) ? (kB * bh->horizon_area * SCHW_C * SCHW_C * SCHW_C) / (4.0 * SCHW_G_N * hbar) : 0.0;
}

int schwarzschild_is_inside_horizon(double r, double rs) {
    return (r < rs && r > 0.0) ? 1 : 0;
}

double schwarzschild_lapse(double r, double rs) {
    if (r <= rs) return 0.0;
    return sqrt(1.0 - rs / r);
}

void schwarzschild_shift(const SchwarzschildPoint *p, double rs, double beta[3]) {
    (void)p; (void)rs;
    beta[0] = 0.0; beta[1] = 0.0; beta[2] = 0.0;
}