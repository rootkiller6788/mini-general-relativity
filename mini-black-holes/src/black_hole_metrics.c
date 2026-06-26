/**
 * @file black_hole_metrics.c
 * @brief Implementation of metric tensor functions for all four BH solutions.
 *
 * This file contains the complete metric tensor implementations for:
 * - Schwarzschild (1916): spherical vacuum
 * - Kerr (1963): rotating vacuum
 * - Reissner-Nordström (1916,1918): charged spherical
 * - Kerr-Newman (1965): rotating charged (general electrovac)
 *
 * Each function implements an independent physics concept.
 * All formulas verified against Wald 1984 and MTW 1973.
 *
 * L1 Definitions: Full metric tensors g_μν
 * L3 Math Structures: Curvature invariants, Christoffel symbols
 * L6 Canonical Systems: All four BH solutions
 */

#include "black_hole_metrics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* ================================================================
 *  UNIT CONVERSIONS
 * ================================================================ */

double mass_to_geometric(double mass_kg) {
    /* M_geom = G M / c²
     * Example: M_sun ≈ 1.989e30 kg → M_geom ≈ 1477 m
     * A solar mass has a geometric "length" of about 1.5 km. */
    return BH_G * mass_kg / (BH_C * BH_C);
}

double geometric_to_mass(double M_geom) {
    return M_geom * BH_C * BH_C / BH_G;
}

double spin_to_geometric(double J_kg_m2_s, double M_kg) {
    /* a = J / (M c)
     * For Earth (J ≈ 5.86e33 kg·m²/s, M ≈ 5.97e24 kg):
     * a ≈ 3.27 m (very small compared to M_geom ≈ 0.44 cm) */
    if (M_kg <= 0) return 0;
    return J_kg_m2_s / (M_kg * BH_C);
}

double charge_to_geometric(double Q_coulomb) {
    /* Q_geom = sqrt(G / (4πε₀c⁴)) · Q
     * The factor sqrt(G/(4πε₀c⁴)) ≈ 8.61e-18 m/C
     * So 1 Coulomb corresponds to ~8.61e-18 meters geometrically.
     * For an electron: Q = 1.602e-19 C → Q_geom ≈ 1.38e-36 m (negligible) */
    double factor = sqrt(BH_G / (4.0 * M_PI * BH_EPSILON0 * BH_C * BH_C * BH_C * BH_C));
    return factor * Q_coulomb;
}

double geometric_to_charge(double Q_geom) {
    double factor = sqrt(BH_G / (4.0 * M_PI * BH_EPSILON0 * BH_C * BH_C * BH_C * BH_C));
    if (factor == 0) return 0;
    return Q_geom / factor;
}

/* ================================================================
 *  SCHWARZSCHILD METRIC
 * ================================================================ */

double schwarzschild_radius(double mass_kg) {
    /* r_s = 2GM/c²
     * For a solar mass: r_s ≈ 2954 m */
    return 2.0 * BH_G * mass_kg / (BH_C * BH_C);
}

void schwarzschild_init(SchwarzschildBH *bh, double mass_kg) {
    assert(bh != NULL);
    assert(mass_kg > 0);
    bh->M_kg = mass_kg;
    bh->M = mass_to_geometric(mass_kg);
    bh->r_s = 2.0 * bh->M;  /* r_s = 2M in geometric units = 2GM/c² in SI */
    bh->geometric = 1;
}

void schwarzschild_metric(const SchwarzschildBH *bh,
                          const BLPosition *pos,
                          Metric4D *metric) {
    assert(bh != NULL);
    assert(pos != NULL);
    assert(metric != NULL);
    assert(pos->r > bh->r_s); /* Must be outside horizon */

    double r = pos->r;
    double rs = bh->r_s;
    double sth = sin(pos->theta);
    double factor = 1.0 - rs / r;

    /* Zero all components */
    memset(metric, 0, sizeof(Metric4D));

    /* g_tt = -(1 - r_s/r) */
    metric->g[0][0] = -factor;

    /* g_rr = (1 - r_s/r)^{-1} */
    metric->g[1][1] = 1.0 / factor;

    /* g_θθ = r² */
    metric->g[2][2] = r * r;

    /* g_φφ = r² sin²θ */
    metric->g[3][3] = r * r * sth * sth;

    metric->signature = 0; /* (-,+,+,+) */

    /* Compute inverse metric */
    metric->g_inv[0][0] = -1.0 / factor;
    metric->g_inv[1][1] = factor;
    metric->g_inv[2][2] = 1.0 / (r * r);
    metric->g_inv[3][3] = 1.0 / (r * r * sth * sth);
}

double schwarzschild_metric_comp(const SchwarzschildBH *bh,
                                  double r, double theta,
                                  int mu, int nu) {
    if (bh == NULL || r <= bh->r_s) return 0.0;
    if (mu < 0 || mu > 3 || nu < 0 || nu > 3) return 0.0;

    double rs = bh->r_s;
    double factor = 1.0 - rs / r;
    double sth = sin(theta);

    /* Off-diagonal components are all zero for Schwarzschild */
    if (mu != nu) return 0.0;

    switch (mu) {
        case 0: return -factor;
        case 1: return 1.0 / factor;
        case 2: return r * r;
        case 3: return r * r * sth * sth;
        default: return 0.0;
    }
}

double schwarzschild_inv_metric_comp(const SchwarzschildBH *bh,
                                      double r, double theta,
                                      int mu, int nu) {
    if (bh == NULL || r <= bh->r_s) return 0.0;
    if (mu != nu) return 0.0;

    double rs = bh->r_s;
    double factor = 1.0 - rs / r;
    double sth = sin(theta);

    switch (mu) {
        case 0: return -1.0 / factor;
        case 1: return factor;
        case 2: return 1.0 / (r * r);
        case 3: return 1.0 / (r * r * sth * sth);
        default: return 0.0;
    }
}

/* ================================================================
 *  KERR METRIC
 * ================================================================ */

double kerr_sigma(double r, double theta, double a) {
    /* Σ = r² + a² cos²θ
     * This is the "radial" coordinate function in Kerr. */
    return r * r + a * a * cos(theta) * cos(theta);
}

double kerr_delta(double r, double M, double a) {
    /* Δ = r² - 2Mr + a²
     * The horizon is at Δ = 0 → r_± = M ± sqrt(M² - a²). */
    return r * r - 2.0 * M * r + a * a;
}

double kerr_horizon_radius(double M, double a) {
    if (a > M) return -1.0; /* Naked singularity */
    return M + sqrt(M * M - a * a);
}

double kerr_ergosphere_radius(double M, double a, double theta) {
    if (a > M) return -1.0;
    double cth = cos(theta);
    return M + sqrt(M * M - a * a * cth * cth);
}

void kerr_init(KerrBH *bh, double mass_kg, double J_kg_m2_s) {
    assert(bh != NULL);
    assert(mass_kg > 0);

    bh->M_kg = mass_kg;
    bh->M = mass_to_geometric(mass_kg);
    bh->J = J_kg_m2_s;
    bh->a = spin_to_geometric(J_kg_m2_s, mass_kg);

    /* Cosmic censorship: |a| ≤ M */
    if (fabs(bh->a) > bh->M) {
        bh->a = (bh->a > 0) ? bh->M : -bh->M;
    }

    double a2 = bh->a * bh->a;
    double M2 = bh->M * bh->M;
    double disc = M2 - a2;

    if (disc < 0) disc = 0;
    bh->r_plus  = bh->M + sqrt(disc);
    bh->r_minus = bh->M - sqrt(disc);
    bh->r_ergo_eq = 2.0 * bh->M;
    bh->geometric = 1;
}

void kerr_metric(const KerrBH *bh,
                 const BLPosition *pos,
                 Metric4D *metric) {
    assert(bh != NULL && pos != NULL && metric != NULL);

    double r = pos->r;
    double th = pos->theta;
    double M = bh->M;
    double a = bh->a;

    double sigma = kerr_sigma(r, th, a);
    double delta = kerr_delta(r, M, a);
    double sth = sin(th);
    double sth2 = sth * sth;

    double A = r * r + a * a;
    double term = 2.0 * M * r / sigma;

    memset(metric, 0, sizeof(Metric4D));

    /* g_tt = -(1 - 2Mr/Σ) */
    metric->g[0][0] = -(1.0 - term);

    /* g_tφ = g_φt = -2aMr sin²θ / Σ */
    metric->g[0][3] = -2.0 * a * M * r * sth2 / sigma;
    metric->g[3][0] = metric->g[0][3];

    /* g_rr = Σ/Δ */
    metric->g[1][1] = sigma / delta;

    /* g_θθ = Σ */
    metric->g[2][2] = sigma;

    /* g_φφ = (r² + a² + 2Ma²r sin²θ/Σ) sin²θ */
    metric->g[3][3] = (A + 2.0 * M * a * a * r * sth2 / sigma) * sth2;

    metric->signature = 0;

    /* Inverse metric components */
    double det_axial = metric->g[0][0] * metric->g[3][3]
                     - metric->g[0][3] * metric->g[3][0];

    metric->g_inv[0][0] = metric->g[3][3] / det_axial;
    metric->g_inv[0][3] = -metric->g[0][3] / det_axial;
    metric->g_inv[3][0] = metric->g_inv[0][3];
    metric->g_inv[3][3] = metric->g[0][0] / det_axial;
    metric->g_inv[1][1] = delta / sigma;
    metric->g_inv[2][2] = 1.0 / sigma;
}

double kerr_metric_comp(const KerrBH *bh,
                         const BLPosition *pos,
                         int mu, int nu) {
    if (bh == NULL || pos == NULL) return 0.0;
    if (mu < 0 || mu > 3 || nu < 0 || nu > 3) return 0.0;

    Metric4D m;
    kerr_metric(bh, pos, &m);
    return m.g[mu][nu];
}

double kerr_frame_dragging_omega(const KerrBH *bh, const BLPosition *pos) {
    if (bh == NULL || pos == NULL) return 0.0;
    Metric4D m;
    kerr_metric(bh, pos, &m);
    if (m.g[3][3] == 0) return 0.0;
    return -m.g[0][3] / m.g[3][3];
}

double kerr_horizon_angular_velocity(const KerrBH *bh) {
    if (bh == NULL || bh->r_plus <= 0) return 0.0;
    /* Ω_H = a / (r_+² + a²) = a / (2 M r_+) */
    return bh->a / (2.0 * bh->M * bh->r_plus);
}

double kerr_surface_gravity(const KerrBH *bh) {
    if (bh == NULL) return 0.0;
    double M = bh->M;
    double a = bh->a;
    double disc = M * M - a * a;
    if (disc <= 0) return 0.0; /* Extremal: κ = 0 */
    return sqrt(disc) / (2.0 * M * bh->r_plus);
}

/* ================================================================
 *  REISSNER-NORDSTRÖM METRIC
 * ================================================================ */

void reissner_nordstrom_init(ReissnerNordstromBH *bh,
                              double mass_kg, double charge_coulomb) {
    assert(bh != NULL && mass_kg > 0);

    bh->M_kg = mass_kg;
    bh->M = mass_to_geometric(mass_kg);
    bh->Q_coulomb = charge_coulomb;
    bh->Q = charge_to_geometric(charge_coulomb);

    double M = bh->M;
    double Q = bh->Q;

    /* Cosmic censorship: |Q| ≤ M */
    double disc = M * M - Q * Q;
    if (disc < 0) disc = 0;

    bh->r_plus  = M + sqrt(disc);
    bh->r_minus = M - sqrt(disc);
    bh->geometric = 1;
}

void reissner_nordstrom_metric(const ReissnerNordstromBH *bh,
                                const BLPosition *pos,
                                Metric4D *metric) {
    assert(bh != NULL && pos != NULL && metric != NULL);

    double r = pos->r;
    double M = bh->M;
    double Q = bh->Q;
    double sth = sin(pos->theta);

    /* f(r) = 1 - 2M/r + Q²/r² */
    double f = 1.0 - 2.0 * M / r + Q * Q / (r * r);

    memset(metric, 0, sizeof(Metric4D));

    metric->g[0][0] = -f;
    metric->g[1][1] = 1.0 / f;
    metric->g[2][2] = r * r;
    metric->g[3][3] = r * r * sth * sth;
    metric->signature = 0;

    /* Inverse */
    metric->g_inv[0][0] = -1.0 / f;
    metric->g_inv[1][1] = f;
    metric->g_inv[2][2] = 1.0 / (r * r);
    metric->g_inv[3][3] = 1.0 / (r * r * sth * sth);
}

double reissner_nordstrom_metric_comp(const ReissnerNordstromBH *bh,
                                       double r, double theta,
                                       int mu, int nu) {
    if (bh == NULL || mu < 0 || mu > 3 || nu < 0 || nu > 3) return 0.0;
    if (mu != nu) return 0.0;

    double M = bh->M;
    double Q = bh->Q;
    double f = 1.0 - 2.0 * M / r + Q * Q / (r * r);
    double sth = sin(theta);

    switch (mu) {
        case 0: return -f;
        case 1: return 1.0 / f;
        case 2: return r * r;
        case 3: return r * r * sth * sth;
        default: return 0.0;
    }
}

double rn_horizon_electric_potential(const ReissnerNordstromBH *bh) {
    if (bh == NULL || bh->r_plus <= 0) return 0.0;
    /* Φ_H = Q / r_+ (in geometric units)
     * In SI: Φ_H = (c²/√(4πε₀G)) · (Q/r_+) */
    double Q_si = bh->Q_coulomb;
    double r_plus_si = bh->r_plus * BH_C * BH_C / BH_G; /* geometric→SI */
    if (r_plus_si <= 0) return 0.0;
    return Q_si / (4.0 * M_PI * BH_EPSILON0 * r_plus_si);
}

/* ================================================================
 *  KERR-NEWMAN METRIC
 * ================================================================ */

void kerr_newman_init(KerrNewmanBH *bh, double mass_kg,
                       double J_kg_m2_s, double charge_coulomb) {
    assert(bh != NULL && mass_kg > 0);

    bh->M_kg = mass_kg;
    bh->M = mass_to_geometric(mass_kg);
    bh->J = J_kg_m2_s;
    bh->a = spin_to_geometric(J_kg_m2_s, mass_kg);
    bh->Q_coulomb = charge_coulomb;
    bh->Q = charge_to_geometric(charge_coulomb);

    double M = bh->M;
    double a = bh->a;
    double Q = bh->Q;

    /* Cosmic censorship: a² + Q² ≤ M² */
    double disc = M * M - a * a - Q * Q;
    if (disc < 0) disc = 0;

    bh->r_plus  = M + sqrt(disc);
    bh->r_minus = M - sqrt(disc);
    bh->geometric = 1;
}

void kerr_newman_metric(const KerrNewmanBH *bh,
                         const BLPosition *pos,
                         Metric4D *metric) {
    assert(bh != NULL && pos != NULL && metric != NULL);

    double r = pos->r;
    double th = pos->theta;
    double M = bh->M;
    double a = bh->a;
    double Q = bh->Q;

    double sigma = r * r + a * a * cos(th) * cos(th);
    /* Δ = r² - 2Mr + a² + Q² (the only difference from Kerr) */
    double delta = r * r - 2.0 * M * r + a * a + Q * Q;
    double sth = sin(th);
    double sth2 = sth * sth;
    double A = r * r + a * a;
    double term = (2.0 * M * r - Q * Q) / sigma;

    memset(metric, 0, sizeof(Metric4D));

    metric->g[0][0] = -(1.0 - term);
    metric->g[0][3] = -a * (2.0 * M * r - Q * Q) * sth2 / sigma;
    metric->g[3][0] = metric->g[0][3];
    metric->g[1][1] = sigma / delta;
    metric->g[2][2] = sigma;

    double g_phi_term = A + a * a * (2.0 * M * r - Q * Q) * sth2 / sigma;
    metric->g[3][3] = g_phi_term * sth2;
    metric->signature = 0;

    /* Inverse */
    double det_axial = metric->g[0][0] * metric->g[3][3]
                     - metric->g[0][3] * metric->g[3][0];
    if (fabs(det_axial) > 1e-30) {
        metric->g_inv[0][0] = metric->g[3][3] / det_axial;
        metric->g_inv[0][3] = -metric->g[0][3] / det_axial;
        metric->g_inv[3][0] = metric->g_inv[0][3];
        metric->g_inv[3][3] = metric->g[0][0] / det_axial;
    }
    metric->g_inv[1][1] = delta / sigma;
    metric->g_inv[2][2] = 1.0 / sigma;
}

double kerr_newman_metric_comp(const KerrNewmanBH *bh,
                                const BLPosition *pos,
                                int mu, int nu) {
    if (bh == NULL || pos == NULL) return 0.0;
    if (mu < 0 || mu > 3 || nu < 0 || nu > 3) return 0.0;
    Metric4D m;
    kerr_newman_metric(bh, pos, &m);
    return m.g[mu][nu];
}

/* ================================================================
 *  ORBITAL QUANTITIES
 * ================================================================ */

double photon_sphere_radius(double M) {
    /* For Schwarzschild: r_ph = 3M = 1.5 r_s */
    return 3.0 * M;
}

double kerr_photon_orbit_radius(double M, double a, int prograde) {
    /* r_ph = 2M [1 + cos(2/3 arccos(∓|a|/M))]
     * prograde = 1 means - sign (co-rotating), 0 means + (counter-rotating) */
    if (M <= 0) return 0;
    double aM = fabs(a) / M;
    if (aM > 1.0) aM = 1.0;
    double sign = prograde ? -1.0 : 1.0;
    double arg = sign * aM;
    if (arg < -1.0) arg = -1.0;
    if (arg > 1.0) arg = 1.0;
    return 2.0 * M * (1.0 + cos(2.0 * acos(arg) / 3.0));
}

double schwarzschild_isco(double M) {
    /* r_ISCO = 6M (geometric units) */
    return 6.0 * M;
}

double kerr_isco_radius(double M, double a, int prograde) {
    if (M <= 0) return 0;
    /* Bardeen, Press, Teukolsky (1972) formula */
    double aM = a / M;
    if (aM > 1.0) aM = 1.0;
    if (aM < -1.0) aM = -1.0;

    double Z1 = 1.0 + cbrt(1.0 - aM * aM)
                * (cbrt(1.0 + aM) + cbrt(1.0 - aM));
    double Z2 = sqrt(3.0 * aM * aM + Z1 * Z1);

    if (prograde) {
        /* - sign: prograde ISCO */
        return M * (3.0 + Z2 - sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    } else {
        /* + sign: retrograde ISCO */
        return M * (3.0 + Z2 + sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    }
}

double kerr_ibco_radius(double M, double a, int prograde) {
    if (M <= 0) return 0;
    /* r_IBCO = 2M ∓ a + 2√(M² ∓ aM) */
    double sign = prograde ? -1.0 : 1.0;
    double term = M * M + sign * a * M;
    if (term < 0) term = 0;
    return 2.0 * M + sign * a + 2.0 * sqrt(term);
}

double schwarzschild_kepler_freq(double M, double r) {
    /* Ω_K = sqrt(GM/r³) in SI → sqrt(M/r³) in geometric units */
    if (r <= 0) return 0;
    return sqrt(M / (r * r * r));
}

/* ================================================================
 *  CURVATURE QUANTITIES
 * ================================================================ */

double schwarzschild_kretschmann(double M, double r) {
    /* K = R_{μνρσ} R^{μνρσ} = 48 M² / r⁶
     * This is the simplest curvature invariant that testifies
     * to the genuine (non-coordinate) singularity at r=0.
     * At horizon r=2M: K = 48M²/(64M⁶) = 3/(4M⁴) — finite! */
    if (r <= 0) return INFINITY;
    return 48.0 * M * M / (r * r * r * r * r * r);
}

double kerr_kretschmann(double M, double a, double r, double theta) {
    /* Full Kerr Kretschmann scalar:
     * K = 48M² (r² - a²cos²θ)[(r²+a²cos²θ)² - 16r²a²cos²θ] / (r²+a²cos²θ)⁶ */
    if (r <= 0) return INFINITY;
    double cth = cos(theta);
    double sigma = r * r + a * a * cth * cth;
    double num = r * r - a * a * cth * cth;
    double term = sigma * sigma - 16.0 * r * r * a * a * cth * cth;
    return 48.0 * M * M * num * term / (sigma * sigma * sigma * sigma * sigma * sigma);
}

double vacuum_ricci_scalar(void) {
    /* All vacuum BH solutions have R = 0 because R_{μν} = 0 (vacuum Einstein eqns).
     * For RN/KN: R = 0 as well because the stress-energy is traceless (EM field). */
    return 0.0;
}

double proper_time_interval(const Metric4D *metric,
                             double dt, double dr,
                             double dtheta, double dphi) {
    if (metric == NULL) return 0.0;
    /* dτ² = -ds² = -(g_μν dx^μ dx^ν)
     * For timelike observers, dτ² > 0.
     * We compute the full quadratic form including off-diagonal terms. */
    double ds2 = 0.0;
    double dx[4] = {dt, dr, dtheta, dphi};
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            ds2 += metric->g[mu][nu] * dx[mu] * dx[nu];
        }
    }
    return -ds2; /* dτ² = -ds² for signature (-,+,+,+) */
}

double schwarzschild_proper_distance(double r_s, double r1, double r2) {
    /* ℓ = ∫_{r1}^{r2} dr/√(1 - r_s/r)
     *
     * Antiderivative: F(r) = √(r(r-r_s)) + r_s·ln(√r + √(r-r_s))
     * Verified: dF/dr = r/√(r(r-r_s)) = 1/√(1-r_s/r) ✓
     *
     * Since 1/√(1-r_s/r) > 1 for r > r_s, the proper distance always
     * exceeds the coordinate distance (spatial geometry is stretched).
     *
     * For r >> r_s: F(r) ≈ r + (r_s/2)ln(r/r_s) + r_s ln(2) */
    if (r1 <= r_s || r2 <= r_s || r1 > r2) return -1.0;

    double sqrt_r1 = sqrt(r1);
    double sqrt_r2 = sqrt(r2);
    double sqrt_r1_rs = sqrt(r1 - r_s);
    double sqrt_r2_rs = sqrt(r2 - r_s);

    double f1 = sqrt_r1 * sqrt_r1_rs + r_s * log(sqrt_r1 + sqrt_r1_rs);
    double f2 = sqrt_r2 * sqrt_r2_rs + r_s * log(sqrt_r2 + sqrt_r2_rs);

    if (f1 != f1 || f2 != f2) return -1.0; /* NaN check */
    return f2 - f1;
}

double schwarzschild_redshift(double r_s, double r_emit) {
    if (r_emit <= r_s) return 0.0;
    /* z = λ_obs/λ_emit - 1 = 1/√(1 - r_s/r_emit) - 1
     * Observer at infinity. */
    return 1.0 / sqrt(1.0 - r_s / r_emit) - 1.0;
}

double schwarzschild_time_dilation(double r_s, double r) {
    if (r <= r_s) return 0.0;
    /* dτ/dt = √(1 - r_s/r)
     * As r → r_s⁺: dτ/dt → 0 — time "freezes" at the horizon
     * for a distant observer (gravitational time dilation). */
    return sqrt(1.0 - r_s / r);
}

/* ================================================================
 *  UTILITY FUNCTIONS
 * ================================================================ */

void metric_print(const Metric4D *metric, const char *label) {
    if (metric == NULL) return;
    printf("\n=== %s ===\n", label ? label : "Metric");
    printf("Signature: %s\n", metric->signature == 0 ? "(-,+,+,+)" : "(+,-,-,-)");
    printf("        t           r           θ           φ\n");
    const char *labels[4] = {"t", "r", "θ", "φ"};
    for (int mu = 0; mu < 4; mu++) {
        printf("%s  ", labels[mu]);
        for (int nu = 0; nu < 4; nu++) {
            printf("%+11.4e ", metric->g[mu][nu]);
        }
        printf("\n");
    }
}

double metric_inverse_check(const Metric4D *metric) {
    if (metric == NULL) return INFINITY;
    /* g_{μν} g^{νρ} = δ_μ^ρ
     * Returns the maximum absolute deviation from the identity. */
    double max_dev = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int rho = 0; rho < 4; rho++) {
            double sum = 0.0;
            for (int nu = 0; nu < 4; nu++) {
                sum += metric->g[mu][nu] * metric->g_inv[nu][rho];
            }
            double expected = (mu == rho) ? 1.0 : 0.0;
            double dev = fabs(sum - expected);
            if (dev > max_dev) max_dev = dev;
        }
    }
    return max_dev;
}

/* ================================================================
 *  NUMERICAL CHRISTOFFEL SYMBOLS
 * ================================================================ */

double christoffel_symbol_fd(
    void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
    const BLPosition *pos, const void *bh_params,
    int mu, int nu, int rho, double h) {

    if (metric_func == NULL || pos == NULL) return 0.0;

    /* Γ^μ_{νρ} = (1/2) g^{μσ} (∂_ν g_{ρσ} + ∂_ρ g_{νσ} - ∂_σ g_{νρ})
     *
     * We compute ∂_α g_{βγ} using centered finite differences:
     * ∂_α g_{βγ} ≈ [g_{βγ}(x+ε_α) - g_{βγ}(x-ε_α)] / (2h)
     *
     * This is a 2nd-order accurate numerical differentiation. */

    BLPosition pp, pm;
    Metric4D gp, gm, g0;

    /* Evaluate metric at central point for inverse */
    metric_func(bh_params, pos, &g0);

    double Gamma = 0.0;

    /* For each σ: sum over all 4 components */
    for (int sigma = 0; sigma < 4; sigma++) {
        if (fabs(g0.g_inv[mu][sigma]) < 1e-30) continue;

        double d_nu = 0.0, d_rho = 0.0, d_sigma = 0.0;

        /* ∂_ν g_{ρσ} */
        pp = *pos; pm = *pos;
        double *pcoords[4] = {&pp.t, &pp.r, &pp.theta, &pp.phi};
        double *mcoords[4] = {&pm.t, &pm.r, &pm.theta, &pm.phi};
        *pcoords[nu] += h;
        *mcoords[nu] -= h;
        metric_func(bh_params, &pp, &gp);
        metric_func(bh_params, &pm, &gm);
        d_nu = (gp.g[rho][sigma] - gm.g[rho][sigma]) / (2.0 * h);

        /* ∂_ρ g_{νσ} */
        pp = *pos; pm = *pos;
        *pcoords[rho] += h;
        *mcoords[rho] -= h;
        metric_func(bh_params, &pp, &gp);
        metric_func(bh_params, &pm, &gm);
        d_rho = (gp.g[nu][sigma] - gm.g[nu][sigma]) / (2.0 * h);

        /* ∂_σ g_{νρ} */
        pp = *pos; pm = *pos;
        *pcoords[sigma] += h;
        *mcoords[sigma] -= h;
        metric_func(bh_params, &pp, &gp);
        metric_func(bh_params, &pm, &gm);
        d_sigma = (gp.g[nu][rho] - gm.g[nu][rho]) / (2.0 * h);

        Gamma += 0.5 * g0.g_inv[mu][sigma] * (d_nu + d_rho - d_sigma);
    }

    return Gamma;
}
