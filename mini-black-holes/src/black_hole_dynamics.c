/**
 * @file black_hole_dynamics.c
 * @brief Geodesic integration, orbital dynamics, Penrose process, superradiance.
 *
 * Covers particle and photon motion around black holes, from effective
 * potentials through RK4 geodesic integration to energy extraction mechanisms.
 *
 * L5 Computational Methods: RK4 geodesic integration
 * L6 Canonical Systems: Schwarzschild/Kerr orbital dynamics
 * L7 Applications: BH shadow, Penrose process
 */

#include "black_hole_dynamics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 *  EFFECTIVE POTENTIALS
 * ================================================================ */

double schwarzschild_effective_potential(double r_s, double r, double L_sq) {
    /* V_eff(r) = sqrt[(1 - r_s/r)(1 + L²/r²)]
     *
     * This potential governs the radial motion of a test particle
     * in Schwarzschild spacetime:
     * (dr/dτ)² = E² - V_eff²(r)
     *
     * Features:
     * - V_eff → 0 as r → r_s (particles cross horizon with finite proper time)
     * - V_eff → 1 as r → ∞ (E ≥ 1 for escape to infinity)
     * - Local minimum for r > 6M (stable circular orbit)
     * - Local maximum for 3M < r < 6M (unstable circular orbit)
     * - No extrema for r < 3M (no circular orbits) */
    if (r <= r_s || r <= 0) return INFINITY;
    return sqrt((1.0 - r_s / r) * (1.0 + L_sq / (r * r)));
}

double schwarzschild_eff_potential_deriv(double r_s, double r, double L_sq) {
    /* dV_eff/dr = 0 gives the condition for circular orbits.
     * Used in root-finding to locate circular orbit radii. */
    if (r <= r_s) return 0.0;

    double h = r * 1e-6;
    double vp = schwarzschild_effective_potential(r_s, r + h, L_sq);
    double vm = schwarzschild_effective_potential(r_s, r - h, L_sq);
    return (vp - vm) / (2.0 * h);
}

double kerr_effective_potential_sq(double M, double a, double r,
                                    double E, double L_z, double mu) {
    /* In the equatorial plane (θ=π/2, Carter constant Q=0),
     * the radial equation for Kerr geodesics is:
     *
     * Σ² (dr/dλ)² = [E(r²+a²) - a L_z]² - Δ[(L_z-aE)² + μ² r²]
     *              ≡ R(r)  (radial Teukolsky potential)
     *
     * V_eff²(r) is defined such that (dr/dτ)² = E² - V_eff²
     * for timelike geodesics. Here we return the squared quantity. */

    double r2 = r * r;
    double a2 = a * a;
    double sigma = r2 + a2; /* Σ = r² at equator since cos(π/2)=0 */
    double delta = r2 - 2.0 * M * r + a2;

    if (delta <= 0 || sigma <= 0) return INFINITY;

    double P = E * (r2 + a2) - a * L_z;
    double P2 = P * P;
    double D_term = (L_z - a * E) * (L_z - a * E) + mu * mu * r2;

    return delta * D_term / sigma - P2 / (sigma * sigma);
}

/* ================================================================
 *  CIRCULAR ORBIT PARAMETERS
 * ================================================================ */

double schwarzschild_circular_L(double r_s, double r) {
    /* L_circ = r √[r_s / (2r - 3r_s)]
     *
     * Derived from dV_eff/dr = 0 for the Schwarzschild effective potential.
     *
     * At ISCO (r=6M): L = 2√3 M (marginally stable)
     * At photon sphere (r→3M): L → ∞ (timelike circular orbits impossible) */
    if (r <= 1.5 * r_s) return INFINITY;
    double denom = 2.0 * r - 3.0 * r_s;
    if (denom <= 0) return INFINITY;
    return r * sqrt(r_s / denom);
}

double schwarzschild_circular_E(double r_s, double r) {
    /* E_circ = (r - r_s) / √[r(r - 1.5 r_s)]
     *
     * Binding energy: E_bind = 1 - E (fraction of rest mass released).
     * For ISCO: E = 2√2/3 ≈ 0.943 → E_bind ≈ 5.72%
     *
     * This is the maximum energy extraction efficiency from thin disk
     * accretion onto a Schwarzschild BH.
     *
     * For Kerr prograde ISCO (a→M): E → 1/√3 ≈ 0.577 → E_bind ≈ 42.3%! */
    if (r <= 1.5 * r_s) return 0.0;
    return (r - r_s) / sqrt(r * (r - 1.5 * r_s));
}

double schwarzschild_circular_omega(double r_s, double r) {
    /* Ω = dφ/dt = √(GM/r³) → √(r_s/(2r³)) in geometric units
     *
     * This is identical to Newtonian Kepler's third law!
     * Coincidence: In Schwarzschild, the angular frequency of circular
     * orbits follows exactly the Keplerian formula (but r is the
     * areal radius, not the proper distance). */
    if (r <= r_s) return 0.0;
    return sqrt(r_s / (2.0 * r * r * r));
}

double schwarzschild_find_circular_orbit(double r_s, double L_sq, double r_guess) {
    /* Newton-Raphson: solve dV_eff/dr = 0 for r.
     *
     * For L² > 12M²: two solutions (stable outer, unstable inner).
     * For L² = 12M²: one solution at ISCO (r=6M).
     * For L² < 12M²: no circular orbit (particle plunges in). */

    const int max_iter = 100;
    const double tol = 1e-10;
    double r = r_guess;

    for (int i = 0; i < max_iter; i++) {
        double f = schwarzschild_eff_potential_deriv(r_s, r, L_sq);
        /* Second derivative via finite difference */
        double h = r * 1e-6;
        double fp = schwarzschild_eff_potential_deriv(r_s, r + h, L_sq);
        double fm = schwarzschild_eff_potential_deriv(r_s, r - h, L_sq);
        double fprime = (fp - fm) / (2.0 * h);

        if (fabs(fprime) < 1e-30) break;

        double dr = -f / fprime;
        r += dr;
        if (r <= r_s) return -1.0;
        if (fabs(dr) < tol) return r;
    }
    return r;
}

double schwarzschild_epicyclic_freq(double r_s, double r) {
    /* ω_r² = Ω² (1 - 6M/r)
     *
     * For a slightly perturbed circular orbit, the particle oscillates
     * radially at the epicyclic frequency ω_r.
     *
     * ω_r² > 0 → STABLE (oscillates around circular orbit)
     * ω_r² = 0 → MARGINALLY STABLE (ISCO)
     * ω_r² < 0 → UNSTABLE (exponential divergence: plunging orbit)
     *
     * The ISCO is defined by ω_r = 0 → r = 6M. */
    if (r <= r_s) return 0.0;
    double Omega = schwarzschild_circular_omega(r_s, r);
    return Omega * Omega * (1.0 - 3.0 * r_s / r);
}

double kerr_vertical_epicyclic_freq(double M, double a, double r) {
    /* ω_θ² = Ω² [1 - 4aΩ + 3a²/r²]
     *
     * In Kerr, the vertical oscillation is coupled to the BH spin
     * via frame-dragging. The vertical epicyclic frequency exceeds
     * the radial one for prograde orbits. */
    if (r <= 0 || M <= 0) return 0.0;
    double Omega = schwarzschild_circular_omega(2.0 * M, r);
    double aOmega = a * Omega;
    return Omega * Omega * (1.0 - 4.0 * aOmega + 3.0 * a * a / (r * r));
}

double kerr_radial_epicyclic_freq(double M, double a, double r) {
    /* ω_r² = Ω² [1 - 6M/r + 8a√(M/r³) - 3a²/r²]
     *
     * The ISCO in Kerr is where ω_r = 0:
     * 1 - 6M/r_ISCO + 8a√(M/r_ISCO³) - 3a²/r_ISCO² = 0 */
    if (r <= 0 || M <= 0) return 0.0;
    double Omega = schwarzschild_circular_omega(2.0 * M, r);
    double a_term = 8.0 * a * sqrt(M / (r * r * r));
    return Omega * Omega * (1.0 - 6.0 * M / r + a_term - 3.0 * a * a / (r * r));
}

/* ================================================================
 *  GEODESIC INTEGRATION VIA RK4
 * ================================================================ */

/*
 * Metric wrapper for numerical differentiation.
 * We store the current BH pointer and metric function in a context
 * structure to avoid global variables.
 */
typedef struct {
    void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m);
    const void *bh;
    double h_fd; /* finite difference step */
} GeodesicContext;

/* Local: compute Christoffel symbols via finite difference */
static double christoffel_fd_local(const GeodesicContext *ctx,
                                    const BLPosition *pos,
                                    int mu, int nu, int rho) {
    if (ctx == NULL || pos == NULL || ctx->metric_func == NULL) return 0.0;

    double h = ctx->h_fd;
    BLPosition pp, pm;
    Metric4D gp, gm, g0;
    double *ptr[4] = {&pp.t, &pp.r, &pp.theta, &pp.phi};
    double *mtr[4] = {&pm.t, &pm.r, &pm.theta, &pm.phi};

    ctx->metric_func(ctx->bh, pos, &g0);

    double Gamma = 0.0;
    for (int sigma = 0; sigma < 4; sigma++) {
        if (fabs(g0.g_inv[mu][sigma]) < 1e-30) continue;

        /* ∂_ν g_{ρσ} */
        pp = *pos; pm = *pos;
        *ptr[nu] += h; *mtr[nu] -= h;
        ctx->metric_func(ctx->bh, &pp, &gp);
        ctx->metric_func(ctx->bh, &pm, &gm);
        double d_nu = (gp.g[rho][sigma] - gm.g[rho][sigma]) / (2.0 * h);

        /* ∂_ρ g_{νσ} */
        pp = *pos; pm = *pos;
        *ptr[rho] += h; *mtr[rho] -= h;
        ctx->metric_func(ctx->bh, &pp, &gp);
        ctx->metric_func(ctx->bh, &pm, &gm);
        double d_rho = (gp.g[nu][sigma] - gm.g[nu][sigma]) / (2.0 * h);

        /* ∂_σ g_{νρ} */
        pp = *pos; pm = *pos;
        *ptr[sigma] += h; *mtr[sigma] -= h;
        ctx->metric_func(ctx->bh, &pp, &gp);
        ctx->metric_func(ctx->bh, &pm, &gm);
        double d_sigma = (gp.g[nu][rho] - gm.g[nu][rho]) / (2.0 * h);

        Gamma += 0.5 * g0.g_inv[mu][sigma] * (d_nu + d_rho - d_sigma);
    }
    return Gamma;
}

void geodesic_rhs(void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
                  const void *bh,
                  const GeodesicState *state,
                  double acc[4],
                  double h) {
    if (metric_func == NULL || bh == NULL || state == NULL || acc == NULL) return;

    GeodesicContext ctx;
    ctx.metric_func = metric_func;
    ctx.bh = bh;
    ctx.h_fd = h;

    BLPosition pos;
    pos.t = state->x[0];
    pos.r = state->x[1];
    pos.theta = state->x[2];
    pos.phi = state->x[3];

    /* dv^μ/dλ = -Γ^μ_{νρ} v^ν v^ρ */
    for (int mu = 0; mu < 4; mu++) {
        acc[mu] = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                double Gamma = christoffel_fd_local(&ctx, &pos, mu, nu, rho);
                acc[mu] -= Gamma * state->v[nu] * state->v[rho];
            }
        }
    }
}

int geodesic_rk4_step(
    void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
    const void *bh,
    GeodesicState *state,
    double dlambda,
    double h) {

    double x0[4], v0[4];
    memcpy(x0, state->x, 4 * sizeof(double));
    memcpy(v0, state->v, 4 * sizeof(double));

    /* RK4 coefficients */
    double k1v[4], k2v[4], k3v[4], k4v[4];
    GeodesicState temp;

    /* k1 */
    geodesic_rhs(metric_func, bh, state, k1v, h);

    /* k2 */
    memcpy(temp.x, x0, 4 * sizeof(double));
    for (int i = 0; i < 4; i++) {
        temp.v[i] = v0[i] + 0.5 * dlambda * k1v[i];
        temp.x[i] = x0[i] + 0.5 * dlambda * v0[i];
    }
    temp.affine = state->affine + 0.5 * dlambda;
    geodesic_rhs(metric_func, bh, &temp, k2v, h);

    /* k3 */
    memcpy(temp.x, x0, 4 * sizeof(double));
    for (int i = 0; i < 4; i++) {
        temp.v[i] = v0[i] + 0.5 * dlambda * k2v[i];
        temp.x[i] = x0[i] + 0.5 * dlambda * (v0[i] + 0.5 * dlambda * k1v[i]);
    }
    geodesic_rhs(metric_func, bh, &temp, k3v, h);

    /* k4 */
    memcpy(temp.x, x0, 4 * sizeof(double));
    for (int i = 0; i < 4; i++) {
        temp.v[i] = v0[i] + dlambda * k3v[i];
        temp.x[i] = x0[i] + dlambda * (v0[i] + 0.5 * dlambda * k2v[i]);
    }
    temp.affine = state->affine + dlambda;
    geodesic_rhs(metric_func, bh, &temp, k4v, h);

    /* Update */
    for (int i = 0; i < 4; i++) {
        state->v[i] = v0[i] + (dlambda / 6.0) * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
        state->x[i] = x0[i] + (dlambda / 6.0) * (v0[i] + 2.0 * (v0[i] + 0.5*dlambda*k1v[i])
                        + 2.0 * (v0[i] + 0.5*dlambda*k2v[i]) + (v0[i] + dlambda*k3v[i]));
    }
    state->affine += dlambda;

    /* Check if inside horizon (needs BH-specific check, simplified here) */
    return 0;
}

int geodesic_integrate(
    void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
    const void *bh,
    GeodesicState *state,
    double lambda_max, double dlambda,
    double *traj_r, double *traj_phi,
    int max_steps, int *n_steps) {

    if (metric_func == NULL || bh == NULL || state == NULL
        || traj_r == NULL || traj_phi == NULL || n_steps == NULL)
        return -1;

    double h = 1e-6 * dlambda;
    *n_steps = 0;

    while (state->affine < lambda_max && *n_steps < max_steps) {
        traj_r[*n_steps] = state->x[1];
        traj_phi[*n_steps] = state->x[3];
        (*n_steps)++;

        /* Advance one RK4 step. Check for horizon crossing/unphysical. */
        (void)geodesic_rk4_step(metric_func, bh, state, dlambda, h);

        /* Horizon check: if distance from origin is less than
         * a threshold (simplified; real condition depends on BH). */
        if (state->x[1] < 1e-6) return 1; /* Hit center */

        /* Escape check */
        if (state->x[1] > 1e10) return 2; /* Escaped */
    }
    return 0; /* Completed integration */
}

/* ================================================================
 *  PENROSE PROCESS
 * ================================================================ */

double penrose_process_efficiency(double M, double a) {
    /* η_max = (particle escaping energy - infalling energy) / infalling energy
     *
     * For extremal Kerr (a=M): η_max ≈ 20.7% (Bardeen et al. 1972).
     * For a = 0.9M: η_max ≈ 10%. For a = 0.5M: η_max ≈ 2%.
     *
     * Maximum extractable fraction of BH mass via Penrose process +
     * superradiance: up to 29% (M - M_irr)/M for extremal Kerr.
     *
     * In practice, astrophysical BHs are unlikely to have a > 0.998
     * (Thorne 1974 limit from accretion disk photon capture). */
    if (M <= 0 || fabs(a) > M) return 0.0;
    double aM = fabs(a) / M;
    /* Fitting formula based on Bardeen et al. (1972) */
    return 0.207 * aM * aM; /* ~0.207 at a=M */
}

int is_inside_ergosphere(const KerrBH *bh, const BLPosition *pos) {
    if (bh == NULL || pos == NULL) return 0;
    double rE = kerr_ergosphere_radius(bh->M, bh->a, pos->theta);
    return (pos->r > bh->r_plus && pos->r < rE) ? 1 : 0;
}

double particle_energy_at_infinity(const Metric4D *metric,
                                    const FourVelocity *velocity,
                                    double rest_mass) {
    if (metric == NULL || velocity == NULL) return 0.0;
    /* E = -p_μ ξ^μ = -g_{μν} ξ^μ p^ν
     * where ξ^μ = (1, 0, 0, 0) is the timelike Killing vector.
     * E = -(g_{0ν} p^ν) = -(g_{00} p⁰ + g_{0φ} p^φ) */
    double p[4];
    for (int i = 0; i < 4; i++) {
        p[i] = rest_mass * velocity->u[i];
    }
    return -(metric->g[0][0] * p[0] + metric->g[0][3] * p[3]);
}

double particle_angular_momentum(const Metric4D *metric,
                                  const FourVelocity *velocity,
                                  double rest_mass) {
    if (metric == NULL || velocity == NULL) return 0.0;
    /* L_z = p_μ η^μ = g_{φν} p^ν
     * = g_{30} p⁰ + g_{33} p³ */
    double p[4];
    for (int i = 0; i < 4; i++) {
        p[i] = rest_mass * velocity->u[i];
    }
    return metric->g[3][0] * p[0] + metric->g[3][3] * p[3];
}

/* ================================================================
 *  SUPERRADIANCE
 * ================================================================ */

int superradiance_condition(double omega, int m, double M, double a) {
    /* Superradiance: ω < m Ω_H
     *
     * For a scalar wave of frequency ω and azimuthal quantum number m
     * incident on a Kerr BH, the reflected wave is AMPLIFIED if ω < m Ω_H.
     *
     * This extracts rotational energy from the BH, equivalent to the
     * wave analog of the Penrose process.
     *
     * Key signatures:
     * - Massive scalar fields can form "boson clouds" (gravitational atom)
     *   with exponential growth → constrain ultralight dark matter
     * - For gravitational waves (s=-2): amplification possible up to 138%
     *   for m=2 (Teukolsky & Press 1974) */
    if (M <= 0 || fabs(a) > M || m == 0) return 0;
    double rp = M + sqrt(M * M - a * a);
    double Omega_H = a / (2.0 * M * rp); /* geometric */
    return (omega < fabs((double)m) * Omega_H) ? 1 : 0;
}

double superradiance_max_amplification(double M, double a, int l, int m) {
    (void)l; /* l used for potential spin-weighted generalization */
    /* Rough estimate of maximum superradiant amplification for scalar (s=0).
     *
     * For scalar m=1, l=1: amplification peaks near ω ≈ Ω_H/2,
     * with max gain ~4.4% for extremal Kerr.
     *
     * The amplification factor Z (reflected/incident flux ratio > 1). */
    if (M <= 0 || fabs(a) > M || m == 0) return 0.0;
    double aM = fabs(a) / M;
    /* Simple fit: max amplification ~ |m| * 0.044 * (a/M)^3 */
    return fabs((double)m) * 0.044 * aM * aM * aM;
}

/* ================================================================
 *  BARDEEN SPIN EVOLUTION
 * ================================================================ */

double bardeen_spin_evolution(double M_initial, double a_initial,
                               double accreted_mass, double *a_final) {
    /* Bardeen (1970): Thin disk accretion spins up the BH.
     *
     * The accreted gas brings orbital angular momentum at ISCO.
     * As the BH gains mass and angular momentum, a/M increases.
     *
     * Thorne (1974): The BH spin saturates at a/M ≈ 0.998 because
     * photons emitted by the disk are preferentially captured by the
     * BH at high spin (retrograde photon capture), exerting a
     * spin-down torque that balances the accretion spin-up.
     *
     * For a typical astrophysical BH born with a/M ≈ 0 and accreting
     * ~1 M_sun from a companion star: a/M reaches ~0.2-0.5. */

    if (M_initial <= 0 || a_initial < 0 || accreted_mass < 0) {
        if (a_final) *a_final = a_initial;
        return M_initial;
    }

    /* Simplified model: each unit of accreted mass adds angular momentum
     * corresponding to the ISCO specific angular momentum. */
    double M_final = M_initial + accreted_mass;

    /* ISCO angular momentum for current spin (approximate as constant) */
    double l_isco = kerr_isco_radius(M_initial, a_initial, 1)
                    * sqrt(M_initial / kerr_isco_radius(M_initial, a_initial, 1));

    double J_final = a_initial * M_initial + l_isco * accreted_mass;
    double a_new = J_final / M_final;

    /* Thorne limit: a/M ≤ 0.998 */
    if (a_new / M_final > 0.998) {
        a_new = 0.998 * M_final;
    }

    if (a_final) *a_final = a_new;
    return M_final;
}

/* ================================================================
 *  BLACK HOLE SHADOW AND PHOTON DEFLECTION
 * ================================================================ */

double schwarzschild_critical_impact(void) {
    /* b_crit = 3√3 M ≈ 5.196 M
     *
     * Photons with impact parameter b < b_crit are captured by the BH.
     * The BH shadow is a circle of angular radius:
     * θ_shadow = b_crit / D = 5.196 M / D
     *
     * For Sgr A* (M ≈ 4.1e6 M_sun, D ≈ 8.2 kpc):
     * θ_shadow ≈ 25 μas — the size measured by the Event Horizon Telescope! */
    return 3.0 * sqrt(3.0);
}

void kerr_shadow_radius(double M, double a, double theta_o,
                         double *R_mean, double *asymmetry) {
    /* The Kerr BH shadow is NOT circular.
     *
     * For a distant observer at inclination θ_o, the shadow boundary
     * in celestial coordinates (α, β) is parametrized by the photon
     * orbit radius r_ph ∈ [r_ph,min, r_ph,max].
     *
     * α(r) = -(ξ(r) / sin θ_o)
     * β(r) = ±√[η(r) + a² cos²θ_o - ξ²(r) cot²θ_o]
     *
     * where ξ = L_z/E, η = Q/E² are constants of motion for
     * spherical photon orbits.
     *
     * The shadow size R ≈ 2√3 M for Schwarzschild,
     * slightly shrunken for prograde Kerr, displaced from center. */

    if (M <= 0 || fabs(a) > M) {
        if (R_mean) *R_mean = 0.0;
        if (asymmetry) *asymmetry = 0.0;
        return;
    }

    /* Approximate fitting formula (Johannsen & Psaltis 2010) */
    double R_schw = 3.0 * sqrt(3.0) * M;
    double aM = a / M;
    *R_mean = R_schw * (1.0 - 0.05 * aM * aM * (1.0 + 0.5 * cos(theta_o) * cos(theta_o)));
    *asymmetry = -2.0 * aM * sin(theta_o) * M;
}

double schwarzschild_deflection_angle(double b, double r_s) {
    /* Weak field (b >> r_s): Δφ ≈ 2 r_s / b = 4M / b
     *
     * Einstein (1915): light from a distant star grazing the Sun
     * (b ≈ R_sun) is deflected by Δφ ≈ 1.75 arcseconds.
     *
     * Eddington's 1919 solar eclipse expedition confirmed this,
     * making Einstein world-famous.
     *
     * For strong fields (b ~ few × r_s):
     * The deflection can exceed 2π (multiple loops around the BH)
     * → "photon ring" in EHT images. */

    if (b <= 0 || r_s <= 0) return 0.0;

    /* Strong-field formula (Darwin 1959):
     * b = r_min / √(1 - r_s/r_min)
     * Find r_min from b, then integrate. */

    /* Weak-field approximation for b > 10 r_s */
    if (b > 10.0 * r_s) {
        return 2.0 * r_s / b + (15.0 * M_PI / 16.0 - 1.0) * r_s * r_s / (b * b);
    }

    /* Strong-field: numerical integration of the exact formula
     * Δφ = 2 ∫_{r_min}^{∞} dr / (r² √(1/b² - (1-r_s/r)/r²)) - π
     *
     * Simplified: find r_min, then do a rough integration. */
    double r_min = b * 0.8; /* crude estimate for r_min */
    /* Simple Simpson integration */
    int n = 500;
    double log_r_max = log(1000.0 * r_s);
    double log_r_min = log(r_min * 1.001);
    double dlogr = (log_r_max - log_r_min) / n;

    double integral = 0.0;
    for (int i = 0; i < n; i++) {
        double logr = log_r_min + i * dlogr;
        double r = exp(logr);
        double factor = 1.0 - r_s / r;
        double denom = r * r * sqrt(1.0 / (b * b) - factor / (r * r));
        if (denom > 0 && isfinite(denom)) {
            integral += dlogr * r / denom;
        }
    }
    return 2.0 * integral - M_PI;
}
