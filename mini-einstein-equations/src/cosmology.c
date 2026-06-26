/**
 * @file cosmology.c
 * @brief FLRW cosmology: Friedmann equations, scale factor, observables.
 *
 * Reference: Wald (1984), Ch.5; Carroll (2004), Ch.8
 *            Dodelson & Schmidt "Modern Cosmology" (2020)
 *            Weinberg "Cosmology" (2008)
 *
 * The FLRW metric is the unique homogeneous, isotropic solution
 * to the Einstein equations (Weyl's postulate + cosmological principle).
 */

#include "cosmology.h"
#include <math.h>
#include <string.h>

/* Conversion: 1 Mpc = 3.085677581e22 m */
#define MPC_TO_M  3.085677581e22
/* Conversion: km to m */
#define KM_TO_M   1000.0
/* Seconds per Gyr */
#define SEC_PER_GYR (1e9 * 365.25 * 24.0 * 3600.0)

/* ========================================================================
 * L1: Cosmological parameters
 * ========================================================================*/

void cosmology_params_init_lcdm(CosmologicalParams *p,
                                 double H0, double Omega_m, double Omega_L)
{
    p->H0 = H0;
    p->Omega_m = Omega_m;
    p->Omega_r = PLANCK_OMEGA_R;
    p->Omega_L = Omega_L;
    p->Omega_k = 1.0 - Omega_m - PLANCK_OMEGA_R - Omega_L;

    /* H0 in SI: (km/s/Mpc) → s^{-1} */
    p->H0_si = H0 * KM_TO_M / MPC_TO_M;
    /* H0/c in Mpc^{-1} */
    p->H0_over_c = H0 / (C_LIGHT * 1e-3); /* c in km/s ≈ 299792.458 km/s */
}

void cosmology_params_init_planck(CosmologicalParams *p)
{
    cosmology_params_init_lcdm(p,
        PLANCK_H0, PLANCK_OMEGA_M, PLANCK_OMEGA_L);
}

/* ========================================================================
 * L4: Friedmann equations
 * ========================================================================*/

/**
 * First Friedmann equation: Hubble parameter as function of scale factor.
 *
 *   H^2(a) = H_0^2 [ Ω_r a^{-4} + Ω_m a^{-3} + Ω_k a^{-2} + Ω_Λ ]
 *
 * The terms represent radiation, matter, curvature, and dark energy
 * respectively.
 *
 * Knowledge point: The Friedmann equation is the energy constraint of
 * GR applied to the FLRW metric. It determines the expansion history
 * of the universe.
 *
 * At early times (a→0): radiation dominates (∝ a^{-4})
 * Matter era (a~10^{-4}-0.5): matter dominates (∝ a^{-3})
 * Late times (a~1): dark energy dominates (constant)
 */
double friedmann_hubble(const CosmologicalParams *p, double a)
{
    double a_inv  = 1.0 / a;
    double a_inv2 = a_inv * a_inv;
    double a_inv3 = a_inv2 * a_inv;
    double a_inv4 = a_inv2 * a_inv2;

    double H2_over_H02 = p->Omega_r * a_inv4
                       + p->Omega_m * a_inv3
                       + p->Omega_k * a_inv2
                       + p->Omega_L;

    if (H2_over_H02 < 0.0) return 0.0; /* curvature dominated, collapsing */
    return p->H0 * sqrt(H2_over_H02);
}

/**
 * Second Friedmann equation (acceleration):
 *
 *   ä/a = -H_0^2 [ Ω_r a^{-4} + (1/2) Ω_m a^{-3} - Ω_Λ ]
 *
 * Knowledge point: The acceleration equation shows that radiation and
 * matter decelerate expansion, while dark energy (Ω_Λ > 0) accelerates it.
 *
 * The transition from deceleration to acceleration occurs when:
 *   Ω_r a^{-4} + (1/2) Ω_m a^{-3} = Ω_Λ
 *
 * For Planck ΛCDM, this occurs at a ≈ 0.75 (z ≈ 0.33).
 */
double friedmann_acceleration(const CosmologicalParams *p, double a)
{
    double a_inv  = 1.0 / a;
    double a_inv3 = a_inv * a_inv * a_inv;
    double a_inv4 = a_inv3 * a_inv;

    double term = p->Omega_r * a_inv4
                + 0.5 * p->Omega_m * a_inv3
                - p->Omega_L;

    return -p->H0 * p->H0 * term;
}

/**
 * Deceleration parameter: q(a) = -a ä / ȧ²
 *
 * q > 0: decelerating
 * q < 0: accelerating
 *
 * Present value q_0 ≈ -0.55 (Planck 2018) → universe is accelerating.
 */
double deceleration_parameter(const CosmologicalParams *p, double a)
{
    double H = friedmann_hubble(p, a);
    if (H < 1e-30) return 0.0;

    double a_ddot = friedmann_acceleration(p, a) * a;
    return -a * a_ddot / (H * H * a * a); /* = -(ä/a) / H^2 */
}

/**
 * Effective equation of state:
 *   w_eff = p_tot / ρ_tot = (Ω_r/3 - Ω_Λ) / (Ω_r + Ω_m + Ω_Λ)  at a=1
 *
 * For ΛCDM at present: w_eff ≈ -0.7 (dominated by Λ)
 * Radiation era (a→0): w_eff → 1/3
 * Matter era: w_eff → 0
 * Λ-dominated future: w_eff → -1
 */
double effective_eos_parameter(const CosmologicalParams *p, double a)
{
    double a_inv  = 1.0 / a;
    double a_inv2 = a_inv * a_inv;
    double a_inv3 = a_inv2 * a_inv;
    double a_inv4 = a_inv2 * a_inv2;

    double rho_r = p->Omega_r * a_inv4;
    double rho_m = p->Omega_m * a_inv3;
    double rho_L = p->Omega_L;

    double rho_tot = rho_r + rho_m + rho_L;
    if (rho_tot < 1e-30) return -1.0;

    /* p_tot = w_r*rho_r + w_m*rho_m + w_L*rho_L
     * w_r=1/3, w_m=0, w_L=-1 */
    double p_tot = (1.0/3.0) * rho_r + 0.0 * rho_m - 1.0 * rho_L;

    return p_tot / rho_tot;
}

/* ========================================================================
 * L5: Scale factor evolution
 * ========================================================================*/

/**
 * RHS: da/dt = a H(a)
 *
 * This ODE is the basis for cosmological integration.
 */
double scale_factor_ode(double a, const CosmologicalParams *p)
{
    return a * friedmann_hubble(p, a);
}

/**
 * Integrate scale factor using RK4.
 *
 * Knowledge point: The Friedmann equation is a first-order ODE in a(t).
 * Numerical integration is needed because the analytic solution involves
 * elliptic integrals for the general Ω_i case.
 *
 * Special cases with analytic solutions:
 *   Ω_m=1, Ω_Λ=0, Ω_r=0: a(t) ∝ t^{2/3} (Einstein-de Sitter)
 *   Ω_Λ=1, others=0: a(t) ∝ exp(H_0 t) (de Sitter)
 *   Ω_r=1, others=0: a(t) ∝ t^{1/2}
 */
void scale_factor_rk4_integrate(const CosmologicalParams *p,
                                 double a0, double t_max, int n_steps,
                                 double *t_out, double *a_out)
{
    double dt = t_max / n_steps;
    double a = a0;
    double t = 0.0;

    t_out[0] = t;
    a_out[0] = a;

    for (int i = 0; i < n_steps; i++) {
        /* RK4 for da/dt = a H(a) */
        double k1 = dt * scale_factor_ode(a, p);
        double k2 = dt * scale_factor_ode(a + 0.5*k1, p);
        double k3 = dt * scale_factor_ode(a + 0.5*k2, p);
        double k4 = dt * scale_factor_ode(a + k3, p);

        a += (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        t += dt;

        if (a < 0.0) a = 0.0; /* collapse */

        t_out[i + 1] = t;
        a_out[i + 1] = a;
    }
}

/**
 * Find cosmic time when a = a_target by integrating from a0.
 *
 * Uses the integral form:
 *   t = ∫_{a0}^{a_target} da / (a H(a))
 *
 * Knowledge point: The cosmic time at a given scale factor gives the
 * age of the universe at that epoch. The present age (a=1) is ~13.8 Gyr.
 */
double cosmic_time_at_a(const CosmologicalParams *p, double a_target, double a0)
{
    int n = 10000;
    double da = (a_target - a0) / n;
    double t = 0.0;
    double a = a0;

    for (int i = 0; i < n; i++) {
        /* E(a) = H(a)/H0 = sqrt(Ω_r a^{-4} + Ω_m a^{-3} + Ω_k a^{-2} + Ω_Λ) */
        double a_inv = 1.0 / a;
        double a_inv2 = a_inv * a_inv;
        double a_inv3 = a_inv2 * a_inv;
        double a_inv4 = a_inv2 * a_inv2;
        double E_sq = p->Omega_r * a_inv4 + p->Omega_m * a_inv3
                    + p->Omega_k * a_inv2 + p->Omega_L;
        if (E_sq < 1e-30) break;
        /* dt = da / (a * H(a)) = da / (a * H0 * sqrt(E_sq)) = (1/H0) * da/(a*E) */
        t += fabs(da) / (a * sqrt(E_sq));
        a += da;
    }

    return t; /* dimensionless, in units of Hubble time 1/H0 */
}

/* ========================================================================
 * L6: FLRW metric
 * ========================================================================*/

void flrw_metric(Metric *m, double a, double r, double k)
{
    /* Initialize all to zero */
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            m->g[mu][nu] = 0.0;
        }
    }

    double a2 = a * a;
    double denom = 1.0 - k * r * r;

    /* g_{00} = -1 */
    m->g[0][0] = -1.0;

    /* g_{11} = a^2 / (1 - k r^2) */
    m->g[1][1] = a2 / denom;

    /* g_{22} = a^2 r^2 */
    m->g[2][2] = a2 * r * r;

    /* g_{33} = a^2 r^2 sin^2(theta) — theta is the 3rd coordinate */
    /* Note: we set g_{33} assuming theta = π/2 for equatorial plane.
     * In general: g_{33} = a^2 r^2 sin^2(θ). */
    m->g[3][3] = a2 * r * r; /* equatorial approximation */

    metric_compute_inverse(m);
}

/* ========================================================================
 * L7: Cosmological observables
 * ========================================================================*/

/**
 * Age of universe: t_0 = ∫_0^1 da / (a H(a))
 *
 * For Planck ΛCDM: t_0 ≈ 13.8 Gyr
 */
double age_of_universe(const CosmologicalParams *p)
{
    /* Integrate from a_min to 1 */
    double a_min = 1e-8;
    double t_hubble = cosmic_time_at_a(p, 1.0, a_min);

    /* t_hubble is dimensionless (in Hubble time units 1/H0).
     * Convert to seconds: age_sec = t_hubble / H0_si
     * Convert to Gyr: age_Gyr = age_sec / SEC_PER_GYR */
    double age_sec = t_hubble / p->H0_si;
    return age_sec / SEC_PER_GYR;
}

double lookback_time(const CosmologicalParams *p, double a)
{
    double a_min = 1e-8;
    double t_to_a = cosmic_time_at_a(p, a, a_min);
    double t_to_now = cosmic_time_at_a(p, 1.0, a_min);
    double dt_hubble = t_to_now - t_to_a;
    double dt_sec = dt_hubble / p->H0_si;
    return dt_sec / SEC_PER_GYR;
}

/**
 * Comoving distance:
 *   χ(a) = ∫_a^1 da' / (a'^2 H(a'))
 *
 * In units of Hubble distance c/H_0.
 */
double comoving_distance(const CosmologicalParams *p, double a)
{
    int n = 10000;
    double da = (1.0 - a) / n;
    double chi = 0.0;
    double a_val = a;

    for (int i = 0; i < n; i++) {
        /* E(a) = H(a)/H0 */
        double a_inv = 1.0 / a_val;
        double a_inv2 = a_inv * a_inv;
        double a_inv3 = a_inv2 * a_inv;
        double a_inv4 = a_inv2 * a_inv2;
        double E_sq = p->Omega_r * a_inv4 + p->Omega_m * a_inv3
                    + p->Omega_k * a_inv2 + p->Omega_L;
        if (E_sq < 1e-30) break;
        double E = sqrt(E_sq);
        /* χ = ∫ da / (a^2 H(a)/H0) = ∫ da / (a^2 E(a)) */
        chi += da / (a_val * a_val * E);
        a_val += da;
    }

    return chi; /* in units of c/H0 */
}

double luminosity_distance(const CosmologicalParams *p, double a)
{
    double z = 1.0 / a - 1.0;  /* redshift */
    double chi = comoving_distance(p, a);

    /* d_L = (1+z) * d_M where d_M depends on curvature */
    double Omega_k = p->Omega_k;
    double H0_dist = C_LIGHT / (p->H0 * 1000.0); /* Hubble distance in meters */
    double d_M;

    if (fabs(Omega_k) < 1e-8) {
        d_M = chi * H0_dist;
    } else if (Omega_k > 0) {
        double sqrt_k = sqrt(Omega_k);
        d_M = H0_dist * sinh(sqrt_k * chi) / sqrt_k;
    } else {
        double sqrt_mk = sqrt(-Omega_k);
        d_M = H0_dist * sin(sqrt_mk * chi) / sqrt_mk;
    }

    return (1.0 + z) * d_M;
}

double angular_diameter_distance(const CosmologicalParams *p, double a)
{
    double z = 1.0 / a - 1.0;
    double d_L = luminosity_distance(p, a);
    return d_L / ((1.0 + z) * (1.0 + z));
}

/**
 * Critical density:
 *   ρ_c = 3 H_0^2 / (8π G)
 *
 * For H_0 = 70 km/s/Mpc: ρ_c ≈ 9.2 × 10^{-27} kg/m^3
 * Equivalent to ~5.5 protons per cubic meter.
 */
double critical_density(const CosmologicalParams *p)
{
    double H0_si = p->H0 * KM_TO_M / MPC_TO_M;
    return 3.0 * H0_si * H0_si / (8.0 * M_PI * G_NEWTON);
}

/* ========================================================================
 * L8: De Sitter inflation
 * ========================================================================*/

/**
 * de Sitter scale factor: a(t) = a_0 exp(H t)
 *
 * The de Sitter universe is the maximally symmetric vacuum solution
 * with positive cosmological constant. It describes:
 *   - The asymptotic future of ΛCDM
 *   - Cosmic inflation in the slow-roll approximation
 *   - The event horizon physics of black holes (near-horizon geometry)
 *
 * Knowledge point: Exponential expansion solves the horizon and flatness
 * problems of standard Big Bang cosmology (Guth 1981, Linde 1982).
 */
double de_sitter_scale_factor(double a0, double H, double t)
{
    return a0 * exp(H * t);
}

/**
 * Number of e-folds:
 *   N = ln(a_final / a_initial)
 *     = H Δt  (for de Sitter)
 *
 * Standard inflation requires N ≥ 50-60 e-folds to solve the horizon
 * and flatness problems.
 */
double inflation_e_folds(double a_initial, double a_final, double H, double dt)
{
    /* Two equivalent formulas for de Sitter: N = ln(a_f/a_i) = H Δt */
    double N_from_ratio = log(a_final / a_initial);
    double N_from_time = H * dt;
    /* Return average of the two (they're equal for exact de Sitter) */
    (void)N_from_time;
    return N_from_ratio;
}
