/**
 * @file cosmology.h
 * @brief FLRW cosmology: Friedmann equations, scale factor evolution.
 *
 * Reference: Wald (1984), Ch.5; Carroll (2004), Ch.8
 *            Dodelson & Schmidt "Modern Cosmology" (2020)
 *            Weinberg "Cosmology" Ch.1-3
 *
 * Knowledge map:
 *   L4: Friedmann equations
 *   L6: FLRW metric, scale factor evolution
 *   L7: Cosmological parameters (H_0, Omega_m, Omega_Lambda, Omega_k)
 *   L7: Age of the universe, luminosity distance
 *   L8: Inflation basics (de Sitter expansion)
 */

#ifndef MINI_EINSTEIN_COSMOLOGY_H
#define MINI_EINSTEIN_COSMOLOGY_H

#include "einstein.h"

/* ---------------------------------------------------------------------------
 * L1: FLRW metric parameters
 * -------------------------------------------------------------------------*/

/** FLRW cosmological model parameters */
typedef struct {
    double H0;         /**< Hubble constant (km/s/Mpc) */
    double Omega_m;    /**< Matter density parameter */
    double Omega_r;    /**< Radiation density parameter */
    double Omega_L;    /**< Dark energy (cosmological constant) density parameter */
    double Omega_k;    /**< Curvature density parameter (= 1 - Omega_m - Omega_r - Omega_L) */
    double H0_si;      /**< Hubble constant in SI (s^{-1}) */
    double H0_over_c;  /**< H0 / c in Mpc^{-1} */
} CosmologicalParams;

/* Planck 2018 best-fit values */
#define PLANCK_H0        67.66    /**< km/s/Mpc */
#define PLANCK_OMEGA_M   0.3111
#define PLANCK_OMEGA_L   0.6889
#define PLANCK_OMEGA_R   9.182e-5

/**
 * Initialize cosmological parameters from standard Lambda-CDM.
 * Sets Omega_r = PLANCK_OMEGA_R, Omega_k = 1 - Omega_m - Omega_r - Omega_L.
 */
void cosmology_params_init_lcdm(CosmologicalParams *p,
                                 double H0, double Omega_m, double Omega_L);

/**
 * Initialize empty FLRW params with Planck 2018 defaults.
 */
void cosmology_params_init_planck(CosmologicalParams *p);

/* ---------------------------------------------------------------------------
 * L4: Friedmann equations
 * -------------------------------------------------------------------------*/

/**
 * First Friedmann equation (energy constraint):
 *   H^2(a) = H_0^2 * [Omega_r/a^4 + Omega_m/a^3 + Omega_k/a^2 + Omega_L]
 *
 * where H(a) = (da/dt)/a is the Hubble parameter at scale factor a.
 *
 * @param p   cosmological parameters
 * @param a   scale factor (> 0)
 * @return H(a) in km/s/Mpc
 */
double friedmann_hubble(const CosmologicalParams *p, double a);

/**
 * Second Friedmann equation (acceleration):
 *   d^2a/dt^2 / a = -H_0^2 * [Omega_r/a^4 + (1/2)*Omega_m/a^3 - Omega_L]
 *
 * @param p   cosmological parameters
 * @param a   scale factor
 * @return (d^2a/dt^2) / a in (km/s/Mpc)^2
 */
double friedmann_acceleration(const CosmologicalParams *p, double a);

/**
 * Deceleration parameter:
 *   q(a) = -a * (d^2a/dt^2) / (da/dt)^2
 *
 * For Lambda-CDM, q crosses through 0 at late times (accelerating now).
 */
double deceleration_parameter(const CosmologicalParams *p, double a);

/**
 * Equation of state parameter w(a):
 *   w_eff = p_eff / rho_eff
 * For mixture: w = (Omega_r/3 - Omega_L) / (Omega_r + Omega_m + Omega_L) at a=1.
 */
double effective_eos_parameter(const CosmologicalParams *p, double a);

/* ---------------------------------------------------------------------------
 * L5: Scale factor evolution (numerical integration)
 * -------------------------------------------------------------------------*/

/**
 * ODE right-hand side for scale factor:
 *   da/dt = a * H(a)
 *
 * State is single variable a(t).
 */
double scale_factor_ode(double a, const CosmologicalParams *p);

/**
 * Integrate scale factor from initial a0 to final a_final using
 * fourth-order Runge-Kutta.
 *
 * @param p         cosmological parameters
 * @param a0        initial scale factor
 * @param t_max     maximum cosmic time to integrate (in Hubble times 1/H0)
 * @param n_steps   number of integration steps
 * @param t_out     output time array (n_steps+1 values)
 * @param a_out     output scale factor array (n_steps+1 values)
 */
void scale_factor_rk4_integrate(const CosmologicalParams *p,
                                 double a0, double t_max, int n_steps,
                                 double *t_out, double *a_out);

/**
 * Find cosmic time when a equals given value.
 *
 * @param p       cosmological parameters
 * @param a_target target scale factor
 * @param a0      initial scale factor (for integration start)
 * @return cosmic time in Hubble times (1/H0)
 */
double cosmic_time_at_a(const CosmologicalParams *p, double a_target, double a0);

/* ---------------------------------------------------------------------------
 * L6: FLRW metric components
 * -------------------------------------------------------------------------*/

/**
 * Fill FLRW metric at comoving coordinates (t, r, theta, phi)
 * with scale factor a(t):
 *
 *   ds^2 = -dt^2 + a^2(t) * [dr^2/(1 - k*r^2) + r^2 dOmega^2]
 *
 * where k = -Omega_k * H_0^2 / c^2 (curvature parameter).
 *
 * @param m    output metric
 * @param a    scale factor
 * @param r    comoving radial coordinate
 * @param k    curvature parameter (>0 closed, <0 open, =0 flat)
 */
void flrw_metric(Metric *m, double a, double r, double k);

/* ---------------------------------------------------------------------------
 * L7: Cosmological observables
 * -------------------------------------------------------------------------*/

/**
 * Age of the universe (cosmic time at a=1 from a=0):
 *   t_0 = ∫_{0}^{1} da / (a * H(a))
 *
 * @param p  cosmological parameters
 * @return age of universe in Gyr
 */
double age_of_universe(const CosmologicalParams *p);

/**
 * Lookback time to a source at scale factor a:
 *   t_lookback = t_0 - t(a)
 */
double lookback_time(const CosmologicalParams *p, double a);

/**
 * Comoving distance to source at scale factor a:
 *   chi(a) = ∫_{a}^{1} da' / (a'^2 * H(a'))
 */
double comoving_distance(const CosmologicalParams *p, double a);

/**
 * Luminosity distance:
 *   d_L = (1+z) * d_M
 * where d_M is the transverse comoving distance.
 */
double luminosity_distance(const CosmologicalParams *p, double a);

/**
 * Angular diameter distance:
 *   d_A = d_M / (1+z)
 */
double angular_diameter_distance(const CosmologicalParams *p, double a);

/**
 * Critical density of the universe at present time:
 *   rho_c = 3 * H_0^2 / (8 * pi * G)
 */
double critical_density(const CosmologicalParams *p);

/* ---------------------------------------------------------------------------
 * L8: Inflation — de Sitter expansion
 * -------------------------------------------------------------------------*/

/**
 * de Sitter scale factor: a(t) = a_0 * exp(H * t)
 * Used as prototype for cosmic inflation model.
 *
 * @param a0    scale factor at t=0
 * @param H     Hubble constant during inflation (s^{-1})
 * @param t     proper time
 * @return a(t)
 */
double de_sitter_scale_factor(double a0, double H, double t);

/**
 * Number of e-folds of inflation:
 *   N = ln(a_final / a_initial) = H * Delta_t
 */
double inflation_e_folds(double a_initial, double a_final, double H, double dt);

#endif /* MINI_EINSTEIN_COSMOLOGY_H */
