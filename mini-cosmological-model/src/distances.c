/** @file distances.c
 *  @brief Cosmological distance measures: comoving, angular diameter,
 *         luminosity, volume, horizons, and distance modulus.
 *
 *  In FLRW cosmology there are multiple distinct distance measures,
 *  each relevant to different observational probes.
 *
 *  Etherington reciprocity theorem (1933):
 *    D_L = (1+z)^2 * D_A
 *  This is a purely geometric result independent of the field equations.
 *
 *  Reference: Hogg (1999) "Distance measures in cosmology"
 *             Dodelson & Schmidt (2020) Ch.2
 *             Weinberg (2008) "Cosmology"
 */

#include "distances.h"
#include "friedmann.h"
#include <math.h>

/* === Comoving line-of-sight distance ===
 * D_C(z) = D_H * integral_0^z dz'/E(z')
 * where D_H = c/H0 is the Hubble distance (~ 4.45 Gpc for Planck 2018).
 */
double distances_comoving_line_of_sight(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;

    double D_H = COSMO_C_LIGHT / model->H0_si;  /* Hubble distance [m] */
    int n = model->n_z_steps;
    if (n < 10) n = 1000;

    double dz = z / n;
    double integral = 0.0;
    int i;

    /* Simpson integration over [0, z] */
    for (i = 0; i < n; i++) {
        double z_i = i * dz;
        double z_mid = z_i + 0.5 * dz;
        double z_next = z_i + dz;

        double f_i   = 1.0 / friedmann_E_z(model, z_i);
        double f_mid = 1.0 / friedmann_E_z(model, z_mid);
        double f_next = 1.0 / friedmann_E_z(model, z_next);

        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return D_H * integral;
}

/* === Transverse comoving distance ===
 * D_M(z) for curved universes:
 *   Omega_k = 0: D_M = D_C
 *   Omega_k > 0: D_M = D_H/sqrt(Omega_k) * sinh(sqrt(Omega_k) * D_C/D_H)
 *   Omega_k < 0: D_M = D_H/sqrt(|Omega_k|) * sin(sqrt(|Omega_k|) * D_C/D_H)
 */
double distances_transverse_comoving(const CosmoModel *model, double z)
{
    if (!model) return -1.0;

    double D_C = distances_comoving_line_of_sight(model, z);
    if (D_C < 0.0) return -1.0;

    double D_H = COSMO_C_LIGHT / model->H0_si;
    double Omega_k = model->Omega_k;

    if (fabs(Omega_k) < 1.0e-8) {
        /* Flat universe: D_M = D_C */
        return D_C;
    }

    double sqrt_abs_ok = sqrt(fabs(Omega_k));
    double x = sqrt_abs_ok * D_C / D_H;

    if (Omega_k > 0.0) {
        /* Spherical: D_M = (D_H/sqrt(Omega_k)) * sinh(x) */
        return D_H / sqrt_abs_ok * sinh(x);
    } else {
        /* Hyperbolic: D_M = (D_H/sqrt(|Omega_k|)) * sin(x) */
        return D_H / sqrt_abs_ok * sin(x);
    }
}

/* === Angular diameter distance ===
 * D_A(z) = D_M(z) / (1+z)
 * = physical_transverse_size / angular_size
 *
 * D_A peaks at z ~ 1.6 for LCDM (objects beyond this appear larger).
 */
double distances_angular_diameter(const CosmoModel *model, double z)
{
    if (!model || z < -0.999) return -1.0;
    double D_M = distances_transverse_comoving(model, z);
    if (D_M < 0.0) return -1.0;
    return D_M / (1.0 + z);
}

/* === Luminosity distance ===
 * D_L(z) = (1+z) * D_M(z) = (1+z)^2 * D_A(z)
 * Flux = L / (4*pi*D_L^2)
 *
 * SNIa Hubble diagram uses D_L to constrain Omega_m and Omega_Lambda.
 */
double distances_luminosity(const CosmoModel *model, double z)
{
    if (!model || z < -0.999) return -1.0;
    double D_M = distances_transverse_comoving(model, z);
    if (D_M < 0.0) return -1.0;
    return D_M * (1.0 + z);
}

/* === Comoving volume element ===
 * dV_c/(dOmega dz) = D_H * (1+z)^2 * D_A^2 / E(z)
 */
double distances_comoving_volume_element(const CosmoModel *model, double z)
{
    if (!model || z < -0.999) return -1.0;
    double D_A = distances_angular_diameter(model, z);
    double E_z = friedmann_E_z(model, z);
    double D_H = COSMO_C_LIGHT / model->H0_si;

    if (D_A < 0.0 || E_z <= 0.0) return 0.0;
    return D_H * (1.0 + z) * (1.0 + z) * D_A * D_A / E_z;
}

/* === Total comoving volume out to redshift z ===
 * V_c(z) = (4*pi/3) * D_M^3  for flat universe (exact).
 * For curved universes, we integrate the volume element.
 */
double distances_comoving_volume(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;

    if (fabs(model->Omega_k) < 1.0e-8) {
        /* Flat universe: V_c = (4*pi/3) * D_C^3 */
        double D_C = distances_comoving_line_of_sight(model, z);
        return (4.0 * M_PI / 3.0) * D_C * D_C * D_C;
    }

    /* Curved universe: integrate numerically */
    int n = model->n_z_steps;
    if (n < 10) n = 1000;
    double dz = z / n;
    double V = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double z_i = i * dz;
        double z_mid = z_i + 0.5 * dz;
        double z_next = z_i + dz;

        double dV_i   = distances_comoving_volume_element(model, z_i);
        double dV_mid = distances_comoving_volume_element(model, z_mid);
        double dV_next = distances_comoving_volume_element(model, z_next);

        V += (dz / 6.0) * (dV_i + 4.0 * dV_mid + dV_next) * (4.0 * M_PI);
    }

    return V;
}

/* === Distance modulus ===
 * mu = m - M = 5*log10(D_L/pc)
 * D_L in pc = D_L [m] / (3.085677581e16 m/pc)
 * mu = 5*log10(D_L [m]) - 5*log10(3.085677581e16)
 *    = 5*log10(D_L [m]) - 82.441...
 *
 * More commonly: mu = 5*log10(D_L [Mpc]) + 25
 */
double distances_distance_modulus(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;
    double D_L = distances_luminosity(model, z);
    if (D_L <= 0.0) return -1.0;

    /* D_L in Mpc */
    double D_L_Mpc = D_L / COSMO_MPC;
    return 5.0 * log10(D_L_Mpc) + 25.0;
}

/* === Light travel distance ===
 * D_LT = c * (t_0 - t_z) = c * lookback_time
 */
double distances_light_travel(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;

    double age0 = friedmann_cosmic_time(model, 0.0);
    double age_z = friedmann_cosmic_time(model, z);

    if (age0 < 0.0 || age_z < 0.0) return -1.0;
    return COSMO_C_LIGHT * (age0 - age_z);
}

/* === Lookback time in Gyr === */
double distances_lookback_time_Gyr(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;

    double age0 = friedmann_cosmic_time(model, 0.0);
    double age_z = friedmann_cosmic_time(model, z);

    if (age0 < 0.0 || age_z < 0.0) return -1.0;
    return (age0 - age_z) / COSMO_GYR_SEC;
}

/* === Age of universe at redshift z in Gyr === */
double distances_age_Gyr(const CosmoModel *model, double z)
{
    if (!model) return -1.0;
    double t = friedmann_cosmic_time(model, z);
    if (t < 0.0) return -1.0;
    return t / COSMO_GYR_SEC;
}


/* === Particle horizon ===
 * d_p(z) = c * integral_0^t dt'/a(t')
 *        = D_H * integral_z^inf dz'/E(z')
 *
 * The particle horizon is the maximum comoving distance from which
 * information could have propagated to the observer since the Big Bang.
 * It solves the horizon problem without inflation.
 */
double distances_particle_horizon(const CosmoModel *model, double z)
{
    if (!model) return -1.0;

    double D_H = COSMO_C_LIGHT / model->H0_si;
    int n = model->n_z_steps;
    if (n < 10) n = 1000;

    double z_max = 1.0e7;
    double log_zp1_min = log(1.0 + z);
    double log_zp1_max = log(1.0 + z_max);
    double d_logzp1 = (log_zp1_max - log_zp1_min) / n;

    double integral = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double log_zp1_i = log_zp1_min + i * d_logzp1;
        double log_zp1_mid = log_zp1_i + 0.5 * d_logzp1;
        double log_zp1_next = log_zp1_i + d_logzp1;

        double zp1_i   = exp(log_zp1_i);
        double zp1_mid = exp(log_zp1_mid);
        double zp1_next = exp(log_zp1_next);

        double z_i   = zp1_i   - 1.0;
        double z_mid = zp1_mid - 1.0;
        double z_next = zp1_next - 1.0;

        double f_i   = 1.0 / friedmann_E_z(model, z_i);
        double f_mid = 1.0 / friedmann_E_z(model, z_mid);
        double f_next = 1.0 / friedmann_E_z(model, z_next);

        double dz = z_next - z_i;
        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return D_H * integral;
}

/* === Event horizon ===
 * d_e(z) = c * a(t) * integral_t^inf dt'/a(t')
 *        = D_H * integral_{-1}^z dz'/E(z')
 *
 * The event horizon exists only in accelerating universes (w_eff < -1/3).
 * It represents the maximum comoving distance from which light emitted
 * NOW will ever reach the observer.
 *
 * For LCDM: d_e(0) ~ 16.5 Glyr (comoving), ~ 5.0 Gpc (proper).
 */
double distances_event_horizon(const CosmoModel *model)
{
    if (!model) return -1.0;

    /* Check if event horizon exists (requires late-time acceleration) */
    if (model->Omega_Lambda <= 0.0) {
        return INFINITY;  /* No event horizon in decelerating universe */
    }

    double D_H = COSMO_C_LIGHT / model->H0_si;
    int n = model->n_z_steps;
    if (n < 10) n = 1000;

    /* Integrate from z=0 to z=-1 (infinite future) */
    double z_min = -0.9999;  /* a very large, avoid singularity at z=-1 */
    double dz = (0.0 - z_min) / n;

    double integral = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double z_i = z_min + i * dz;
        double z_mid = z_i + 0.5 * dz;
        double z_next = z_i + dz;

        double f_i   = 1.0 / friedmann_E_z(model, z_i);
        double f_mid = 1.0 / friedmann_E_z(model, z_mid);
        double f_next = 1.0 / friedmann_E_z(model, z_next);

        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return D_H * integral;
}

/* === Integrand for E(z)^{-1} === */
double distances_integrand_E_inv(const CosmoModel *model, double z)
{
    if (!model) return 0.0;
    double E = friedmann_E_z(model, z);
    if (E <= 0.0) return 0.0;
    return 1.0 / E;
}

/* === Adaptive Gauss-Legendre integration ===
 * Integrates f(z) from z1 to z2 using 5-point Gauss-Legendre quadrature
 * with adaptive sub-interval refinement.
 *
 * 5-point GL nodes and weights on [-1, 1]:
 *   x_i = 0, +/-sqrt(5-2*sqrt(10/7))/3, +/-sqrt(5+2*sqrt(10/7))/3
 *   w_i = 128/225, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900
 */
double distances_adaptive_integrate(const CosmoModel *model,
                                     double z1, double z2,
                                     double (*f)(const CosmoModel*, double),
                                     double tol, int max_refine)
{
    if (!model || !f || z1 >= z2 || max_refine < 1) return 0.0;

    /* Pre-computed 5-point Gauss-Legendre on [-1,1] */
    static const double x_gl[5] = {
        0.0,
        -0.5384693101056831,
         0.5384693101056831,
        -0.9061798459386640,
         0.9061798459386640
    };
    static const double w_gl[5] = {
        0.5688888888888889,
        0.4786286704993665,
        0.4786286704993665,
        0.2369268850561891,
        0.2369268850561891
    };

    double mid = 0.5 * (z1 + z2);
    double half = 0.5 * (z2 - z1);
    double integral = 0.0;
    int j;

    for (j = 0; j < 5; j++) {
        double z_j = mid + half * x_gl[j];
        integral += w_gl[j] * f(model, z_j);
    }
    integral *= half;

    /* For adaptive refinement, we split the interval and check convergence.
     * This simplified version does a single-level refinement check. */
    if (max_refine > 1 && half > tol) {
        double left  = distances_adaptive_integrate(model, z1, mid, f, tol, max_refine - 1);
        double right = distances_adaptive_integrate(model, mid, z2, f, tol, max_refine - 1);
        double refined = left + right;

        /* Accept refinement if it differs significantly */
        if (fabs(refined - integral) > tol * fabs(refined)) {
            return refined;
        }
    }

    return integral;
}
