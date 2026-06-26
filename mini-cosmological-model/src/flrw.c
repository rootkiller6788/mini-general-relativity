/** @file flrw.c
 *  @brief FLRW metric, geometry, and coordinate computations.
 *
 *  The Friedmann-Lemaitre-Robertson-Walker (FLRW) metric is the most general
 *  metric describing a homogeneous and isotropic universe.
 *
 *  Line element (spherical comoving coordinates):
 *    ds^2 = -c^2 dt^2 + a(t)^2 [dr^2/(1-k*r^2) + r^2(dtheta^2 + sin^2(theta) dphi^2)]
 *
 *  Reference: Wald (1984) Ch.5, Carroll (2004) Ch.8
 *             Dodelson & Schmidt (2020) Ch.2
 */

#include "flrw.h"
#include "friedmann.h"
#include <math.h>

/* === Metric components === */

/* g_tt = -c^2 (time-time component) */
double flrw_g_tt(double c)
{
    return -c * c;
}

/* g_rr = a^2 / (1 - k * r^2) (radial-radial component) */
double flrw_g_rr(double a, double k, double r)
{
    double denom = 1.0 - k * r * r;
    if (fabs(denom) < 1.0e-15) {
        /* Singularity at r = 1/sqrt(k) for k > 0 (closed universe pole) */
        return 1.0 / 1.0e-15;
    }
    return a * a / denom;
}

/* g_thetatheta = a^2 * r^2 (polar-polar component) */
double flrw_g_thetatheta(double a, double r)
{
    return a * a * r * r;
}

/* g_phiphi = a^2 * r^2 * sin^2(theta) (azimuthal-azimuthal component) */
double flrw_g_phiphi(double a, double r, double theta)
{
    return a * a * r * r * sin(theta) * sin(theta);
}

/* === Ricci scalar === */

/* Ricci scalar for FLRW metric:
 * R = 6/c^2 * [a_ddot/a + (a_dot/a)^2 + k*c^2/a^2]
 *
 * This is derived from the Christoffel symbols of the FLRW metric.
 * For a flat universe (k=0) with Lambda dominance:
 *   R = 12 * H^2 / c^2  (de Sitter)
 * For flat matter-dominated:
 *   R = 3 * H0^2 / (c^2 * a^3)  (vanishes at late times)
 */
double flrw_ricci_scalar(double a, double a_dot, double a_ddot,
                          double k, double c)
{
    double H2 = a_dot * a_dot / (a * a);
    double adot_over_a = a_ddot / a;
    double curvature_term = k * c * c / (a * a);

    return 6.0 / (c * c) * (adot_over_a + H2 + curvature_term);
}

/* === Curvature parameter === */

/* Compute k(today) from Omega_k0.
 * Omega_k = -k*c^2/(a^2*H^2) => k = -Omega_k * a^2 * H^2 / c^2
 * At present day: k = -Omega_k0 * H0^2 / c^2 (with a_0 = 1) */
double flrw_k_from_omega_k(double Omega_k0, double H0_si, double c)
{
    return -Omega_k0 * H0_si * H0_si / (c * c);
}

/* === Conformal time === */

/* Conformal time eta is defined by d(eta) = dt / a(t).
 * In terms of redshift: eta(z) = integral_0^z dz' / H(z')
 * We compute eta at redshift z (eta = 0 at Big Bang a=0).
 *
 * The conformal time to recombination (z_recomb ~ 1100) determines
 * the sound horizon scale that sets CMB acoustic peak positions. */
double flrw_conformal_time(const CosmoModel *model, double z)
{
    if (!model) return -1.0;

    int n = model->n_z_steps;
    if (n < 10) n = 1000;

    double z_max = 1.0e6;
    double log_z_min = log(1.0 + z);
    double log_z_max = log(1.0 + z_max);
    double d_logz = (log_z_max - log_z_min) / n;

    double integral = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double log_z_i = log_z_min + i * d_logz;
        double log_z_mid = log_z_i + 0.5 * d_logz;
        double log_z_next = log_z_i + d_logz;

        double z_i   = exp(log_z_i)   - 1.0;
        double z_mid = exp(log_z_mid) - 1.0;
        double z_next = exp(log_z_next) - 1.0;

        double H_i   = friedmann_H_z_si(model, z_i);
        double H_mid = friedmann_H_z_si(model, z_mid);
        double H_next = friedmann_H_z_si(model, z_next);

        double f_i   = (H_i   > 0.0) ? 1.0 / H_i   : 0.0;
        double f_mid = (H_mid > 0.0) ? 1.0 / H_mid : 0.0;
        double f_next = (H_next > 0.0) ? 1.0 / H_next : 0.0;

        double dz = z_next - z_i;
        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return integral;
}

/* === Proper distance === */

/* The proper (physical) distance between two comoving points
 * scales with the scale factor: d_proper = a * r_comoving.
 * Hubble's law: v = H * d_proper for nearby objects. */
double flrw_proper_distance(double a, double r_comoving)
{
    return a * r_comoving;
}

/* === Volume expansion === */

/* The proper volume of a comoving region expands as dV/dt = 3*H*V.
 * This is the kinematic expression of cosmic expansion. */
double flrw_volume_expansion_rate(double H_si, double V_proper)
{
    return 3.0 * H_si * V_proper;
}

/* === Redshift-scale factor relations === */

/* 1 + z = a_obs / a_emit.
 * For observations today: a_obs = 1, so 1 + z = 1 / a_emit. */
double flrw_redshift_from_scale(double a_emit, double a_obs)
{
    if (a_emit <= 0.0) return -1.0;
    return a_obs / a_emit - 1.0;
}

double flrw_scale_from_redshift(double z)
{
    if (z < -1.0) return -1.0;
    return 1.0 / (1.0 + z);
}

/* === Cosmological distances === */

/* Angular diameter distance: D_A = a_emit * r_comoving = r_comoving / (1+z).
 * For a standard ruler of physical size l subtending angle Delta_theta:
 *   D_A = l / Delta_theta */
double flrw_angular_diameter_distance(double r_comoving, double z)
{
    return r_comoving / (1.0 + z);
}

/* Luminosity distance: D_L = (1+z) * r_comoving.
 * For a standard candle of luminosity L and observed flux F:
 *   F = L / (4*pi*D_L^2)
 * Etherington reciprocity: D_L = (1+z)^2 * D_A */
double flrw_luminosity_distance(double r_comoving, double z)
{
    return r_comoving * (1.0 + z);
}

/* Comoving volume element per unit solid angle per unit redshift:
 * dV/(dOmega dz) = D_H * (1+z)^2 * D_A^2 / E(z) */
double flrw_comoving_volume_element(double z, double D_A, double E_z,
                                     double D_H, double dOmega, double dz)
{
    double dV = D_H * (1.0 + z) * (1.0 + z) * D_A * D_A / E_z;
    return dV * dOmega * dz;
}

/* === Cosmographic parameters === */

/* Deceleration parameter q.
 * q = -a*a_ddot / a_dot^2 = -a_ddot / (a * H^2)
 * q > 0: decelerating expansion (matter/radiation dominated)
 * q = 0: coasting
 * q < 0: accelerating expansion (dark energy dominated)
 *
 * For Lambda-CDM today: q_0 = Omega_m/2 - Omega_Lambda ~ -0.53 */
double flrw_deceleration_parameter(double a, double a_dot, double a_ddot)
{
    if (fabs(a_dot) < 1.0e-30) return 0.0;
    return -a * a_ddot / (a_dot * a_dot);
}

/* Jerk parameter j (third derivative of scale factor).
 * j = a^2*a_dddot / a_dot^3
 * For flat LCDM: j_0 = 1 exactly (independent of parameters).
 *
 * The jerk provides a model-independent test of Lambda-CDM:
 * any deviation from j=1 indicates dynamics beyond a cosmological constant. */
double flrw_jerk_parameter(double a, double a_dot, double a_ddot,
                            double a_dddot)
{
    (void)a_ddot;  /* Not needed for flat LCDM (j_0 = 1), but kept for API completeness */
    if (fabs(a_dot) < 1.0e-30) return 0.0;
    return a * a * a_dddot / (a_dot * a_dot * a_dot);
}

/* Snap parameter s (fourth derivative of scale factor).
 * s = a^3*a_ddddot / a_dot^4
 * For LCDM: s_0 = 1 - (9/2)*Omega_m, which is a probe of matter density. */
double flrw_snap_parameter(double a, double a_dot, double a_ddot,
                            double a_dddot, double a_ddddot)
{
    (void)a_ddot;   /* Not needed for flat LCDM, kept for API completeness */
    (void)a_dddot;  /* Not needed for flat LCDM, kept for API completeness */
    if (fabs(a_dot) < 1.0e-30) return 0.0;
    return a * a * a * a_ddddot / (a_dot * a_dot * a_dot * a_dot);
}
