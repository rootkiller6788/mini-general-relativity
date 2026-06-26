#ifndef FLRW_H
#define FLRW_H

#include "cosmo_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* === FLRW metric line element ===
 * ds^2 = -c^2 dt^2 + a(t)^2 [dr^2/(1 - k r^2) + r^2 (dtheta^2 + sin^2 theta dphi^2)]
 *
 * where a(t) is the scale factor, k in {-1,0,+1} is the curvature parameter,
 * and (r, theta, phi) are comoving spherical coordinates.
 *
 * Reference: Wald Ch.5, Carroll Ch.8, Dodelson Ch.2
 * L3 - Mathematical Structures: Robertson-Walker metric, Riemann geometry
 */

double flrw_g_tt(double c);
double flrw_g_rr(double a, double k, double r);
double flrw_g_thetatheta(double a, double r);
double flrw_g_phiphi(double a, double r, double theta);

/* Compute the Ricci scalar for FLRW metric:
 * R = 6/c^2 * (a_ddot/a + (a_dot/a)^2 + k*c^2/a^2)
 * L3 - Riemann curvature, Ricci scalar
 */
double flrw_ricci_scalar(double a, double a_dot, double a_ddot,
                          double k, double c);

/* Determine the curvature parameter k from Omega_k:
 * Omega_k = -k*c^2/(a^2 * H^2)
 * => k = -Omega_k * a^2 * H^2 / c^2
 * At present day (a=1): k = -Omega_k0 * H0^2 / c^2
 * L1 - Curvature parameter
 */
double flrw_k_from_omega_k(double Omega_k0, double H0_si, double c);

/* Compute the conformal time eta:
 * eta = integral_0^t dt'/a(t')
 * d(eta) = dt / a(t)
 * L3 - Conformal time coordinate
 */
double flrw_conformal_time(const CosmoModel *model, double z);

/* Compute proper distance d_proper = a(t) * r_comoving
 * L1 - Proper distance vs comoving distance
 */
double flrw_proper_distance(double a, double r_comoving);

/* Compute the expansion rate of a comoving volume element:
 * dV_proper / dt = 3 * H * V_proper
 * L2 - Cosmic expansion, volume dilution
 */
double flrw_volume_expansion_rate(double H_si, double V_proper);

/* Compute the redshift of a photon emitted at scale factor a_emit:
 * 1 + z = a_obs / a_emit
 * L1 - Cosmological redshift
 */
double flrw_redshift_from_scale(double a_emit, double a_obs);

/* Inverse: scale factor from redshift */
double flrw_scale_from_redshift(double z);

/* Compute the angular diameter distance:
 * D_A = a_emit * r_comoving = r_comoving / (1+z)
 * L6 - Angular diameter distance (canonical)
 */
double flrw_angular_diameter_distance(double r_comoving, double z);

/* Compute the luminosity distance:
 * D_L = (1+z)^2 * D_A = (1+z) * r_comoving
 * L6 - Luminosity distance (canonical), Etherington reciprocity theorem
 */
double flrw_luminosity_distance(double r_comoving, double z);

/* Compute the comoving volume element:
 * dV_comoving = D_H * (1+z)^2 * D_A^2 / E(z) * dOmega * dz
 * where D_H = c/H0 is the Hubble distance
 * L6 - Comoving volume element
 */
double flrw_comoving_volume_element(double z, double D_A, double E_z,
                                     double D_H, double dOmega, double dz);

/* === Spacetime classification === */

/* Compute the deceleration parameter q:
 * q = -a*a_ddot / a_dot^2
 * q > 0: decelerating; q < 0: accelerating
 * L1 - Deceleration parameter
 */
double flrw_deceleration_parameter(double a, double a_dot, double a_ddot);

/* Compute the jerk parameter j (cosmographic expansion):
 * j = a^2*a_dddot / a_dot^3
 * For LCDM: j(t0) = 1
 * L1 - Jerk parameter (cosmography)
 */
double flrw_jerk_parameter(double a, double a_dot, double a_ddot, double a_dddot);

/* Compute the snap parameter s (cosmographic expansion):
 * s = a^3*a_ddddot / a_dot^4
 * L8 - Cosmographic expansion series
 */
double flrw_snap_parameter(double a, double a_dot, double a_ddot,
                            double a_dddot, double a_ddddot);

#ifdef __cplusplus
}
#endif

#endif
