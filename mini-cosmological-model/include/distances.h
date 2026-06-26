#ifndef DISTANCES_H
#define DISTANCES_H

#include "cosmo_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Cosmological distance measures.
 *
 * In FLRW cosmology, there are multiple distinct distance definitions,
 * each measurable by different observational techniques.
 *
 * Core relation (Etherington 1933, reciprocity theorem):
 *   D_L = (1+z)^2 * D_A
 *
 * Reference: Hogg (1999) "Distance measures in cosmology"
 *            Dodelson & Schmidt Ch.2
 *
 * L1 - Distance definitions
 * L2 - Observational cosmology
 * L6 - Canonical distance ladder
 * ============================================================================ */

/* === Comoving line-of-sight distance ===
 * D_C(z) = D_H * integral_0^z [dz'/E(z')]
 * where D_H = c/H0 is the Hubble distance.
 *
 * This is the distance that remains constant between comoving observers
 * in an expanding universe.
 */
double distances_comoving_line_of_sight(const CosmoModel *model, double z);

/* === Transverse comoving distance ===
 * For a flat universe (Omega_k = 0): D_M = D_C
 * For Omega_k > 0 (spherical): D_M = D_H/sqrt(Omega_k) * sinh(sqrt(Omega_k)*D_C/D_H)
 * For Omega_k < 0 (hyperbolic): D_M = D_H/sqrt(|Omega_k|) * sin(sqrt(|Omega_k|)*D_C/D_H)
 *
 * D_M is used in the angular diameter distance and luminosity distance.
 * L6 - Transverse comoving distance
 */
double distances_transverse_comoving(const CosmoModel *model, double z);

/* === Angular diameter distance ===
 * D_A(z) = D_M(z) / (1+z)
 * = physical_size / angular_size
 *
 * At high z, D_A actually decreases with z (objects appear larger),
 * a key test of the FLRW metric.
 * L6 - Angular diameter distance
 */
double distances_angular_diameter(const CosmoModel *model, double z);

/* === Luminosity distance ===
 * D_L(z) = (1+z) * D_M(z)
 * Flux = Luminosity / (4*pi*D_L^2)
 *
 * This is the distance measured from supernova standard candles,
 * providing the primary evidence for cosmic acceleration.
 * L6 - Luminosity distance, SNIa cosmology
 */
double distances_luminosity(const CosmoModel *model, double z);

/* === Comoving volume ===
 * dV_c/dzdOmega = D_H * (1+z)^2 * D_A^2 / E(z)
 * V_c(z) = integral dV_c
 *
 * Used for number counts of galaxies, clusters, etc.
 * L6 - Comoving volume
 */
double distances_comoving_volume_element(const CosmoModel *model, double z);
double distances_comoving_volume(const CosmoModel *model, double z);

/* === Distance modulus ===
 * mu(z) = m - M = 5*log10(D_L / 10 pc)
 *        = 5*log10(D_L/Mpc) + 25
 *
 * This is the difference between apparent and absolute magnitude.
 * L6 - Distance modulus (standard candles)
 */
double distances_distance_modulus(const CosmoModel *model, double z);

/* === Light travel distance ===
 * D_LT(z) = c * lookback_time(z)
 * Not directly observable but useful for conceptual purposes.
 * L1 - Light travel distance
 */
double distances_light_travel(const CosmoModel *model, double z);

/* === Lookback time ===
 * t_L(z) = t_H * integral_0^z [dz'/((1+z')*E(z'))]
 * where t_H = 1/H0 is the Hubble time.
 * L1 - Lookback time
 */
double distances_lookback_time_Gyr(const CosmoModel *model, double z);

/* === Age of universe at redshift z ===
 * t_age(z) = t_H * integral_z^inf [dz'/((1+z')*E(z'))]
 * L6 - Age of the universe
 */
double distances_age_Gyr(const CosmoModel *model, double z);

/* === Particle horizon ===
 * d_p(z) = c * integral_0^t dt'/a(t') = D_H * integral_z^inf [dz'/E(z')]
 *
 * The maximum comoving distance from which light could have reached
 * an observer since the Big Bang.
 * L6 - Particle horizon (canonical)
 */
double distances_particle_horizon(const CosmoModel *model, double z);

/* === Event horizon ===
 * d_e(z) = c * a(t) * integral_t^inf dt'/a(t')
 *        = D_H * integral_{-1}^z [dz'/E(z')]
 *
 * The maximum comoving distance from which light emitted now will
 * ever reach the observer. Non-zero only for accelerating universes.
 * L6 - Event horizon (canonical)
 */
double distances_event_horizon(const CosmoModel *model);

/* === Gaussian integration kernel for distance integrals ===
 *
 * Internal helper: computes the integrand for comoving distance:
 * f(z) = 1/E(z)
 * L5 - Numerical integration kernel
 */
double distances_integrand_E_inv(const CosmoModel *model, double z);

/* === Adaptively integrate a function f(z) from z1 to z2 ===
 * Uses 5-point Gauss-Legendre quadrature on sub-intervals,
 * refining until relative error < tol.
 *
 * L5 - Adaptive numerical quadrature
 */
double distances_adaptive_integrate(const CosmoModel *model,
                                     double z1, double z2,
                                     double (*f)(const CosmoModel*, double),
                                     double tol, int max_refine);

#ifdef __cplusplus
}
#endif

#endif
