#ifndef PERTURBATIONS_H
#define PERTURBATIONS_H

#include "cosmo_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Cosmological perturbation theory.
 *
 * Describes the growth of small density fluctuations delta = delta_rho/rho
 * in the linear regime (|delta| << 1). Nonlinear corrections extend to
 * spherical collapse and Press-Schechter formalism.
 *
 * Key equation (Newtonian gauge, sub-horizon):
 *   delta_ddot + 2*H*delta_dot - 4*pi*G*rho_m*delta = 0
 *
 * Reference: Dodelson & Schmidt Ch.7-8, Peebles (1980)
 * L3 - Perturbation theory
 * L6 - Linear growth factor
 * L8 - Nonlinear corrections
 * ============================================================================ */

/* === Linear growth factor D(a) === */

/* Compute the linear growth factor D(a) normalized to D(a=1) = 1.
 * D(a) satisfies: D'' + (3/a + E'/E)*D' - (3*Omega_m/(2*a^2*E^2))*D = 0
 * where prime = d/da.
 *
 * For matter-dominated flat universe: D(a) = a (exact).
 * For LCDM: numerically integrated.
 *
 * L6 - Linear growth factor (canonical)
 */
double perturbations_growth_factor_D(const CosmoModel *model, double a);

/* Compute the logarithmic growth rate f(a) = d(ln D)/d(ln a).
 * f(a) ~ Omega_m(a)^gamma with gamma ~ 0.55 for LCDM.
 *
 * This is measurable from redshift-space distortions (RSD).
 * L7 - Redshift-space distortions, growth rate measurement
 */
double perturbations_growth_rate_f(const CosmoModel *model, double a);

/* === Matter power spectrum evolution === */

/* Scale the linear matter power spectrum from initial conditions
 * to redshift z: P(k,z) = P(k,0) * [D(z)/D(0)]^2 (linear regime).
 *
 * L5 - Power spectrum time evolution
 */
double perturbations_power_spectrum_evolved(double P_k_ini, double D_ini,
                                             double D_z);

/* === RMS fluctuation sigma(R) === */

/* Compute sigma(R): the RMS matter fluctuation smoothed on scale R:
 * sigma^2(R) = (1/2pi^2) * integral P(k) * W^2(kR) * k^2 dk
 * where W(x) = 3*(sin(x) - x*cos(x))/x^3 (top-hat filter).
 *
 * sigma_8 = sigma(8 Mpc/h) is the canonical normalization of the
 * matter power spectrum.
 *
 * L6 - sigma_8 normalization (canonical)
 */
double perturbations_sigma_R(const CosmoModel *model, double R_Mpc_h);

/* === Spherical collapse model === */

/* Compute the critical linear overdensity delta_c for spherical collapse.
 * For Einstein-de Sitter: delta_c = 1.686 (constant).
 * For LCDM: delta_c(z) has weak cosmology dependence.
 *
 * L8 - Spherical collapse model, delta_c
 */
double perturbations_delta_c(const CosmoModel *model, double z);

/* === Press-Schechter mass function === */

/* Compute the halo mass function dn/dM (number density of halos per unit mass):
 * dn/dM = sqrt(2/pi) * (rho_m/M) * (delta_c/sigma) * |d(ln sigma)/dM|
 *        * exp(-delta_c^2/(2*sigma^2))
 *
 * L8 - Press-Schechter formalism
 */
double perturbations_mass_function_PS(const CosmoModel *model,
                                       double M_Msun_h, double z);

/* === Halo bias === */

/* Compute the linear halo bias b(M) for Press-Schechter:
 * b(M) = 1 + (delta_c^2/sigma^2 - 1)/delta_c
 *
 * L8 - Halo bias (peak-background split)
 */
double perturbations_halo_bias_PS(const CosmoModel *model,
                                   double M_Msun_h, double z);

/* === Nonlinear power spectrum (halofit) === */

/* Compute the nonlinear matter power spectrum using the Smith et al. (2003)
 * halofit prescription:
 *   Delta_NL(k) = f_NL(Delta_L(k))
 * where Delta(k) = k^3*P(k)/(2*pi^2) is the dimensionless power.
 *
 * L7 - Nonlinear power spectrum (halofit)
 * L8 - Smith et al. 2003 fitting functions
 */
double perturbations_nonlinear_power(const CosmoModel *model,
                                      double k_Mpc, double z);

/* === Baryon acoustic oscillations (BAO) === */

/* Compute the BAO wiggle component of the power spectrum.
 * P(k) = P_smooth(k) * [1 + O(k) * exp(-k^2*Sigma_nl^2/2)]
 * where O(k) is the oscillatory BAO feature.
 *
 * L7 - Baryon acoustic oscillations
 */
double perturbations_bao_oscillation(double k_Mpc, double r_drag_Mpc,
                                      double Sigma_nl);

/* === CMB temperature anisotropy spectrum === */

/* Compute the Sachs-Wolfe contribution to C_l^TT:
 * C_l^TT_SW = (2/pi) * integral k^2 dk * P_primordial(k) * j_l^2(k*D_rec)
 * where D_rec is the comoving distance to recombination.
 *
 * L7 - CMB temperature power spectrum (Sachs-Wolfe)
 */
double perturbations_cmb_Cl_TT_SW(const CosmoModel *model, int l);

#ifdef __cplusplus
}
#endif

#endif
