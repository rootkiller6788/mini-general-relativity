#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H

#include "cosmo_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Primordial and matter power spectrum.
 *
 * The primordial power spectrum from inflation:
 *   P_R(k) = A_s * (k/k_pivot)^{n_s - 1}
 *
 * The matter power spectrum at redshift z:
 *   P_m(k,z) = T^2(k) * [D(z)/D(0)]^2 * P_R(k) * (transfer function)
 *
 * Reference: Dodelson & Schmidt Ch.7-8, Eisenstein & Hu (1998)
 * L1 - Primordial power spectrum
 * L3 - Transfer function
 * L6 - LCDM power spectrum
 * ============================================================================ */

/* === Primordial power spectrum === */

/* Compute the dimensionless primordial curvature power spectrum:
 * P_R(k) = A_s * (k/k_pivot)^{n_s - 1}
 * Delta_R^2(k) = k^3 * P_R(k) / (2*pi^2)
 *
 * L1 - Primordial power spectrum definition
 */
double ps_primordial_dimensionless(double k_Mpc, double A_s,
                                    double n_s, double k_pivot);

/* Compute P_R(k) with running spectral index:
 * P_R(k) = A_s * (k/k_pivot)^{n_s - 1 + 0.5*alpha_s*ln(k/k_pivot)}
 *
 * L8 - Running of the spectral index
 */
double ps_primordial_running(double k_Mpc, double A_s, double n_s,
                              double alpha_s, double k_pivot);

/* === Transfer function === */

/* BBKS transfer function (Bardeen, Bond, Kaiser, Szalay 1986):
 * T_BBKS(q) = ln(1+2.34q)/(2.34q)
 *           * [1 + 3.89q + (16.1q)^2 + (5.46q)^3 + (6.71q)^4]^{-1/4}
 * where q = k/(Omega_m * h^2) Mpc^{-1}
 *
 * This is valid for CDM with negligible baryon fraction.
 * L6 - BBKS transfer function
 */
double ps_transfer_BBKS(double k_Mpc, double Omega_m, double h);

/* Eisenstein & Hu (1998) transfer function with baryon acoustic features:
 * T_EH(k) = (Omega_b/Omega_m)*T_b(k) + (Omega_c/Omega_m)*T_cdm(k)
 *
 * Includes BAO wiggles from pre-recombination sound waves.
 * L6 - Eisenstein & Hu transfer function
 * L7 - Baryon acoustic oscillations
 */
double ps_transfer_EH(double k_Mpc, const CosmoModel *model);

/* === Matter power spectrum === */

/* Compute the linear matter power spectrum at z=0:
 * P_lin(k) = A_s * (k/k_pivot)^{n_s-1} * T^2(k) * (k/H0^4)
 * scaled appropriately to return [Mpc^3/h^3].
 *
 * L6 - Linear matter power spectrum
 */
double ps_matter_power_linear(double k_Mpc, const CosmoModel *model);

/* Compute the dimensionless matter power spectrum:
 * Delta_m^2(k) = k^3 * P_m(k) / (2*pi^2)
 *
 * L1 - Dimensionless power spectrum
 */
double ps_matter_dimensionless(double k_Mpc, const CosmoModel *model);

/* === Effective spectral index and running === */

/* Compute the effective spectral index n_eff(k) = d(ln P)/d(ln k)
 * using centered finite differences.
 *
 * L8 - Effective spectral index
 */
double ps_effective_n(double k_Mpc, const CosmoModel *model, double dk_frac);

/* === Correlation function === */

/* Fourier transform of the power spectrum to real space:
 * xi(r) = (1/2pi^2) * integral_0^inf [k^2 * P(k) * sin(kr)/(kr) dk]
 *
 * The correlation function describes the excess probability of
 * finding a pair of galaxies separated by distance r.
 * L6 - Two-point correlation function
 */
double ps_correlation_function(double r_Mpc_h, const CosmoModel *model);

/* === Window function integrals === */

/* Top-hat window function in Fourier space:
 * W_th(x) = 3*(sin(x) - x*cos(x))/x^3
 * L5 - Top-hat filter kernel
 */
double ps_tophat_window(double kR);

/* Gaussian window function in Fourier space:
 * W_gauss(x) = exp(-x^2/2)
 * L5 - Gaussian filter kernel
 */
double ps_gaussian_window(double kR);

/* === Power spectrum convolution === */

/* Convolve the power spectrum with a window function:
 * sigma_R^2 = integral k^2 dk/(2pi^2) * P(k) * |W(kR)|^2
 * L6 - Variance on scale R
 */
double ps_variance_R(double R_Mpc_h, const CosmoModel *model,
                      double (*window)(double));

/* === Nonlinear power spectrum fitting (halofit) === */

/* Compute the nonlinear matter power spectrum using the Halofit model
 * (Smith et al. 2003, Takahashi et al. 2012 revision).
 *
 * P_NL(k) = f_HF(P_lin(k), sigma(R), n_eff)
 *
 * L7 - Halofit nonlinear corrections
 * L8 - Smith et al. 2003 / Takahashi et al. 2012
 */
double ps_nonlinear_halofit(double k_Mpc, double z, const CosmoModel *model);

/* === Power spectrum at arbitrary redshift === */

/* Evolve the linear power spectrum to redshift z:
 * P_lin(k, z) = P_lin(k, 0) * [D(z)/D(0)]^2
 *
 * L5 - Power spectrum redshift evolution
 */
double ps_power_at_z(double k_Mpc, double z, const CosmoModel *model);

#ifdef __cplusplus
}
#endif

#endif
