#ifndef FRIEDMANN_H
#define FRIEDMANN_H

#include "cosmo_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Friedmann equations: the fundamental ODEs governing cosmic expansion.
 *
 * First Friedmann equation (energy constraint):
 *   H^2 = (a_dot/a)^2 = (8*pi*G/3)*rho - k*c^2/a^2 + Lambda*c^2/3
 *
 * Second Friedmann equation (acceleration equation):
 *   a_ddot/a = -(4*pi*G/3)*(rho + 3*P/c^2) + Lambda*c^2/3
 *
 * Continuity equation (energy conservation for each perfect fluid):
 *   rho_dot_i + 3*H*(rho_i + P_i/c^2) = 0
 *
 * Equation of state: P_i = w_i * rho_i * c^2
 *   => rho_i(a) = rho_i0 * a^{-3(1+w_i)}
 *
 * Reference: Dodelson & Schmidt Ch.2, Wald Ch.5, Carroll Ch.8
 * L4 - Fundamental Laws: Friedmann equations
 * ============================================================================ */

/* === Density evolution === */

/* Compute the energy density for a species with equation of state w at scale a.
 * rho(a) = rho0 * a^{-3*(1+w)}
 * L4 - Continuity equation solution
 */
double friedmann_rho_a(double rho0, double a, double w);

/* Compute Omega_i(a) for species i:
 * Omega_i(a) = rho_i(a) / rho_crit(a)
 * = Omega_i0 * a^{-3*(1+w_i)} * (H0/H(a))^2
 * L1 - Density parameter evolution
 */
double friedmann_omega_a(double Omega_i0, double a, double w,
                          double H0_si, double H_si);

/* === Hubble parameter === */

/* Compute E(z) = H(z)/H0 for a multi-component universe:
 * E^2(z) = Omega_r0*(1+z)^4 + Omega_m0*(1+z)^3
 *        + Omega_k0*(1+z)^2 + Omega_Lambda0
 *        + Omega_de0*f_de(z)
 * L4 - Friedmann equation (first), Hubble parameter
 */
double friedmann_E_z(const CosmoModel *model, double z);

/* Compute H(z) [km/s/Mpc] at redshift z */
double friedmann_H_z(const CosmoModel *model, double z);

/* Compute H(z) [s^{-1}] at redshift z */
double friedmann_H_z_si(const CosmoModel *model, double z);

/* === Scale factor evolution ODE system === */

/* Right-hand side of the Friedmann ODE: da/dt = a * H(a).
 * For use with numerical integrators.
 * L5 - ODE formulation for numerical integration
 */
double friedmann_dadt(double a, double H_a_si);

/* Compute dH/da for the scale-factor evolution:
 * dH/da = -3*H/(2*a) * [1 + w_eff(a)]
 * where w_eff(a) = P_tot/(rho_tot*c^2)
 * L4 - Acceleration equation
 */
double friedmann_dH_da(const CosmoModel *model, double a, double H_si);

/* Compute the effective equation of state:
 * w_eff(a) = P_tot / (rho_tot * c^2)
 * = (Omega_r*1/3 + Omega_m*0 + Omega_Lambda*(-1) + Omega_k*(-1/3)) / sum(Omega)
 * L2 - Effective equation of state
 */
double friedmann_w_eff(const CosmoModel *model, double z);

/* === Epoch transitions === */

/* Find the redshift of matter-radiation equality:
 * Omega_m0*(1+z_eq)^3 = Omega_r0*(1+z_eq)^4
 * => 1 + z_eq = Omega_m0 / Omega_r0
 * L6 - Matter-radiation equality (canonical)
 */
double friedmann_z_eq_matter_radiation(const CosmoModel *model);

/* Find the redshift of matter-dark energy equality:
 * Omega_m0*(1+z_acc)^3 = Omega_Lambda0
 * => 1 + z_acc = (Omega_Lambda0/Omega_m0)^{1/3}
 * L6 - Cosmic acceleration onset (canonical)
 */
double friedmann_z_acceleration_onset(const CosmoModel *model);

/* === Analytic solutions (special cases) === */

/* Scale factor in a flat, matter-dominated universe (Einstein-de Sitter):
 * a(t) = (t/t_0)^{2/3}
 * H(t) = 2/(3*t)
 * t_0 = 2/(3*H0)
 * L6 - Einstein-de Sitter universe (canonical)
 */
double friedmann_scale_matter_only(double t, double H0_si);

/* Scale factor in a flat, radiation-dominated universe:
 * a(t) = (t/t_0)^{1/2}
 * H(t) = 1/(2*t)
 * L6 - Radiation-dominated expansion
 */
double friedmann_scale_radiation_only(double t, double H0_si);

/* Scale factor in a de Sitter universe (pure Lambda):
 * a(t) = a_0 * exp(H*(t - t_0))
 * H = constant = sqrt(Lambda*c^2/3)
 * L6 - de Sitter universe (canonical)
 */
double friedmann_scale_de_sitter(double t, double H_constant);

/* === Lookup table generation === */

/* Generate a redshift-array look-up table for H(z), Omega_i(z).
 * Output arrays: z_arr, H_arr, Omega_arr, t_arr each of length n.
 * z_arr is logarithmically spaced from z_min to z_max.
 * L5 - Cosmological lookup table generation
 */
void friedmann_generate_lookup_table(const CosmoModel *model,
                                      double z_min, double z_max,
                                      int n,
                                      double *z_out,
                                      double *H_out,
                                      double *t_out);

/* === ODE system for a(t) integration === */

/* Single-step RK4 integration of the Friedmann equations.
 * Given a_n at time t_n, compute a_{n+1} at t_n + dt.
 *
 * d^2a/dt^2 = a * H^2 + a * dH/dt*dH/da * a_dot
 *           = a*H^2 * (1 - 3/2*(1+w_eff))
 *
 * L5 - RK4 integration method
 */
double friedmann_rk4_step(const CosmoModel *model,
                           double a_n, double t_n, double dt);

/* Full forward integration of a(t) from a_start to a_end.
 * Uses adaptive RK4 with tolerance tol.
 * L5 - Adaptive ODE integration
 */
int friedmann_integrate_scale_factor(const CosmoModel *model,
                                      double a_start, double a_end,
                                      double dt_initial, double tol,
                                      int max_steps,
                                      double *t_out, double *a_out, int n_out);

/* Integrate cosmic time from Big Bang to redshift z.
 * t(z) = integral_z^inf [dz' / ((1+z') * H(z'))]
 * L5 - Cosmic time integration from redshift
 */
double friedmann_cosmic_time(const CosmoModel *model, double z);

#ifdef __cplusplus
}
#endif

#endif
