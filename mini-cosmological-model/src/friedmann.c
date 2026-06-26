/** @file friedmann.c
 *  @brief Friedmann equations: Hubble parameter, scale factor evolution,
 *         density evolution, ODE integration, and epoch transitions.
 */

#include "friedmann.h"
#include <math.h>
#include <string.h>

/* Density evolution for constant-w fluid: rho(a) = rho0 * a^{-3(1+w)} */
double friedmann_rho_a(double rho0, double a, double w)
{
    if (a <= 0.0) return 0.0;
    double exponent = -3.0 * (1.0 + w);
    return rho0 * pow(a, exponent);
}

/* Omega_i(a) = Omega_i0 * a^{-3(1+w_i)} * (H0/H(a))^2 */
double friedmann_omega_a(double Omega_i0, double a, double w,
                          double H0_si, double H_si)
{
    if (a <= 0.0 || H_si <= 0.0) return 0.0;
    double rho_ratio = pow(a, -3.0 * (1.0 + w));
    double H_ratio = (H0_si * H0_si) / (H_si * H_si);
    return Omega_i0 * rho_ratio * H_ratio;
}

/* E(z) = H(z)/H0 = sqrt(Omega_r0*(1+z)^4 + Omega_m0*(1+z)^3
 *                     + Omega_k0*(1+z)^2 + Omega_Lambda0) */
double friedmann_E_z(const CosmoModel *model, double z)
{
    if (!model) return 0.0;
    double zp1  = 1.0 + z;
    double zp1_2 = zp1 * zp1;
    double zp1_3 = zp1_2 * zp1;
    double zp1_4 = zp1_3 * zp1;
    double E2 = model->Omega_r      * zp1_4
              + model->Omega_m      * zp1_3
              + model->Omega_k      * zp1_2
              + model->Omega_Lambda;
    if (E2 < 0.0) E2 = 0.0;
    return sqrt(E2);
}

double friedmann_H_z(const CosmoModel *model, double z)
{
    if (!model) return 0.0;
    return model->H0 * friedmann_E_z(model, z);
}

double friedmann_H_z_si(const CosmoModel *model, double z)
{
    if (!model) return 0.0;
    return friedmann_H_z(model, z) * 1.0e3 / COSMO_MPC;
}

double friedmann_dadt(double a, double H_a_si)
{
    return a * H_a_si;
}


/* dH/da from second Friedmann equation */
double friedmann_dH_da(const CosmoModel *model, double a, double H_si)
{
    if (!model || a <= 0.0 || H_si <= 0.0) return 0.0;
    double z = 1.0 / a - 1.0;
    double w = friedmann_w_eff(model, z);
    return -(3.0 * H_si) / (2.0 * a) * (1.0 + w);
}

/* w_eff(z) = P_tot/(rho_tot*c^2) */
double friedmann_w_eff(const CosmoModel *model, double z)
{
    if (!model) return 0.0;
    double zp1  = 1.0 + z;
    double zp1_2 = zp1 * zp1;
    double zp1_3 = zp1_2 * zp1;
    double zp1_4 = zp1_3 * zp1;
    double E2 = model->Omega_r      * zp1_4
              + model->Omega_m      * zp1_3
              + model->Omega_k      * zp1_2
              + model->Omega_Lambda;
    if (E2 <= 0.0) return 0.0;
    double P_over_rho_c2 = model->Omega_r * zp1_4 * (1.0/3.0)
                         + model->Omega_m * zp1_3 *  0.0
                         + model->Omega_k * zp1_2 * (-1.0/3.0)
                         + model->Omega_Lambda * (-1.0);
    return P_over_rho_c2 / E2;
}

/* Matter-radiation equality redshift */
double friedmann_z_eq_matter_radiation(const CosmoModel *model)
{
    if (!model || model->Omega_r <= 0.0) return -1.0;
    return model->Omega_m / model->Omega_r - 1.0;
}

/* Acceleration onset redshift */
double friedmann_z_acceleration_onset(const CosmoModel *model)
{
    if (!model || model->Omega_m <= 0.0) return -1.0;
    return pow(2.0 * model->Omega_Lambda / model->Omega_m, 1.0/3.0) - 1.0;
}

/* Einstein-de Sitter: a(t) = (t/t_0)^{2/3} */
double friedmann_scale_matter_only(double t, double H0_si)
{
    if (t <= 0.0 || H0_si <= 0.0) return 0.0;
    double t0 = 2.0 / (3.0 * H0_si);
    return pow(t / t0, 2.0/3.0);
}

/* Radiation-dominated: a(t) = sqrt(2*H0*t) */
double friedmann_scale_radiation_only(double t, double H0_si)
{
    if (t <= 0.0 || H0_si <= 0.0) return 0.0;
    return sqrt(2.0 * H0_si * t);
}

/* de Sitter: a(t) = exp(H*t) */
double friedmann_scale_de_sitter(double t, double H_constant)
{
    if (H_constant <= 0.0) return 1.0;
    return exp(H_constant * t);
}


/* Generate lookup table for H(z) and t(z) */
void friedmann_generate_lookup_table(const CosmoModel *model,
                                      double z_min, double z_max,
                                      int n,
                                      double *z_out,
                                      double *H_out,
                                      double *t_out)
{
    if (!model || !z_out || !H_out || !t_out || n < 2) return;
    double log_zmin = log(z_min + 1.0);
    double log_zmax = log(z_max + 1.0);
    double d_logz = (log_zmax - log_zmin) / (n - 1);
    int i;
    for (i = 0; i < n; i++) {
        double log_zp1 = log_zmin + i * d_logz;
        double z = exp(log_zp1) - 1.0;
        z_out[i] = z;
        H_out[i] = friedmann_H_z(model, z);
        t_out[i] = friedmann_cosmic_time(model, z);
    }
}

/* Single RK4 step for da/dt = a*H(a) */
double friedmann_rk4_step(const CosmoModel *model,
                           double a_n, double t_n, double dt)
{
    if (!model || a_n <= 0.0) return a_n;
    (void)t_n;
    double z_n = 1.0 / a_n - 1.0;
    double H_n = friedmann_H_z_si(model, z_n);
    double k1 = friedmann_dadt(a_n, H_n);

    double a2 = a_n + 0.5 * k1 * dt;
    double z2 = 1.0 / a2 - 1.0;
    double H2 = friedmann_H_z_si(model, z2);
    double k2 = friedmann_dadt(a2, H2);

    double a3 = a_n + 0.5 * k2 * dt;
    double z3 = 1.0 / a3 - 1.0;
    double H3 = friedmann_H_z_si(model, z3);
    double k3 = friedmann_dadt(a3, H3);

    double a4 = a_n + k3 * dt;
    double z4 = 1.0 / a4 - 1.0;
    double H4 = friedmann_H_z_si(model, z4);
    double k4 = friedmann_dadt(a4, H4);

    return a_n + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

/* Adaptive RK4 integration of a(t) with step-size control */
int friedmann_integrate_scale_factor(const CosmoModel *model,
                                      double a_start, double a_end,
                                      double dt_initial, double tol,
                                      int max_steps,
                                      double *t_out, double *a_out, int n_out)
{
    if (!model || !t_out || !a_out || n_out < 2) return -1;
    if (a_start <= 0.0 || a_end <= a_start) return -1;

    double a = a_start;
    double t = 0.0;
    double dt = dt_initial;
    int step = 0;
    int out_idx = 0;

    t_out[out_idx] = t;
    a_out[out_idx] = a;
    out_idx++;

    while (a < a_end && step < max_steps && out_idx < n_out) {
        double a_full = friedmann_rk4_step(model, a, t, dt);
        double a_half = friedmann_rk4_step(model, a, t, 0.5 * dt);
        double a_half2 = friedmann_rk4_step(model, a_half, t + 0.5*dt, 0.5*dt);

        double error = fabs(a_half2 - a_full);
        double rel_error = (fabs(a_half2) > 1e-30) ? error / fabs(a_half2) : error;

        if (rel_error < tol) {
            a = a_half2;
            t += dt;
            t_out[out_idx] = t;
            a_out[out_idx] = a;
            out_idx++;
            if (rel_error < tol * 0.1) dt *= 2.0;
        } else {
            dt *= 0.5;
        }

        if (dt > 1.0e16) dt = 1.0e16;
        if (dt < 1.0) dt = 1.0;
        step++;
    }

    return out_idx;
}

/* Cosmic time at redshift z: t(z) = integral_z^inf dz'/((1+z')*H(z')) */
double friedmann_cosmic_time(const CosmoModel *model, double z)
{
    if (!model) return -1.0;

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

        double H_i   = friedmann_H_z_si(model, z_i);
        double H_mid = friedmann_H_z_si(model, z_mid);
        double H_next = friedmann_H_z_si(model, z_next);

        double f_i   = (H_i   > 0.0) ? 1.0 / (zp1_i   * H_i)   : 0.0;
        double f_mid = (H_mid > 0.0) ? 1.0 / (zp1_mid * H_mid) : 0.0;
        double f_next = (H_next > 0.0) ? 1.0 / (zp1_next * H_next) : 0.0;

        double dz = z_next - z_i;
        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return integral;
}
