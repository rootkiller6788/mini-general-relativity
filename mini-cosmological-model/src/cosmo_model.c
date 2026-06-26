/** @file cosmo_model.c
 *  @brief Implementation of core cosmological model initialization and utilities.
 *
 *  Reference: Planck Collaboration (2020), A&A 641, A6
 *             Dodelson & Schmidt, Modern Cosmology (2020)
 */

#include "cosmo_model.h"
#include "friedmann.h"
#include <math.h>
#include <string.h>

/* Initialize with Planck 2018 baseline cosmology parameters. */
void cosmo_model_init_planck(CosmoModel *model)
{
    if (!model) return;

    memset(model, 0, sizeof(CosmoModel));

    /* Hubble constant */
    model->H0 = COSMO_H0_PLANCK;
    model->h  = COSMO_H0_PLANCK / 100.0;

    /* Density parameters */
    model->Omega_m      = COSMO_OMEGA_M0;
    model->Omega_b      = COSMO_OMEGA_B0;
    model->Omega_cdm    = COSMO_OMEGA_CDM0;
    model->Omega_r      = COSMO_OMEGA_R0;
    /* Enforce flatness: Omega_Lambda = 1 - Omega_m - Omega_r (neglect Omega_k) */
    model->Omega_Lambda = 1.0 - model->Omega_m - model->Omega_r;
    model->Omega_k = 0.0;

    /* Primordial power spectrum */
    model->A_s      = COSMO_A_S;
    model->n_s      = COSMO_N_S;
    model->k_pivot  = 0.05; /* Mpc^{-1} */

    /* Reionization and neutrinos */
    model->tau_reio  = COSMO_TAU_REIO;
    model->N_eff     = 3.046;
    model->sum_m_nu  = 0.06; /* eV, normal hierarchy minimum */

    /* Derived quantities */
    model->H0_si    = model->H0 * 1.0e3 / COSMO_MPC;
    model->T_cmb0   = COSMO_T_CMB;
    model->rho_crit0 = 3.0 * model->H0_si * model->H0_si / (8.0 * M_PI * COSMO_G_NEWTON);

    /* Integration settings */
    model->n_z_steps = 1000;
}

/* Initialize with user-specified parameters and compute derived quantities. */
void cosmo_model_init(CosmoModel *model, double H0, double Omega_m,
                      double Omega_r, double Omega_Lambda, double w0)
{
    if (!model) return;

    memset(model, 0, sizeof(CosmoModel));

    model->H0          = H0;
    model->h           = H0 / 100.0;
    model->Omega_m     = Omega_m;
    model->Omega_r     = Omega_r;
    model->Omega_Lambda = Omega_Lambda;

    /* Compute curvature from sum rule */
    model->Omega_k = 1.0 - (Omega_m + Omega_r + Omega_Lambda);

    /* Partition matter into baryons and CDM using Planck ratio */
    double f_b = COSMO_OMEGA_B0 / COSMO_OMEGA_M0;
    model->Omega_b   = Omega_m * f_b;
    model->Omega_cdm = Omega_m * (1.0 - f_b);

    /* Default primordial spectrum (Planck 2018) */
    model->A_s     = COSMO_A_S;
    model->n_s     = COSMO_N_S;
    model->k_pivot = 0.05;

    /* Default neutrino and reionization */
    model->N_eff    = 3.046;
    model->sum_m_nu = 0.06;
    model->tau_reio = COSMO_TAU_REIO;

    /* Derived quantities */
    model->H0_si    = H0 * 1.0e3 / COSMO_MPC;
    model->T_cmb0   = COSMO_T_CMB;
    model->rho_crit0 = 3.0 * model->H0_si * model->H0_si / (8.0 * M_PI * COSMO_G_NEWTON);

    model->n_z_steps = 1000;

    /* We store w0 but do not yet implement general w(z) models */
    (void)w0;
}

/* Validate physical consistency of the cosmological model. */
int cosmo_model_validate(const CosmoModel *model)
{
    if (!model) return 0;

    /* Check H0 is positive */
    if (model->H0 <= 0.0) return 0;

    /* Check CMB temperature is positive */
    if (model->T_cmb0 <= 0.0) return 0;

    /* Check density parameters are in physically reasonable range */
    if (model->Omega_m < -1e-6 || model->Omega_m > 2.0 + 1e-6) return 0;
    if (model->Omega_r < -1e-6 || model->Omega_r > 2.0 + 1e-6) return 0;
    if (model->Omega_Lambda < -2.0 || model->Omega_Lambda > 5.0) return 0;

    /* Sum rule: Omega_total = 1 */
    double sum = model->Omega_m + model->Omega_r + model->Omega_Lambda;
    if (fabs(sum - 1.0) > 0.1) return 0;

    /* Derived quantities must be positive */
    if (model->H0_si <= 0.0 || model->rho_crit0 <= 0.0) return 0;

    /* Primordial spectrum sanity checks */
    if (model->A_s <= 0.0 || model->n_s < 0.5 || model->n_s > 1.5) return 0;
    if (model->k_pivot <= 0.0) return 0;

    return 1;
}

/* Compute critical density at given H. */
double cosmo_rho_critical(double H_si)
{
    return 3.0 * H_si * H_si / (8.0 * M_PI * COSMO_G_NEWTON);
}

/* Integrate the age of the universe using dt/dz = 1/[H(z)*(1+z)].
 * We integrate from z to infinity (a=0) for cosmic time at redshift z. */
double cosmo_age_of_universe(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;

    /* Use Simpson integration with n_z_steps points in log(1+z). */
    int n = model->n_z_steps;
    if (n < 10) n = 1000;

    double z_max = 1.0e6;  /* effectively infinite redshift (radiation era) */
    double log_z_min = log(1.0 + z);
    double log_z_max = log(1.0 + z_max);
    double d_logz = (log_z_max - log_z_min) / n;

    double integral = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double log_z_i = log_z_min + i * d_logz;
        double log_z_next = log_z_i + d_logz;
        double log_z_mid = 0.5 * (log_z_i + log_z_next);

        double z_i   = exp(log_z_i)   - 1.0;
        double z_mid = exp(log_z_mid) - 1.0;
        double z_next = exp(log_z_next) - 1.0;

        /* dt/dz factor at each point */
        double H_i   = friedmann_H_z_si(model, z_i);
        double H_mid = friedmann_H_z_si(model, z_mid);
        double H_next = friedmann_H_z_si(model, z_next);

        double f_i   = (H_i   > 0.0) ? 1.0 / ((1.0 + z_i)   * H_i)   : 0.0;
        double f_mid = (H_mid > 0.0) ? 1.0 / ((1.0 + z_mid) * H_mid) : 0.0;
        double f_next = (H_next > 0.0) ? 1.0 / ((1.0 + z_next) * H_next) : 0.0;

        /* Simpson's rule on the sub-interval */
        double dz = z_next - z_i;
        integral += (dz / 6.0) * (f_i + 4.0 * f_mid + f_next);
    }

    return integral / COSMO_GYR_SEC; /* Convert to Gyr */
}

/* Compute lookback time from z=0 to z. */
double cosmo_lookback_time(const CosmoModel *model, double z)
{
    if (!model || z < 0.0) return -1.0;

    double age0 = cosmo_age_of_universe(model, 0.0);
    double age_z = cosmo_age_of_universe(model, z);

    if (age0 < 0.0 || age_z < 0.0) return -1.0;
    return age0 - age_z;
}
