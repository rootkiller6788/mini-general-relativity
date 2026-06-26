#ifndef COSMO_MODEL_H
#define COSMO_MODEL_H

#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define COSMO_C_LIGHT         2.99792458e8
#define COSMO_G_NEWTON        6.67430e-11
#define COSMO_M_PLANCK_RED    2.435e18
#define COSMO_M_SOLAR         1.98847e30
#define COSMO_MPC             3.085677581e22
#define COSMO_YEAR_SEC        3.15576e7
#define COSMO_GYR_SEC         3.15576e16
#define COSMO_SIGMA_SB        5.670374419e-8
#define COSMO_K_BOLTZMANN     1.380649e-23
#define COSMO_SIGMA_THOMSON   6.6524587321e-29
#define COSMO_M_PROTON        1.67262192369e-27
#define COSMO_M_ELECTRON      9.1093837015e-31

#define COSMO_H0_PLANCK       67.4
#define COSMO_H0_SI           (COSMO_H0_PLANCK * 1.0e3 / COSMO_MPC)
#define COSMO_T_CMB           2.7255
#define COSMO_OMEGA_M0        0.315
#define COSMO_OMEGA_B0        0.0493
#define COSMO_OMEGA_CDM0      0.2657
#define COSMO_OMEGA_LAMBDA0   0.685
#define COSMO_OMEGA_R0        9.146e-5
#define COSMO_OMEGA_K0        0.0
#define COSMO_A_S             2.100e-9
#define COSMO_N_S             0.9649
#define COSMO_TAU_REIO        0.054
#define COSMO_R_DRAG          147.09

typedef enum {
    COSMO_SPECIES_RADIATION = 0,
    COSMO_SPECIES_MATTER    = 1,
    COSMO_SPECIES_CURVATURE = 2,
    COSMO_SPECIES_LAMBDA    = 3,
    COSMO_SPECIES_DARK_ENERGY = 4,
    COSMO_SPECIES_COUNT     = 5
} CosmoSpecies;

typedef enum {
    COSMO_K_HYPERBOLIC  = -1,
    COSMO_K_FLAT        =  0,
    COSMO_K_SPHERICAL   = +1
} CosmoCurvature;

typedef enum {
    COSMO_EPOCH_INFLATION   = 0,
    COSMO_EPOCH_RADIATION   = 1,
    COSMO_EPOCH_MATTER      = 2,
    COSMO_EPOCH_LAMBDA      = 3,
    COSMO_EPOCH_COUNT       = 4
} CosmoEpoch;

typedef struct {
    double H0;
    double h;
    double Omega_m;
    double Omega_b;
    double Omega_cdm;
    double Omega_r;
    double Omega_Lambda;
    double Omega_k;
    double A_s;
    double n_s;
    double k_pivot;
    double tau_reio;
    double N_eff;
    double sum_m_nu;
    double H0_si;
    double rho_crit0;
    double T_cmb0;
    int    n_z_steps;
} CosmoModel;

typedef struct {
    double z;
    double a;
    double t;
    double H;
    double H_si;
    double rho_crit;
    double rho[5];
    double Omega[5];
    double w_eff;
    double q;
    double j;
} CosmoState;

void cosmo_model_init_planck(CosmoModel *model);
void cosmo_model_init(CosmoModel *model, double H0, double Omega_m,
                      double Omega_r, double Omega_Lambda, double w0);
int cosmo_model_validate(const CosmoModel *model);
double cosmo_rho_critical(double H_si);
double cosmo_age_of_universe(const CosmoModel *model, double z);
double cosmo_lookback_time(const CosmoModel *model, double z);

#ifdef __cplusplus
}
#endif

#endif
