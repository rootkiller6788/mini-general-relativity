/**
 * @file cosmic_evolution.c
 * @brief Example: Track cosmic expansion history from early universe to today.
 *
 * Computes H(z), Omega_i(z), w_eff(z), and epoch transitions across
 * cosmic history. Visualizes the radiation-matter-Lambda transitions.
 *
 * Usage: make examples && ./build/examples/cosmic_evolution
 */

#include <stdio.h>
#include <math.h>
#include "../include/cosmo_model.h"
#include "../include/friedmann.h"
#include "../include/distances.h"

int main(void)
{
    printf("===========================================\n");
    printf("  Cosmic Expansion History\n");
    printf("===========================================\n\n");

    CosmoModel model;
    cosmo_model_init_planck(&model);

    /* Key epochs */
    double z_eq_mr = friedmann_z_eq_matter_radiation(&model);
    double z_acc = friedmann_z_acceleration_onset(&model);

    printf("Cosmology: Planck 2018 LCDM\n");
    printf("  H0 = %.1f km/s/Mpc\n", model.H0);
    printf("  Omega_m0 = %.4f\n", model.Omega_m);
    printf("  Omega_r0 = %.2e\n", model.Omega_r);
    printf("  Omega_Lambda0 = %.4f\n\n", model.Omega_Lambda);

    printf("Epoch transitions:\n");
    printf("  Matter-radiation equality: z_eq = %.1f\n", z_eq_mr);
    printf("  (t_eq ~ %.1f kyr after Big Bang)\n",
           cosmo_age_of_universe(&model, z_eq_mr) * 1000.0);
    printf("  Acceleration onset: z_acc = %.3f\n", z_acc);
    printf("  (t_acc ~ %.2f Gyr after Big Bang)\n",
           cosmo_age_of_universe(&model, z_acc));
    printf("  Age today: t_0 = %.3f Gyr\n\n",
           cosmo_age_of_universe(&model, 0.0));

    printf("Cosmic history table (logarithmic redshift sampling):\n");
    printf("%10s  %8s  %10s  %10s  %10s  %10s  %10s\n",
           "z", "t[Gyr]", "H(z)", "Omega_m", "Omega_r", "Omega_L", "w_eff");
    printf("----------  --------  ----------  ----------  ----------  ----------  ----------\n");

    double z_epochs[] = {1.0e7, 1.0e6, 1.0e5, 1.0e4, 5000, z_eq_mr,
                         1000, 100, 10, 5, z_acc, 1.0, 0.5, 0.1, 0.01, 0.0};
    int n_epoch = sizeof(z_epochs) / sizeof(z_epochs[0]);
    int i;

    for (i = 0; i < n_epoch; i++) {
        double z = z_epochs[i];
        double t = cosmo_age_of_universe(&model, z);
        double H = friedmann_H_z(&model, z);
        double w_eff = friedmann_w_eff(&model, z);

        double E_z = friedmann_E_z(&model, z);
        double Om_m_z = model.Omega_m * pow(1.0+z, 3.0) / (E_z * E_z);
        double Om_r_z = model.Omega_r * pow(1.0+z, 4.0) / (E_z * E_z);
        double Om_L_z = model.Omega_Lambda / (E_z * E_z);

        printf("%10.1f  %8.3f  %10.1f  %10.6f  %10.6f  %10.6f  %10.4f",
               z, t, H, Om_m_z, Om_r_z, Om_L_z, w_eff);

        /* Mark key transitions */
        if (fabs(z - z_eq_mr) < 1.0) printf("  <- Matter-radiation equality");
        if (fabs(z - z_acc) < 0.01) printf("  <- Acceleration onset");
        printf("\n");
    }

    printf("\nObservations:\n");
    printf("  1. Early universe (z>>%.0f): radiation dominates\n", z_eq_mr);
    printf("     w_eff ~ 1/3, a(t) ~ t^{1/2}\n\n");
    printf("  2. Matter era (%.0f > z > %.2f): matter dominates\n", z_eq_mr, z_acc);
    printf("     w_eff ~ 0, a(t) ~ t^{2/3}\n\n");
    printf("  3. Late universe (z < %.2f): dark energy dominates\n", z_acc);
    printf("     w_eff < -1/3, expansion accelerates\n\n");
    printf("  4. Future (z -> -1): pure de Sitter, w_eff -> -1\n");
    printf("     a(t) ~ exp(H*t), exponential expansion\n");
    printf("===========================================\n");

    return 0;
}
