/**
 * Example: Light Deflection by the Sun (Eddington 1919)
 *
 * Computes the GR prediction for starlight deflection during
 * a solar eclipse and compares with the Newtonian prediction.
 * The factor-of-2 difference was the first experimental
 * confirmation of General Relativity.
 *
 * Reference: Dyson, Eddington & Davidson (1920)
 * "A Determination of the Deflection of Light by the 
 *  Sun's Gravitational Field from Observations made at 
 *  the Total Eclipse of May 29, 1919"
 */

#include <stdio.h>
#include <math.h>
#include "schwarzschild_observables.h"
#include "schwarzschild_defs.h"

int main(void) {
    printf("============================================\n");
    printf("Light Deflection by the Sun (Eddington 1919)\n");
    printf("============================================\n\n");

    double M_sun = SCHW_M_SUN;
    double R_sun = 6.957e8;     /* Solar radius [m] */
    double c     = SCHW_C;
    double G     = SCHW_G_N;

    /* GR deflection: alpha = 4*G*M/(b*c^2) */
    double alpha_GR = schwarzschild_light_deflection_angle(R_sun, M_sun);
    double alpha_GR_arcsec = alpha_GR * 180.0/M_PI * 3600.0;

    /* Newtonian deflection (half of GR) */
    double alpha_Newton = 2.0 * G * M_sun / (R_sun * c * c);
    double alpha_Newton_arcsec = alpha_Newton * 180.0/M_PI * 3600.0;

    printf("Impact parameter (solar limb): %.3e m\n", R_sun);
    printf("\nPredicted deflection angles:\n");
    printf("  Newtonian:  %.4e rad  =  %.4f arcsec\n", alpha_Newton, alpha_Newton_arcsec);
    printf("  Einstein:   %.4e rad  =  %.4f arcsec\n", alpha_GR, alpha_GR_arcsec);
    printf("  Ratio GR/Newton: %.4f\n", alpha_GR/alpha_Newton);
    printf("\nObserved (Eddington 1919):\n");
    printf("  Sobral (Brazil):    1.98 +/- 0.16 arcsec\n");
    printf("  Principe (Africa):  1.61 +/- 0.40 arcsec\n");
    printf("  GR prediction:      1.75 arcsec\n\n");

    /* Deflection at different impact parameters */
    printf("Deflection vs Impact Parameter:\n");
    printf("  b/R_sun    deflection [arcsec]\n");
    printf("  -------------------------------\n");
    for (int i = 1; i <= 10; i++) {
        double b = i * R_sun;
        double alpha = schwarzschild_light_deflection_angle(b, M_sun);
        printf("  %5.1f      %8.4f\n", (double)i, alpha*180.0/M_PI*3600.0);
    }

    /* Solar limb deflection for different stars */
    printf("\nNote: The 1.75 arcsec deflection at the solar limb\n");
    printf("decreases as 1/b with increasing impact parameter.\n");
    printf("Observations require a total solar eclipse because\n");
    printf("the Sun's brightness overwhelms starlight.\n");
    printf("The 1919 Eddington expedition confirmed GR and made\n");
    printf("Einstein world-famous.\n");

    return 0;
}
