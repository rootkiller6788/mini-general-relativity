/**
 * Example: Schwarzschild Black Hole Properties
 *
 * Computes key properties of Schwarzschild black holes
 * across a range of astrophysical masses, from stellar-mass
 * to supermassive, including the M87* black hole imaged
 * by the Event Horizon Telescope.
 *
 * Reference: Wald (1984) Ch.6, Ch.12; EHT Collaboration (2019)
 */

#include <stdio.h>
#include <math.h>
#include "schwarzschild_defs.h"
#include "schwarzschild_geodesics.h"
#include "schwarzschild_observables.h"
#include "schwarzschild_extensions.h"

int main(void) {
    printf("=====================================================\n");
    printf("Schwarzschild Black Hole Properties\n");
    printf("=====================================================\n\n");

    double M_sun = SCHW_M_SUN;
    double masses[] = {
        3.0 * M_sun,        /* Stellar BH */
        10.0 * M_sun,       /* Typical X-ray binary BH */
        30.0 * M_sun,       /* GW150914-like */
        4.1e6 * M_sun,      /* Sagittarius A* (Milky Way center) */
        6.5e9 * M_sun       /* M87* (EHT 2019) */
    };
    const char *labels[] = {
        "3 Msun (stellar)", "10 Msun (XRB)", "30 Msun (GW150914)",
        "4.1e6 Msun (Sgr A*)", "6.5e9 Msun (M87*)"
    };

    printf("%-25s %10s %10s %10s %10s %10s\n",
           "Mass", "rs [km]", "rho [kg/m3]", "T_H [K]", "t_evap [yr]", "S [J/K]");
    printf("-----------------------------------------------------------"
           "------------------------------------\n");

    for (int i = 0; i < 5; i++) {
        double M = masses[i];
        double rs = schwarzschild_radius(M);
        double rs_km = rs / 1e3;
        double V = (4.0/3.0) * M_PI * rs * rs * rs;
        double rho = M / V;
        double T_H = schwarzschild_hawking_temperature(M);
        double t_evap = schwarzschild_evaporation_time(M) / (365.25 * 86400.0);
        double S = schwarzschild_bekenstein_hawking_entropy(M);

        printf("%-25s %10.3e %10.3e %10.3e %10.3e %10.3e\n",
               labels[i], rs_km, rho, T_H, t_evap, S);
    }

    printf("\nKey Radii (in units of rs = 2GM/c^2):\n");
    double rs_ref = 1.0e4;  /* reference Schwarzschild radius */
    printf("  Photon sphere:     r_ph  = %.1f rs\n",
           schwarzschild_photon_sphere_radius(rs_ref) / rs_ref);
    printf("  ISCO:              r_isco = %.1f rs\n",
           schwarzschild_isco_radius(rs_ref) / rs_ref);
    printf("  Marginally bound:  r_mb  = %.1f rs\n",
           schwarzschild_marginally_bound_radius(rs_ref) / rs_ref);
    printf("  Photon capture:    b_crit = %.3f rs\n",
           schwarzschild_critical_impact_parameter(rs_ref) / rs_ref);

    printf("\nAccretion Efficiency:\n");
    printf("  Schwarzschild:  eta = 1 - sqrt(8/9) = %.1f%%\n",
           100.0 * schwarzschild_accretion_efficiency());

    printf("\nShadow Angular Radius:\n");
    double M_M87 = 6.5e9 * M_sun;
    double D_M87 = 16.8e6 * SCHW_PC;
    double theta_shadow = schwarzschild_shadow_angular_radius(M_M87, D_M87);
    double theta_uas = theta_shadow * 180.0/M_PI * 3600.0 * 1e6;
    printf("  M87* at 16.8 Mpc:  theta_shadow = %.1f uas\n", theta_uas);
    printf("  (EHT measured ~42 uas, consistent with 6.5e9 Msun)\n");

    printf("\nPhysical Insights:\n");
    printf("  - Stellar BHs have T_H << CMB (2.7 K), net mass gain from CMB\n");
    printf("  - Supermassive BHs have very low density (< water for M87*!)\n");
    printf("  - t_evap >> age of universe for all astrophysical BHs\n");
    printf("  - Only primordial BHs with M < 1e12 kg could have evaporated\n");
    return 0;
}
