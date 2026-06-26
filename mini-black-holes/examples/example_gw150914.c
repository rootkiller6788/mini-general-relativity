/**
 * example_gw150914.c
 *
 * Demonstrates GW150914 — the first gravitational wave detection:
 * - Binary BH parameters: chirp mass, final mass/spin
 * - Area theorem verification
 * - QNM ringdown frequency and damping time
 * - Radiated energy in gravitational waves
 *
 * L7 Applications: GW data analysis
 * L4 Fundamental Laws: Area theorem test
 *
 * Reference: Abbott et al., PRL 116, 061102 (2016)
 *
 * Build: gcc -I../include example_gw150914.c -L../build -lblackholes -lm -o example_gw150914
 */

#include "black_hole_metrics.h"
#include "black_hole_thermodynamics.h"
#include "black_hole_dynamics.h"
#include "black_hole_waves.h"
#include "black_hole_advanced.h"
#include <stdio.h>

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   GW150914 — The First Gravitational Wave Detection  ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* GW150914 parameters (LIGO-Virgo 2016) */
    double m1 = 36.0 * BH_SOLAR_MASS;
    double m2 = 29.0 * BH_SOLAR_MASS;
    double M_total = m1 + m2;

    printf("══════════════════════════════════════════════\n");
    printf("  Progenitor Parameters\n");
    printf("══════════════════════════════════════════════\n");
    printf("  m₁ = %.0f M_sun\n", m1 / BH_SOLAR_MASS);
    printf("  m₂ = %.0f M_sun\n", m2 / BH_SOLAR_MASS);
    printf("  M_total = %.0f M_sun\n", M_total / BH_SOLAR_MASS);

    double M_chirp = chirp_mass(m1, m2);
    printf("  M_chirp = %.1f M_sun\n", M_chirp / BH_SOLAR_MASS);
    printf("  (LIGO measured: 28.1^{+1.8}_{-1.5} M_sun)\n\n");

    printf("══════════════════════════════════════════════\n");
    printf("  Merger and Remnant\n");
    printf("══════════════════════════════════════════════\n");

    double chi1 = 0.3, chi2 = 0.3;
    double M_final = final_mass_after_merger(m1, m2, chi1, chi2);
    double chi_final = final_spin_after_merger(m1, m2, chi1, chi2);
    double E_rad_frac = radiated_energy_fraction(m1, m2, chi1, chi2);
    double E_rad = E_rad_frac * M_total * BH_C * BH_C;

    printf("  Final mass:     %.1f M_sun\n", M_final / BH_SOLAR_MASS);
    printf("  (LIGO measured: 62.2 ± 0.4 M_sun)\n");
    printf("  Final spin:     χ = %.3f\n", chi_final);
    printf("  (LIGO measured: 0.68 ± 0.03)\n");
    printf("  Radiated:       %.2f M_sun c²\n", (M_total - M_final) / BH_SOLAR_MASS);
    printf("  Radiated E:     %.2e J\n", E_rad);
    printf("  Radiated fraction: %.1f%%\n", E_rad_frac * 100.0);

    double L_peak = gw_luminosity_binary(m1, m2,
                     6.0 * mass_to_geometric(M_total)); /* near ISCO */
    printf("  Peak L_GW:      %.2e W\n", L_peak);
    printf("  (> all EM luminosity in observable universe!)\n\n");

    printf("══════════════════════════════════════════════\n");
    printf("  Area Theorem Test (Hawking 1971)\n");
    printf("══════════════════════════════════════════════\n");

    double M1g = mass_to_geometric(m1);
    double M2g = mass_to_geometric(m2);
    double Mfg = mass_to_geometric(M_final);

    int holds = area_theorem_merger_check(M1g, chi1*M1g, 0,
                                            M2g, chi2*M2g, 0,
                                            Mfg, chi_final*Mfg, 0);
    double A1 = kerr_horizon_area(M1g, chi1*M1g);
    double A2 = kerr_horizon_area(M2g, chi2*M2g);
    double Af = kerr_horizon_area(Mfg, chi_final*Mfg);

    printf("  A₁ = %.2e m²\n", A1);
    printf("  A₂ = %.2e m²\n", A2);
    printf("  A₁ + A₂ = %.2e m²\n", A1 + A2);
    printf("  A_final = %.2e m²\n", Af);
    printf("  ΔA = %.2e m²  %s\n", Af - (A1+A2),
           holds ? "(≥ 0 ✓ — Area Theorem holds)" :
                   "(< 0 — Area Theorem VIOLATED!)");

    printf("\n══════════════════════════════════════════════\n");
    printf("  Ringdown (Quasinormal Modes)\n");
    printf("══════════════════════════════════════════════\n");

    QuasinormalMode qnm = schwarzschild_qnm_fundamental(2, 2);
    double f_ring = qnm_frequency_hz(&qnm, M_final);
    double tau_ring = qnm_damping_time(&qnm, M_final);

    printf("  l=m=2 fundamental QNM:\n");
    printf("  Mω = %.5f %+.5f i\n", qnm.omega_R, qnm.omega_I);
    printf("  f = %.1f Hz\n", f_ring);
    printf("  τ = %.2f ms\n", tau_ring * 1000.0);
    printf("  Q = %.1f (≈ %d visible cycles)\n",
           qnm_quality_factor(&qnm), (int)qnm_quality_factor(&qnm));
    printf("  (LIGO observed: f ≈ 251 Hz, τ ≈ 4.0 ms)\n\n");

    printf("══════════════════════════════════════════════\n");
    printf("  No-Hair Theorem Check\n");
    printf("══════════════════════════════════════════════\n");

    KerrMultipoles mp;
    kerr_multipoles(Mfg, chi_final * Mfg, &mp);
    double nh_ratio = no_hair_test_quadrupole(mp.M2, Mfg, mp.J);
    printf("  M₂ / (J²/M) = %.3f  (Kerr prediction: -1.0)\n", nh_ratio);
    printf("  (LIGO-Virgo 2019: -1.0 ± 0.3)\n\n");

    return 0;
}
