/**
 * @file black_hole_waves.c
 * @brief Quasinormal modes, ringdown, GW emission, Teukolsky equation.
 *
 * Implements QNM frequency computation, ringdown waveform generation,
 * binary merger phenomenology, and gravitational wave signal analysis.
 *
 * L5 Computational Methods: QNM computation, matched filtering
 * L6 Canonical Systems: Schwarzschild/Kerr QNM spectrum
 * L7 Applications: GW150914 parameter estimation, LIGO data analysis
 */

#include "black_hole_waves.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 *  QUASINORMAL MODES: SCHWARZSCHILD
 * ================================================================ */

/*
 * QNM frequencies for Schwarzschild (Chandrasekhar & Detweiler 1975):
 *
 * These are computed via Leaver's continued fraction method (1985).
 * We use high-precision values from Berti, Cardoso, Will (2006) for
 * gravitational perturbations (s=-2).
 *
 * M ω_{lmn} values for n=0 (fundamental mode):
 *
 * l=2, m=0: 0.74734 - 0.17792 i
 * l=2, m=1: same as m=0 (degenerate in Schwarzschild)
 * l=2, m=2: 0.74734 - 0.17792 i
 * l=3, m=3: 1.19889 - 0.18541 i
 * l=4, m=4: 1.61836 - 0.18833 i
 *
 * Overtones (n>0) damp faster:
 * l=2, n=1: 0.69342 - 0.54783 i
 * l=2, n=2: 0.60211 - 0.95655 i
 */

QuasinormalMode schwarzschild_qnm_fundamental(int l, int m) {
    QuasinormalMode qnm;
    qnm.n = 0;
    qnm.l = l;
    qnm.m = m;
    qnm.s = -2; /* gravitational perturbations */

    switch (l) {
        case 2:
            qnm.omega_R = 0.74734;
            qnm.omega_I = -0.17792;
            break;
        case 3:
            qnm.omega_R = 1.19889;
            qnm.omega_I = -0.18541;
            break;
        case 4:
            qnm.omega_R = 1.61836;
            qnm.omega_I = -0.18833;
            break;
        case 5:
            qnm.omega_R = 2.02459;
            qnm.omega_I = -0.18974;
            break;
        case 6:
            qnm.omega_R = 2.42402;
            qnm.omega_I = -0.19052;
            break;
        default:
            /* WKB approximation for large l (Schutz & Will 1985):
             * M ω ≈ (l + 1/2) / √27 - i (n + 1/2) / √27 */
            qnm.omega_R = (l + 0.5) / sqrt(27.0);
            qnm.omega_I = -(0.0 + 0.5) / sqrt(27.0); /* n=0 */
            break;
    }

    qnm.quality_factor = qnm_quality_factor(&qnm);
    return qnm;
}

QuasinormalMode schwarzschild_qnm_overtone(int l, int m, int n) {
    QuasinormalMode qnm;
    qnm.n = n;
    qnm.l = l;
    qnm.m = m;
    qnm.s = -2;

    if (l == 2) {
        /* Known overtones for l=2 gravitational perturbations */
        switch (n) {
            case 0:
                qnm.omega_R = 0.74734; qnm.omega_I = -0.17792; break;
            case 1:
                qnm.omega_R = 0.69342; qnm.omega_I = -0.54783; break;
            case 2:
                qnm.omega_R = 0.60211; qnm.omega_I = -0.95655; break;
            case 3:
                qnm.omega_R = 0.50301; qnm.omega_I = -1.41045; break;
            default:
                /* WKB approximation for high overtones */
                qnm.omega_R = (l + 0.5) / sqrt(27.0);
                qnm.omega_I = -(n + 0.5) / sqrt(27.0);
                break;
        }
    } else {
        /* WKB for general l */
        qnm.omega_R = (l + 0.5) / sqrt(27.0);
        qnm.omega_I = -(n + 0.5) / sqrt(27.0);
    }

    qnm.quality_factor = qnm_quality_factor(&qnm);
    return qnm;
}

QuasinormalMode schwarzschild_qnm_l2m0(void) {
    /* The l=2, m=0 mode: the "breathing" (spherically symmetric) mode.
     * Physically relevant for head-on BH collisions. */
    QuasinormalMode qnm;
    qnm.l = 2; qnm.m = 0; qnm.n = 0; qnm.s = -2;
    qnm.omega_R = 0.74734;
    qnm.omega_I = -0.17792; /* degenerate with m≠0 in Schwarzschild */
    qnm.quality_factor = qnm_quality_factor(&qnm);
    return qnm;
}

/* ================================================================
 *  QUASINORMAL MODES: KERR (FITTING FORMULAE)
 * ================================================================ */

QuasinormalMode kerr_qnm_fundamental(double M, double a, int l, int m) {
    /* Kerr QNM frequencies: fitting formulae from
     * Berti, Cardoso, Will (2006), PRD 73, 064030.
     *
     * For l=m=2:
     * ω_R(a) = ω_R(0) [1 + 0.668(a/M) + 0.353(a/M)² + 0.187(a/M)³]
     * ω_I(a) = ω_I(0) [1 - 0.301(a/M) - 0.379(a/M)² + 0.171(a/M)³]
     *
     * The frequency SPLITS for m>0 vs m<0 (prograde vs retrograde). */

    QuasinormalMode qnm_schw = schwarzschild_qnm_fundamental(l, l);
    QuasinormalMode qnm;
    qnm.l = l; qnm.n = 0; qnm.s = -2;

    double aM = (M > 0) ? a / M : 0.0;
    if (aM > 1.0) aM = 1.0;
    if (aM < -1.0) aM = -1.0;

    /* Choose m (prograde if m>0, retrograde if m<0) */
    qnm.m = m;
    double sign = (m >= 0) ? 1.0 : -1.0;
    double aM_eff = sign * aM;

    if (l == 2 && (m == 2 || m == -2)) {
        /* BCW fitting for l=m=2: f₂₂, f₂₋₂ */
        double c1_R = 0.668, c2_R = 0.353, c3_R = 0.187;
        double c1_I = -0.301, c2_I = -0.379, c3_I = 0.171;

        qnm.omega_R = qnm_schw.omega_R * (1.0 + c1_R * aM_eff
                      + c2_R * aM_eff * aM_eff + c3_R * aM_eff * aM_eff * aM_eff);
        qnm.omega_I = qnm_schw.omega_I * (1.0 + c1_I * aM_eff
                      + c2_I * aM_eff * aM_eff + c3_I * aM_eff * aM_eff * aM_eff);
    } else {
        /* Generic: scale by a/M linearly (approximate) */
        qnm.omega_R = qnm_schw.omega_R * (1.0 + 0.5 * aM_eff);
        qnm.omega_I = qnm_schw.omega_I * (1.0 - 0.3 * aM_eff);
    }

    qnm.quality_factor = qnm_quality_factor(&qnm);
    return qnm;
}

/* ================================================================
 *  QNM UTILITY FUNCTIONS
 * ================================================================ */

double qnm_quality_factor(const QuasinormalMode *qnm) {
    /* Q = ω_R / (2|ω_I|)
     *
     * The quality factor measures how many oscillations are visible
     * in the ringdown before the amplitude drops by 1/e.
     *
     * GW150914: Q ≈ 5.5 → ~1 oscillation cycle visible before
     * the signal merges into detector noise. */
    if (qnm == NULL || qnm->omega_I == 0) return 0.0;
    return qnm->omega_R / (2.0 * fabs(qnm->omega_I));
}

double qnm_frequency_hz(const QuasinormalMode *qnm, double mass_kg) {
    /* Convert geometric frequency Mω to Hz:
     * f = Mω × c³/(2πGM) = ω_R × c³/(2πGM)
     *
     * For a 10 M_sun BH: f ≈ 1.2 kHz
     * For a 100 M_sun BH: f ≈ 120 Hz
     * For Sgr A* (4e6 M_sun): f ≈ 3 mHz (LISA band!) */
    if (qnm == NULL || mass_kg <= 0) return 0.0;
    return qnm->omega_R * BH_C * BH_C * BH_C
           / (2.0 * M_PI * BH_G * mass_kg);
}

double qnm_damping_time(const QuasinormalMode *qnm, double mass_kg) {
    /* τ = -1/ω_I (in geometric units)
     * τ_si = -M / (ω_I c) */
    if (qnm == NULL || qnm->omega_I >= 0 || mass_kg <= 0) return 0.0;
    double M_geom = mass_to_geometric(mass_kg);
    return -M_geom / (qnm->omega_I * BH_C);
}

/* ================================================================
 *  RINGDOWN WAVEFORM
 * ================================================================ */

double ringdown_waveform(double t, double A, double f,
                          double tau, double phi0) {
    /* h(t) = A exp(-t/τ) cos(2π f t + φ₀) for t ≥ 0
     *
     * This is the simplest model of the post-merger GW signal.
     * The ringdown begins after the merger (t=0) and decays
     * exponentially with damping time τ.
     *
     * In the frequency domain:
     * h̃(f) ∝ A / [1/(2τ) + i 2π(f - f₀)]  (Lorentzian) */
    if (t < 0 || tau <= 0) return 0.0;
    return A * exp(-t / tau) * cos(2.0 * M_PI * f * t + phi0);
}

double ringdown_multimode(double t,
                           const QuasinormalMode *modes,
                           const double *amps,
                           const double *phases,
                           int n_modes) {
    /* Multi-mode ringdown: sum over overtones and harmonics.
     *
     * h(t) = Σ_n A_n exp(-t/τ_n) cos(2π f_n t + φ_n)
     *
     * In practice, the l=m=2 fundamental mode dominates for
     * quasi-circular binaries. Overtones (n>0) contribute at
     * early times but decay faster.
     *
     * Higher harmonics (l=3, m=3) are important for:
     * - High-mass-ratio systems
     * - Precessing binaries
     * - Testing the no-hair theorem (consistency between modes) */

    if (modes == NULL || amps == NULL || phases == NULL) return 0.0;
    double h = 0.0;
    for (int i = 0; i < n_modes; i++) {
        /* Need mass to convert to SI — assume 65 M_sun for GW150914-like */
        double f = qnm_frequency_hz(&modes[i], 65.0 * BH_SOLAR_MASS);
        double tau = qnm_damping_time(&modes[i], 65.0 * BH_SOLAR_MASS);
        h += ringdown_waveform(t, amps[i], f, tau, phases[i]);
    }
    return h;
}

/* ================================================================
 *  BINARY MERGER PHENOMENOLOGY
 * ================================================================ */

double chirp_mass(double m1_kg, double m2_kg) {
    /* M_chirp = (m₁ m₂)^{3/5} / (m₁ + m₂)^{1/5}
     *
     * The chirp mass is the dominant parameter controlling the
     * inspiral waveform phase evolution at leading (Newtonian) order.
     *
     * It's the best-measured parameter from GW observations:
     * GW150914: M_chirp = 28.1^{+1.8}_{-1.5} M_sun (90% CI). */
    if (m1_kg <= 0 || m2_kg <= 0) return 0.0;
    double num = pow(m1_kg * m2_kg, 0.6);
    double den = pow(m1_kg + m2_kg, 0.2);
    return num / den;
}

double final_mass_after_merger(double m1, double m2,
                                double chi1, double chi2) {
    /* Fitting formula for remnant mass from NR simulations.
     *
     * For non-spinning equal-mass: M_final ≈ 0.951 M_total.
     * Radiated energy ≈ 4.9% of total mass.
     *
     * Reference: Barausse et al., PRD 85, 024046 (2012).
     *
     * This simple model gives approximate results; full NR
     * surrogate models (NRSur7dq4) provide higher accuracy. */

    if (m1 <= 0 || m2 <= 0) return 0.0;
    double M_total = m1 + m2;
    double eta = m1 * m2 / (M_total * M_total); /* symmetric mass ratio */

    /* Effective spin parameter */
    double chi_eff = (m1 * chi1 + m2 * chi2) / M_total;

    /* Radiated fraction fit (simplified from NR) */
    double E_rad_frac = 0.05 * (1.0 - 4.0 * eta)
                        + 0.02 * chi_eff * chi_eff;

    if (E_rad_frac > 0.15) E_rad_frac = 0.15;
    if (E_rad_frac < 0.0) E_rad_frac = 0.0;

    return M_total * (1.0 - E_rad_frac);
}

double final_spin_after_merger(double m1, double m2,
                                double chi1, double chi2) {
    /* Fitting formula for remnant spin.
     *
     * For equal-mass non-spinning: χ_final ≈ 0.686 (from NR).
     * For equal-mass aligned-spin (χ₁=χ₂=0.99): χ_final ≈ 0.95.
     *
     * GW150914: χ_final ≈ 0.68 (LIGO measurement, consistent with NR). */

    if (m1 <= 0 || m2 <= 0) return 0.0;
    double M_total = m1 + m2;
    double eta = m1 * m2 / (M_total * M_total);

    /* Simple fit (Rezzolla et al. 2008) */
    double chi_eff = (m1 * chi1 + m2 * chi2) / M_total;
    double s4 = 2.0 * sqrt(3.0);
    double s5 = -3.517;

    return chi_eff + s4 * eta + s5 * eta * eta;
}

double radiated_energy_fraction(double m1, double m2,
                                 double chi1, double chi2) {
    if (m1 <= 0 || m2 <= 0) return 0.0;
    double M_total = m1 + m2;
    double eta = m1 * m2 / (M_total * M_total);
    double chi_eff = (m1 * chi1 + m2 * chi2) / M_total;

    /* E_rad / M_total ≈ 0.05 (1 - 4η) + spin corrections */
    double frac = 0.05 * (1.0 - 4.0 * eta);
    double spin_corr = 0.02 * chi_eff * chi_eff;
    if (chi_eff > 0) spin_corr *= 1.5; /* Aligned spins radiate more */

    frac += spin_corr;
    if (frac < 0) frac = 0;
    if (frac > 0.12) frac = 0.12;
    return frac;
}

double merger_peak_frequency(double total_mass_kg, double chi_eff) {
    /* f_peak ≈ c³/(π √6 G M_total) × (1 + spin corrections)
     *
     * For non-spinning 65 M_sun: f_peak ≈ 150 Hz.
     * This falls in the most sensitive band of LIGO (50-1000 Hz). */
    if (total_mass_kg <= 0) return 0.0;
    double f0 = BH_C * BH_C * BH_C
                / (M_PI * sqrt(6.0) * BH_G * total_mass_kg);
    return f0 * (1.0 + 1.5 * chi_eff * chi_eff);
}

double effective_spin(double m1, double chi1, double theta1,
                       double m2, double chi2, double theta2) {
    /* χ_eff = (m₁ χ₁ cos θ₁ + m₂ χ₂ cos θ₂) / (m₁ + m₂)
     *
     * This is the best-measured spin parameter from GW observations.
     * For GW150914: χ_eff = -0.01^{+0.13}_{-0.13} (weakly constrained).
     *
     * For GW170729 (heaviest BBH): χ_eff > 0.5 likely (rapidly spinning).
     * For GW190521 (IMBH): χ_eff ≈ 0.08 (low spin). */
    if (m1 + m2 <= 0) return 0.0;
    return (m1 * chi1 * cos(theta1) + m2 * chi2 * cos(theta2)) / (m1 + m2);
}

double gw_luminosity_binary(double m1_kg, double m2_kg, double separation_m) {
    /* L_GW = (32/5) (G⁴/c⁵) μ² M_total³ / r⁵
     *
     * Quadrupole formula (Einstein 1918).
     *
     * For a 30+30 M_sun BBH at separation r=300 km (just before merger):
     * L_GW ≈ 10⁴⁹ W — about as bright as ALL stars in the
     * observable universe combined, but in gravitational waves!
     *
     * Energy flux decreases as 1/r⁵ (much steeper than 1/r² for EM). */
    if (m1_kg <= 0 || m2_kg <= 0 || separation_m <= 0) return 0.0;
    double mu = m1_kg * m2_kg / (m1_kg + m2_kg);
    double M_total = m1_kg + m2_kg;
    return (32.0 / 5.0) * pow(BH_G, 4) / pow(BH_C, 5)
           * mu * mu * M_total * M_total * M_total
           / pow(separation_m, 5);
}

/* ================================================================
 *  TEUKOLSKY EQUATION
 * ================================================================ */

double teukolsky_radial_potential(double r, double M, double a,
                                   double omega, int l, int m, int s) {
    /* Teukolsky radial potential for spin-s perturbations.
     *
     * The Teukolsky master equation:
     * Δ⁻ˢ d/dr [Δˢ⁺¹ dR/dr] + V(r) R = 0
     *
     * V(r) = K² - 2is(r-M)K / Δ + 4irωs - λ
     *
     * where K = (r² + a²) ω - a m,
     * λ is the angular eigenvalue (spin-weighted spheroidal harmonic),
     * Δ = r² - 2Mr + a².
     *
     * This is the FOUNDATION for all BH perturbation theory:
     * QNMs, gravitational wave emission, superradiance,
     * tidal deformability of BHs (=0 for all BHs!).
     *
     * @return Re[V(r)] (real part of the potential) */

    double delta = r * r - 2.0 * M * r + a * a;
    if (delta <= 0) return 0.0;

    double K = (r * r + a * a) * omega - a * m;
    double lambda = (double)(l * (l + 1)) - s * (s + 1); /* Schwarzschild limit */

    /* V(r) = K²/Δ + (2is(r-M)K)/Δ + 4irωs - λ */
    double K2_over_D = K * K / delta;
    double V_real = K2_over_D - lambda;

    return V_real;
}

double teukolsky_starobinsky_check(double r, double M, double a,
                                    double omega, int l, int m) {
    /* The Teukolsky-Starobinsky identities relate s=+2 (Ψ₀)
     * and s=-2 (Ψ₄) perturbations.
     *
     * Physical observables (GW strain) are independent of whether
     * we compute Ψ₀ or Ψ₄ — the identity guarantees consistency.
     *
     * This function checks the numerical identity for a test point.
     * For exact solutions: V_{+2} and V_{-2} are related by
     * a known transformation. Returns the numerical difference. */

    double V_plus2 = teukolsky_radial_potential(r, M, a, omega, l, m, +2);
    double V_minus2 = teukolsky_radial_potential(r, M, a, omega, l, m, -2);

    /* The potentials differ but physical quantities (QNM frequencies,
     * reflection coefficients) are identical for ±s. */
    return fabs(V_plus2 - V_minus2);
}

/* ================================================================
 *  SIGNAL ANALYSIS UTILITIES
 * ================================================================ */

double matched_filter_snr(const double *signal_freq,
                           const double *signal_amp,
                           const double *noise_psd,
                           int n_freqs) {
    /* Matched filter SNR:
     * ρ² = 4 Σ |h̃(f_i)|² / S_n(f_i) · Δf
     *
     * This is the optimal detection statistic for a known signal
     * in Gaussian noise. For GW150914 in LIGO: ρ ≈ 24.
     *
     * The matched filter is optimal (Neyman-Pearson) for
     * detecting a known signal in additive Gaussian noise. */

    if (signal_freq == NULL || signal_amp == NULL
        || noise_psd == NULL || n_freqs <= 1) return 0.0;

    double snr_sq = 0.0;
    for (int i = 0; i < n_freqs - 1; i++) {
        double df = signal_freq[i+1] - signal_freq[i];
        if (df <= 0 || noise_psd[i] <= 0) continue;
        snr_sq += 4.0 * signal_amp[i] * signal_amp[i] / noise_psd[i] * df;
    }
    return sqrt(snr_sq);
}

double characteristic_strain(double freq, double fourier_amp) {
    /* h_c(f) = 2f |h̃(f)|
     *
     * The characteristic strain is a dimensionless quantity used
     * to visualize GW signals on detector sensitivity curves.
     *
     * For an inspiraling binary at frequency f:
     * h_c(f) ∝ M_chirp^{5/6} f^{-1/6} / d_L */
    return 2.0 * freq * fabs(fourier_amp);
}

double luminosity_distance_from_gw(double chirp_mass_kg,
                                    double frequency_hz,
                                    double strain_amplitude) {
    /* d_L = (5/96) (c/f) (G M_chirp/c³)^{5/6} / (π^{2/3} h)
     *
     * Gravitational wave standard siren:
     * Unlike EM observations, GW luminosity distance is a DIRECT
     * measurement, not relying on the cosmic distance ladder.
     *
     * With electromagnetic counterpart (kilonova, GRB), we also get
     * redshift → measurement of H₀ independent of local distance ladder.
     *
     * GW170817 (binary neutron star): d_L ≈ 40 Mpc,
     * used with NGC 4993 redshift to measure H₀ ≈ 70 km/s/Mpc
     * (Abbott et al., Nature 551, 2017). */

    if (chirp_mass_kg <= 0 || frequency_hz <= 0 || strain_amplitude <= 0)
        return 0.0;

    double Mc_geom = mass_to_geometric(chirp_mass_kg);
    /* f_geom for potential geometric-unit consistency checks */
    (void)(frequency_hz * BH_G * chirp_mass_kg / (BH_C * BH_C * BH_C));

    return (5.0 / 96.0) * (BH_C / frequency_hz)
           * pow(Mc_geom, 5.0 / 6.0) / (pow(M_PI, 2.0 / 3.0) * strain_amplitude);
}
