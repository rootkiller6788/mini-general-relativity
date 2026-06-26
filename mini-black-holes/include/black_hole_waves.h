/**
 * @file black_hole_waves.h
 * @brief Quasinormal modes, ringdown, gravitational wave emission from BH mergers.
 *
 * Reference: Teukolsky 1973; Press & Teukolsky 1973; Kokkotas & Schmidt 1999;
 *            Berti, Cardoso, Will 2006; LIGO-Virgo 2016
 *
 * L1 Definitions: QNM frequencies, ringdown waveform
 * L5 Computational Methods: Continued fraction method for QNM, matched filtering
 * L6 Canonical Systems: Schwarzschild/Kerr QNM spectrum
 * L7 Applications: GW150914 parameter estimation
 */

#ifndef BLACK_HOLE_WAVES_H
#define BLACK_HOLE_WAVES_H

#include "black_hole_metrics.h"

/* ================================================================
 *  QUASINORMAL MODES
 * ================================================================ */

/**
 * @brief A single quasinormal mode (QNM).
 *
 * QNMs are the characteristic "ringing" frequencies of a perturbed BH.
 * They are complex: ω = ω_R + i ω_I.
 * - ω_R: oscillation frequency
 * - ω_I < 0: damping rate (exponential decay)
 *
 * The waveform is h(t) ~ exp(i ω t) = exp(i ω_R t) exp(ω_I t)
 * = exp(i ω_R t) exp(-t/τ) where τ = -1/ω_I is the damping time.
 *
 * QNMs are the gravitational wave analog of atomic spectral lines —
 * they depend ONLY on the BH mass and spin (no-hair theorem).
 */
typedef struct {
    double omega_R;         /**< Real part: oscillation frequency [Hz or geometric] */
    double omega_I;         /**< Imaginary part: damping rate (negative for stable) */
    double quality_factor;  /**< Q = ω_R/(2|ω_I|) — number of oscillations in ringdown */
    int l;                  /**< Spherical harmonic index l (multipole) */
    int m;                  /**< Azimuthal index m */
    int n;                  /**< Overtone number n=0 (fundamental), n=1,2,... */
    int s;                  /**< Spin weight: 0=scalar, -1=EM, -2=gravitational */
} QuasinormalMode;

/**
 * @brief Fundamental (n=0, l=m=2) QNM frequency for Schwarzschild BH.
 *
 * M ω_{220} = 0.74734 + i (-0.17792)   (Chandrasekhar & Detweiler 1975)
 *
 * In SI units:
 * f = 0.74734 × c³/(2πGM) ≈ 12.07 kHz × (M_sun/M)
 * τ = M/(0.17792 c) ≈ 0.554 ms × (M/M_sun)
 *
 * For GW150914 (M_final ≈ 62 M_sun):
 * f ≈ 195 Hz, τ ≈ 4.3 ms — matches LIGO observation!
 */
QuasinormalMode schwarzschild_qnm_fundamental(int l, int m);

/**
 * @brief Fundamental QNM frequency for Kerr BH.
 *
 * Fitting formula (Berti, Cardoso, Will 2006):
 * ω(a) ≈ ω_{Schw} [1 + c₁(a/M) + c₂(a/M)² + ...]
 *
 * The QNM spectrum splits into prograde (m>0) and retrograde (m<0)
 * modes for rotating BHs.
 *
 * For l=m=2, the frequency INCREASES with spin for prograde modes:
 * Mω ≈ 0.75 (a=0) → 1.0 (a=0.9M) → 1.2 (a=0.99M)
 */
QuasinormalMode kerr_qnm_fundamental(double M, double a, int l, int m);

/**
 * @brief QNM frequency for the l=m=2 overtone n.
 * Overtones decay faster but contribute to early ringdown.
 */
QuasinormalMode schwarzschild_qnm_overtone(int l, int m, int n);

/**
 * @brief QNM for l=2, m=0 (spherical, non-rotating) mode.
 */
QuasinormalMode schwarzschild_qnm_l2m0(void);

/**
 * @brief Compute the quality factor Q of a QNM.
 * Q = ω_R / (2|ω_I|) ≈ number of cycles in the ringdown.
 * For GW150914: Q ≈ 5.5 (only a few cycles visible).
 */
double qnm_quality_factor(const QuasinormalMode *qnm);

/**
 * @brief Compute QNM frequency in Hz from geometric units.
 * f = ω_R × c³/(2π G M)
 */
double qnm_frequency_hz(const QuasinormalMode *qnm, double mass_kg);

/**
 * @brief Compute QNM damping time in seconds.
 * τ = -M / (ω_I × c) (geometric → SI)
 */
double qnm_damping_time(const QuasinormalMode *qnm, double mass_kg);

/* ================================================================
 *  RINGDOWN WAVEFORM
 * ================================================================ */

/**
 * @brief Ringdown waveform: h(t) = A exp(-t/τ) cos(2πf t + φ₀).
 *
 * This describes the gravitational wave signal from the final
 * BH after a binary merger. The ringdown is a superposition of
 * QNMs, dominated by the l=m=2 fundamental mode.
 *
 * @param t Time [s]
 * @param A Amplitude
 * @param f Frequency [Hz]
 * @param tau Damping time [s]
 * @param phi0 Initial phase [rad]
 * @return h(t)
 */
double ringdown_waveform(double t, double A, double f,
                          double tau, double phi0);

/**
 * @brief Multi-mode ringdown: sum over QNMs.
 *
 * h(t) = Σ_n A_n exp(-t/τ_n) cos(2πf_n t + φ_n)
 *
 * @param t Time
 * @param modes Array of QNM parameters
 * @param amps Array of amplitudes
 * @param phases Array of phases
 * @param n_modes Number of modes
 * @return h(t)
 */
double ringdown_multimode(double t,
                           const QuasinormalMode *modes,
                           const double *amps,
                           const double *phases,
                           int n_modes);

/* ================================================================
 *  BINARY BLACK HOLE MERGER PHENOMENOLOGY
 * ================================================================ */

/**
 * @brief Chirp mass: M_chirp = (m₁ m₂)^{3/5} / (m₁ + m₂)^{1/5}.
 *
 * The chirp mass is the dominant parameter governing the inspiral
 * waveform. It's the best-measured parameter from GW observations.
 *
 * For GW150914: M_chirp ≈ 28 M_sun.
 */
double chirp_mass(double m1_kg, double m2_kg);

/**
 * @brief Final mass after binary BH merger (fitting formula).
 *
 * Based on NR simulations (Barausse et al. 2012, Healy et al. 2014):
 * M_final ≈ M_total - E_radiated/c²
 *
 * Where radiated energy is ~3-5% of total mass for non-spinning,
 * up to ~10% for optimally spinning configurations.
 *
 * @param m1 First BH mass [kg]
 * @param m2 Second BH mass [kg]
 * @param chi1 Dimensionless spin of first BH
 * @param chi2 Dimensionless spin of second BH
 * @return Final mass [kg]
 */
double final_mass_after_merger(double m1, double m2,
                                double chi1, double chi2);

/**
 * @brief Final spin after binary BH merger (fitting formula).
 *
 * For equal-mass, non-spinning: χ_final ≈ 0.686 (consistent with NR).
 * For GW150914: χ_final ≈ 0.68 (LIGO measurement).
 *
 * @return Dimensionless spin parameter χ = a/M ∈ [0, 1]
 */
double final_spin_after_merger(double m1, double m2,
                                double chi1, double chi2);

/**
 * @brief Radiated energy in GWs during merger (fraction of total mass).
 *
 * E_rad/M_total ≈ 0.05 for non-spinning equal-mass,
 * up to ~0.11 for aligned-spin configurations.
 *
 * For GW150914: E_rad ≈ 3.0 M_sun c² radiated in ~0.2 s!
 * Peak luminosity: ~3.6 × 10⁴⁹ W (> all stars in observable universe).
 */
double radiated_energy_fraction(double m1, double m2,
                                 double chi1, double chi2);

/**
 * @brief Peak GW frequency at merger (approximate).
 *
 * f_peak ≈ c³/(π√6 GM_total) × function(spin)
 *
 * For GW150914 (M_total ≈ 65 M_sun): f_peak ≈ 150 Hz.
 */
double merger_peak_frequency(double total_mass_kg, double chi_eff);

/**
 * @brief Effective spin parameter for binary BH.
 *
 * χ_eff = (m₁ χ₁ cos θ₁ + m₂ χ₂ cos θ₂) / (m₁ + m₂)
 *
 * This is the spin combination best measured by GW detectors.
 * For GW150914: χ_eff ≈ -0.01 (consistent with zero, weakly constrained).
 */
double effective_spin(double m1, double chi1, double theta1,
                       double m2, double chi2, double theta2);

/**
 * @brief GW luminosity from a circular binary (quadrupole formula).
 *
 * L_GW = (32/5) (G⁴/c⁵) (μ² M_total³) / r⁵
 *
 * where μ = m₁m₂/(m₁+m₂) is the reduced mass and r is the separation.
 *
 * This drives the inspiral: energy loss → orbital decay → "chirp".
 */
double gw_luminosity_binary(double m1_kg, double m2_kg, double separation_m);

/* ================================================================
 *  TEUKOLSKY EQUATION
 * ================================================================ */

/**
 * @brief Teukolsky radial potential for spin-s perturbations of Kerr.
 *
 * The Teukolsky equation (1973) is the master equation for
 * linear perturbations of Kerr:
 *
 * Δ^{-s} d/dr [Δ^{s+1} dR/dr] + V(r) R = 0
 *
 * where R(r) is the radial Teukolsky function,
 * s = 0 (scalar), -1 (EM), -2 (gravitational).
 *
 * @param r Radial coordinate
 * @param M BH mass
 * @param a Spin parameter
 * @param omega Mode frequency
 * @param l Spherical harmonic index
 * @param m Azimuthal index
 * @param s Spin weight
 * @return Potential V(r)
 */
double teukolsky_radial_potential(double r, double M, double a,
                                   double omega, int l, int m, int s);

/**
 * @brief Teukolsky-Starobinsky identities: relate s=+2 and s=-2 perturbations.
 *
 * These identities guarantee that the physical degrees of freedom
 * for gravitational perturbations of Kerr are the same whether
 * described by Weyl scalar Ψ₀ (s=2) or Ψ₄ (s=-2).
 *
 * This function checks the identity for a test frequency.
 * @return Numerical violation (should be < 1e-10)
 */
double teukolsky_starobinsky_check(double r, double M, double a,
                                    double omega, int l, int m);

/* ================================================================
 *  SIGNAL ANALYSIS UTILITIES
 * ================================================================ */

/**
 * @brief Matched filter SNR for a GW signal.
 *
 * ρ² = 4 ∫₀^∞ |h̃(f)|² / S_n(f) df
 *
 * where S_n(f) is the detector noise power spectral density.
 * For LIGO design sensitivity, ρ ~ 24 for GW150914.
 *
 * @param signal_freq Array of frequencies
 * @param signal_amp Array of |h̃(f)|
 * @param noise_psd Array of S_n(f)
 * @param n_freqs Number of frequency bins
 * @return SNR ρ
 */
double matched_filter_snr(const double *signal_freq,
                           const double *signal_amp,
                           const double *noise_psd,
                           int n_freqs);

/**
 * @brief Characteristic strain for GW detection.
 * h_c(f) = 2f |h̃(f)| — dimensionless.
 * Compare to detector noise curve to assess detectability.
 */
double characteristic_strain(double freq, double fourier_amp);

/**
 * @brief Luminosity distance from GW amplitude.
 *
 * For a detected GW signal, the luminosity distance can be
 * measured directly (standard siren) without distance ladder:
 *
 * d_L = (5/96) (c/f) (G M_chirp/c³)^{5/6} / (π^{2/3} h)
 *
 * This is a primary science goal of LIGO-Virgo: measuring H₀
 * with GW standard sirens.
 */
double luminosity_distance_from_gw(double chirp_mass_kg,
                                    double frequency_hz,
                                    double strain_amplitude);

#endif /* BLACK_HOLE_WAVES_H */
