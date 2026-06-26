/** @file nr_wave.h
 * @brief Gravitational wave extraction in numerical relativity.
 *
 * References: Newman & Penrose (1962), Baumgarte & Shapiro (2010), Ch. 9.
 *
 * Knowledge Coverage:
 *   L1: Weyl scalar Ψ_4, strain h_+/h_×, spin-weighted harmonics
 *   L2: Newman-Penrose formalism for wave extraction
 *   L3: Spin-weighted harmonics, tetrad formalism
 *   L5: Interpolation to extraction spheres, mode decomposition
 *   L7: Gravitational wave templates (LIGO/Virgo)
 */

#ifndef NR_WAVE_H
#define NR_WAVE_H

#include "nr_grid.h"
#include "nr_adm.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum l-mode for wave extraction */
#define NR_WAVE_LMAX 8

/** Maximum number of time samples for wave strain time series */
#define NR_WAVE_MAX_TIME 4096

/* ===========================================================================
 * L1: Weyl Scalar & Strain Modes
 *
 * Complex numbers nr_complex_t are defined in nr_grid.h.
 * =========================================================================== */

/**
 * L1: Weyl scalar Ψ_4 mode coefficients.
 *
 * Ψ_4(t, θ, φ) = Σ_{l=2}^{Lmax} Σ_{m=-l}^{l} Ψ_4^{lm}(t) _{-2}Y_{lm}(θ, φ)
 */
typedef struct {
    int          lmax;
    nr_complex_t psi4_mode[(NR_WAVE_LMAX+1)*(2*NR_WAVE_LMAX+1)];
} nr_psi4_modes_t;

/**
 * L1: Gravitational wave strain mode time series.
 * h_{lm}(t) = h_+^{lm}(t) − i h_×^{lm}(t)
 */
typedef struct {
    int          lmax;
    int          ntimes;
    double       dt;
    double       t0;
    nr_complex_t *hlm_data;
} nr_strain_series_t;

/* ===========================================================================
 * L2: Newman-Penrose Tetrad
 *
 * The NP tetrad on a spatial hypersurface:
 *   l^μ = (1/√2) (n^μ + s^μ)     outgoing null (approximate)
 *   n^μ = (1/√2) (n^μ − s^μ)     ingoing null (approximate)
 *   m^μ = (1/√2) (e_1^μ + i e_2^μ)  complex null on S^2
 *
 * where n^μ is the timelike unit normal to Σ_t,
 * s^i is the unit radial vector on Σ_t,
 * and (e_1^i, e_2^i) span the tangent space to the extraction sphere.
 * =========================================================================== */

/**
 * Construct the NP tetrad at a grid point for wave extraction.
 *
 * @param s_i[3]     Outward unit normal to extraction sphere (in)
 * @param e1_i[3]    First tangent vector to sphere (in)
 * @param e2_i[3]    Second tangent vector to sphere (in)
 * @param mi[3]      Output: complex m^i vector (complex spatial components)
 *
 * Reference: Baumgarte & Shapiro (2010), Sec. 9.2
 * Complexity: O(1).
 */
void nr_wave_tetrad_m(const double s_i[3],
                       const double e1_i[3],
                       const double e2_i[3],
                       nr_complex_t mi[3]);

/* ===========================================================================
 * L2: Ψ_4 Computation
 * =========================================================================== */

/**
 * Compute Ψ_4 at a single point from 3+1 variables.
 *
 * Ψ_4 = C_{αβγδ} n^α m̄^β n^γ m̄^δ
 *
 * In 3+1 decomposition (using Gauss-Codazzi relations):
 *   Ψ_4 = (R_{ij} − K K_{ij} + K_{ik} K^k_j) m̄^i m̄^j
 *         + i ε_i^{kl} D_k K_{jl} m̄^i m̄^j
 *
 * where m^i is the complex null vector from nr_wave_tetrad_m.
 *
 * @param gamma     3-metric
 * @param K         Extrinsic curvature
 * @param i,j,k     Grid indices
 * @param dx        Grid spacing
 * @param mi[3]     Complex m^i vector
 * @return          Ψ_4 at the point
 *
 * Reference: Baker, Campanelli, Lousto, PRD 65, 044001 (2002), Eq. (10)
 * Complexity: O(1).
 */
nr_complex_t nr_wave_psi4_point(const nr_sym_tensor3_t *gamma,
                                 const nr_sym_tensor3_t *K,
                                 int i, int j, int k, double dx,
                                 const nr_complex_t mi[3]);

/**
 * Compute Ψ_4 on an extraction sphere at radius r_ext.
 *
 * Evaluates Ψ_4 at (N_θ × N_φ) uniformly spaced points on the sphere
 * of radius r_ext centered at (cx, cy, cz). Uses trilinear interpolation
 * to get 3+1 variables at each angular point.
 *
 * @param gamma     3-metric
 * @param K         Extrinsic curvature
 * @param grid      Grid structure
 * @param dx        Grid spacing
 * @param r_ext     Extraction radius
 * @param cx,cy,cz  Center of extraction sphere
 * @param n_theta   Polar resolution
 * @param n_phi     Azimuthal resolution
 * @param psi4_grid Output: Ψ_4 values on the angular grid [n_theta × n_phi]
 *
 * Complexity: O(N_θ × N_φ).
 */
void nr_wave_psi4_sphere(const nr_sym_tensor3_t *gamma,
                          const nr_sym_tensor3_t *K,
                          const nr_grid_t *grid, double dx,
                          double r_ext,
                          double cx, double cy, double cz,
                          int n_theta, int n_phi,
                          nr_complex_t *psi4_grid);

/* ===========================================================================
 * L3: Spin-Weighted Spherical Harmonic Decomposition
 * =========================================================================== */

/**
 * Compute spin-weighted spherical harmonic _{-2}Y_{lm}(θ, φ).
 *
 * The spin-weight s = −2 harmonic is used for Ψ_4 decomposition:
 *   Ψ_4^{lm} = ∫ Ψ_4(θ, φ) _{-2}Ȳ_{lm}(θ, φ) dΩ
 *
 * Computed via recurrence relations from the scalar Y_{lm}.
 *
 * @param l     Harmonic degree (l ≥ 2)
 * @param m     Harmonic order (−l ≤ m ≤ l)
 * @param theta Polar angle
 * @param phi   Azimuthal angle
 * @return      _{-2}Y_{lm}(θ, φ)
 *
 * Reference: Goldberg et al., JMP 8, 2155 (1967)
 * Complexity: O(1) using pre-computed coefficients.
 */
nr_complex_t nr_wave_spin_weighted_Ylm(int l, int m,
                                        double theta, double phi);

/**
 * Decompose Ψ_4 on the extraction sphere into spin-weight −2 harmonics.
 *
 * Ψ_4^{lm} = ∫ Ψ_4(θ, φ) _{-2}Ȳ_{lm}(θ, φ) dΩ
 *
 * Uses trapezoidal integration on the angular grid.
 *
 * @param psi4_grid    Ψ_4 on angular grid [n_theta × n_phi]
 * @param n_theta, n_phi  Angular grid dimensions
 * @param lmax         Maximum l for decomposition
 * @param modes        Output: mode coefficients Ψ_4^{lm}
 *
 * Complexity: O(N_θ × N_φ × Lmax²).
 */
void nr_wave_decompose_psi4(const nr_complex_t *psi4_grid,
                             int n_theta, int n_phi, int lmax,
                             nr_psi4_modes_t *modes);

/* ===========================================================================
 * L5: Strain Computation (Double Time Integration)
 * =========================================================================== */

/**
 * Accumulate Ψ_4^{lm}(t) into a strain time series via fixed-frequency
 * integration (FFI) to remove secular drift.
 *
 * In the frequency domain:
 *   h̃_{lm}(ω) = −Ψ̃_4^{lm}(ω) / ω²   (for ω ≠ 0)
 *
 * Low-frequency components below ω_cut are suppressed.
 *
 * @param psi4_lm     Current Ψ_4^{lm} (complex)
 * @param l, m        Mode indices
 * @param t           Current time
 * @param dt          Time step
 * @param omega_cut   Cutoff frequency for FFI
 * @param series      Accumulating strain series (modified)
 *
 * Complexity: O(1) per call — uses running integrators.
 */
void nr_wave_integrate_psi4(nr_complex_t psi4_lm, int l, int m,
                             double t, double dt, double omega_cut,
                             nr_strain_series_t *series);

/**
 * Compute the plus and cross strain at (θ, φ) from mode coefficients.
 *
 * h_+(θ, φ) − i h_×(θ, φ) = Σ_{lm} h_{lm} _{-2}Y_{lm}(θ, φ)
 *
 * @param series     Strain mode time series
 * @param time_idx   Time index into the series
 * @param theta, phi Sky direction
 * @param hp         Output: h_+
 * @param hc         Output: h_×
 *
 * Complexity: O(Lmax²).
 */
void nr_wave_strain_at_angle(const nr_strain_series_t *series,
                              int time_idx, double theta, double phi,
                              double *hp, double *hc);

/**
 * Compute the total radiated energy from the strain modes.
 *
 * dE/dt = (r²/(16π)) Σ_{lm} |dh_{lm}/dt|²
 *
 * E_rad = ∫ (dE/dt) dt
 *
 * @param series     Strain mode time series
 * @return           Total radiated energy (in units where G=c=1)
 *
 * Reference: Baumgarte & Shapiro (2010), Eq. (9.137)
 * Complexity: O(N_times × Lmax²).
 */
double nr_wave_radiated_energy(const nr_strain_series_t *series);

/**
 * Compute the final black hole properties from the strain modes.
 *
 * From the dominant (l=2, m=2) mode, estimate:
 *   - Merger time t_peak (peak of |h_{22}|)
 *   - Ringdown frequency ω_{22} and damping time τ_{22}
 *
 * @param series     Strain mode series
 * @param t_peak     Output: merger time
 * @param omega22    Output: ringdown frequency
 * @param tau22      Output: ringdown damping time
 *
 * Complexity: O(N_times).
 */
void nr_wave_final_properties(const nr_strain_series_t *series,
                               double *t_peak, double *omega22,
                               double *tau22);

/* ===========================================================================
 * L7: GW150914-like Binary Black Hole Merger Analysis
 * =========================================================================== */

/**
 * Analyze a strain time series for LIGO/Virgo-relevant quantities.
 *
 * Computes:
 *   - Peak strain amplitude
 *   - Chirp mass M_chirp
 *   - Signal-to-noise ratio (simplified)
 *   - Total mass estimate
 *
 * @param series     Strain mode series
 * @param distance   Luminosity distance in Mpc
 * @param M_chirp    Output: chirp mass (solar masses)
 * @param M_total    Output: total mass (solar masses)
 * @param SNR        Output: signal-to-noise ratio
 *
 * Reference: Abbott et al., PRL 116, 061102 (2016) — GW150914
 * Complexity: O(N_times × Lmax²).
 */
void nr_wave_ligo_analysis(const nr_strain_series_t *series,
                            double distance,
                            double *M_chirp, double *M_total,
                            double *SNR);

/**
 * Allocate a strain time series.
 * Complexity: O(1).
 */
nr_strain_series_t* nr_strain_series_alloc(int lmax, int max_times);

/** Free a strain time series */
void nr_strain_series_free(nr_strain_series_t *series);

#ifdef __cplusplus
}
#endif

#endif /* NR_WAVE_H */
