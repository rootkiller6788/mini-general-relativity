/**
 * @file    gw_detection.h
 * @brief   GW Detection Theory — Matched Filtering, SNR, Fisher Matrix
 *
 * Reference: Maggiore Vol.1 Ch.7, Cutler & Flanagan (1994),
 *            Moore, Cole & Berry (2015), LIGO-Virgo data analysis guide
 *
 * L2 — Matched filtering, optimal SNR
 * L3 — Fisher information matrix, parameter estimation
 * L5 — Inner product, chi-square discriminator
 * L7 — LIGO/Virgo noise curves, detection statistics
 */

#ifndef GW_DETECTION_H
#define GW_DETECTION_H

#include "gw_core.h"

/* ================================================================
 * L2 — Noise-Weighted Inner Product
 * ================================================================ */

/**
 * Noise-weighted inner product (frequency domain):
 *
 *   (a | b) = 4 Re integral_{f_min}^{f_max} [ a*(f) b(f) / S_n(f) ] df
 *
 * @param a_re, a_im  Real/imag parts of signal a(f), length nf
 * @param b_re, b_im  Real/imag parts of template b(f), length nf
 * @param psd         One-sided noise PSD S_n(f), length nf
 * @param nf          Number of frequency bins
 * @param df          Frequency resolution [Hz]
 * @return Inner product value (dimensionless)
 */
double gw_inner_product(const double *a_re, const double *a_im,
                        const double *b_re, const double *b_im,
                        const double *psd, int nf, double df);

/**
 * Signal-to-noise ratio (SNR) for template h against data d:
 *
 *   SNR = (d | h) / sqrt( (h | h) )
 *
 * Optimal SNR (when h is the true signal): SNR_opt = sqrt( (h|h) )
 *
 * @param d_re, d_im  Data in frequency domain
 * @param h_re, h_im  Template in frequency domain
 * @param psd         PSD
 * @param nf          Number of frequency bins
 * @param df          Frequency resolution [Hz]
 * @return SNR
 */
double gw_snr(const double *d_re, const double *d_im,
              const double *h_re, const double *h_im,
              const double *psd, int nf, double df);

/**
 * Optimal SNR (perfect template match, d = h):
 *
 *   rho_opt = sqrt( (h|h) )
 */
double gw_optimal_snr(const double *h_re, const double *h_im,
                      const double *psd, int nf, double df);

/* ================================================================
 * L3 — Fisher Matrix for Parameter Estimation
 * ================================================================ */

/**
 * Fisher information matrix element for parameters theta_i, theta_j:
 *
 *   Gamma_{ij} = ( partial_i h | partial_j h )
 *
 * Evaluated by finite difference at fiducial parameter values.
 *
 * @param Gamma         Output Fisher matrix [nparam x nparam], row-major
 * @param nparam        Number of parameters
 * @param h_func        Waveform function: h(f; params) -> re[], im[]
 * @param params        Fiducial parameter values [nparam]
 * @param dparams       Finite-difference steps [nparam]
 * @param psd           Noise PSD
 * @param nf, df        Frequency grid
 * @param work_re, ...  Work arrays (6 * nf doubles)
 * @return 0 on success, -1 on error
 */
int gw_fisher_matrix(double *Gamma, int nparam,
    void (*h_func)(double *re, double *im, const double *params,
                   const double *f, int nf, double df),
    const double *params, const double *dparams,
    const double *psd, int nf, double df,
    double *work_re, double *work_im,
    double *work2_re, double *work2_im);

/**
 * Parameter uncertainty from Fisher matrix (1-sigma):
 *
 *   sigma_i = sqrt( (Gamma^{-1})_{ii} )
 *
 * @param sigma     Output uncertainties [nparam]
 * @param Gamma     Fisher matrix [nparam x nparam] row-major
 * @param nparam    Number of parameters
 * @return 0 on success, -1 if matrix singular
 */
int gw_fisher_uncertainties(double *sigma, const double *Gamma, int nparam);

/**
 * Covariance matrix = Fisher^{-1}.
 *
 * @param cov    Output covariance [nparam x nparam] row-major
 * @param Gamma  Fisher matrix
 * @param nparam Number of parameters
 * @return 0 on success, -1 if singular
 */
int gw_fisher_covariance(double *cov, const double *Gamma, int nparam);

/* ================================================================
 * L3 — Matrix Operations for Detection
 * ================================================================ */

/**
 * Invert a symmetric positive-definite matrix (Cholesky-based).
 *
 * @param A_inv  Output inverse [n x n] row-major
 * @param A      Input matrix [n x n] row-major
 * @param n      Dimension
 * @return 0 on success, -1 if not positive-definite
 */
int gw_matrix_invert_cholesky(double *A_inv, const double *A, int n);

/**
 * Solve A x = b for symmetric positive-definite A (Cholesky).
 *
 * @param x  Output solution [n]
 * @param A  SPD matrix [n x n] row-major
 * @param b  Right-hand side [n]
 * @param n  Dimension
 * @return 0 on success, -1 on failure
 */
int gw_cholesky_solve(double *x, const double *A, const double *b, int n);

/* ================================================================
 * L5 — Matched Filter in Time Domain
 * ================================================================ */

/**
 * Matched filter output (time domain) via cross-correlation with template.
 *
 * Output z(t) = integral d(tau) h(tau - t) / S_n integrated appropriately.
 * Simplified: z[k] = sum_j d[j] * h[j-k]
 *
 * @param z          Output matched filter time series [n_data]
 * @param data       Data time series [n_data]
 * @param template   Template time series [n_template] (n_template <= n_data)
 * @param n_data     Length of data
 * @param n_template Length of template
 */
void gw_matched_filter_td(double *z, const double *data,
                          const double *template_,
                          int n_data, int n_template);

/**
 * Peak SNR from matched filter output.
 *
 * @param z     Matched filter output
 * @param n     Length
 * @param peak_index Output: index of peak
 * @return Peak SNR value
 */
double gw_matched_filter_peak(const double *z, int n, int *peak_index);

/* ================================================================
 * L2 — LIGO/Virgo Noise Models
 * ================================================================ */

/**
 * LIGO design sensitivity PSD (approximate analytic fit).
 *
 * S_n(f) for aLIGO design, valid roughly f in [10, 2000] Hz.
 *
 * @param f   Frequency [Hz]
 * @return S_n(f) [Hz^{-1}]
 */
double gw_psd_aligo_design(double f);

/**
 * LIGO O3 actual sensitivity PSD (approximate).
 */
double gw_psd_aligo_o3(double f);

/**
 * LISA sensitivity PSD (approximate analytic fit).
 */
double gw_psd_lisa(double f);

/**
 * Einstein Telescope (ET) sensitivity PSD (approximate).
 */
double gw_psd_et(double f);

/* ================================================================
 * L2 — Horizon Distance
 * ================================================================ */

/**
 * Horizon distance: the maximum distance at which an optimally oriented
 * equal-mass binary would produce SNR = 8.
 *
 * Computed by root-finding: find D_L such that SNR(D_L) = 8.
 *
 * @param Mc     Chirp mass in detector frame [kg]
 * @param psd_func PSD function pointer
 * @param f_min  Low-frequency cutoff [Hz]
 * @param f_max  High-frequency cutoff [Hz]
 * @param nf     Frequency grid points
 * @return Horizon distance [m]
 */
double gw_horizon_distance(double Mc,
    double (*psd_func)(double f),
    double f_min, double f_max, int nf);

/* ================================================================
 * L7 — False-Alarm and Detection Statistics
 * ================================================================ */

/**
 * False-alarm probability (FAP) for a given SNR threshold,
 * assuming Gaussian noise:
 *
 *   FAP(SNR_thr) = exp(-SNR_thr^2 / 2)
 *
 * (Single-template, known phase case.)
 */
double gw_false_alarm_probability(double snr_threshold);

/**
 * False-alarm rate (FAR) given FAP and number of independent trials:
 *
 *   FAR = FAP * N_trials / T_obs
 */
double gw_false_alarm_rate(double fap, double N_trials, double T_obs);

#endif /* GW_DETECTION_H */
