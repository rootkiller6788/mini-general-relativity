/**
 * @file    gw_detection.c
 * @brief   GW Detection — Matched Filtering, SNR, Fisher Matrix, Noise Curves
 *
 * L2 — Noise-weighted inner product, SNR
 * L3 — Fisher matrix, parameter estimation
 * L5 — Matched filter, Cholesky decomposition
 * L7 — LIGO/Virgo/LISA/ET noise curves, detection statistics
 */

#include "gw_detection.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * L2 — Noise-Weighted Inner Product
 *
 * (a|b) = 4 Re integral_{f_min}^{f_max} [ a*(f) b(f) / S_n(f) ] df
 *
 * Discrete: (a|b) ≈ 4 sum_j Re( a_j* b_j ) / S_n(f_j) * df
 * ================================================================ */

double gw_inner_product(const double *a_re, const double *a_im,
                        const double *b_re, const double *b_im,
                        const double *psd, int nf, double df) {
    double sum = 0.0;
    for (int i = 0; i < nf; i++) {
        if (psd[i] <= 0.0) continue;
        double re_part = a_re[i]*b_re[i] + a_im[i]*b_im[i];
        double im_part = a_re[i]*b_im[i] - a_im[i]*b_re[i];
        /* Re(a* b) = a_re*b_re + a_im*b_im */
        sum += re_part / psd[i];
    }
    return 4.0 * sum * df;
}

/* ================================================================
 * L2 — SNR
 *
 * SNR = (d|h) / sqrt((h|h))
 *
 * For matched filter, the SNR is the overlap between data and template,
 * normalized by the template norm.
 * ================================================================ */

double gw_snr(const double *d_re, const double *d_im,
              const double *h_re, const double *h_im,
              const double *psd, int nf, double df) {
    double d_dot_h = gw_inner_product(d_re, d_im, h_re, h_im, psd, nf, df);
    double h_dot_h = gw_inner_product(h_re, h_im, h_re, h_im, psd, nf, df);

    if (h_dot_h <= 0.0) return 0.0;
    return d_dot_h / sqrt(h_dot_h);
}

double gw_optimal_snr(const double *h_re, const double *h_im,
                      const double *psd, int nf, double df) {
    double h_dot_h = gw_inner_product(h_re, h_im, h_re, h_im, psd, nf, df);
    return sqrt(h_dot_h);
}

/* ================================================================
 * L3 — Fisher Matrix
 *
 * Gamma_{ij} = (partial_i h | partial_j h)
 *
 * Numerical differentiation: central finite difference
 *   partial_i h ≈ [ h(params + d_i) - h(params - d_i) ] / (2 * dparam_i)
 * ================================================================ */

int gw_fisher_matrix(double *Gamma, int nparam,
    void (*h_func)(double *re, double *im, const double *params,
                   const double *f, int nf, double df),
    const double *params, const double *dparams,
    const double *psd, int nf, double df,
    double *work_re, double *work_im,
    double *work2_re, double *work2_im) {
    if (nparam <= 0) return -1;

    /* Allocate arrays for partial derivatives */
    double **dh_re = (double**)malloc(nparam * sizeof(double*));
    double **dh_im = (double**)malloc(nparam * sizeof(double*));
    if (!dh_re || !dh_im) {
        free(dh_re); free(dh_im);
        return -1;
    }
    for (int i = 0; i < nparam; i++) {
        dh_re[i] = (double*)calloc(nf, sizeof(double));
        dh_im[i] = (double*)calloc(nf, sizeof(double));
        if (!dh_re[i] || !dh_im[i]) {
            for (int j = 0; j <= i; j++) {
                free(dh_re[j]); free(dh_im[j]);
            }
            free(dh_re); free(dh_im);
            return -1;
        }
    }

    /* Dummy f array for h_func (waveform functions take f array) */
    double *f_arr = (double*)malloc(nf * sizeof(double));
    if (!f_arr) {
        for (int i = 0; i < nparam; i++) { free(dh_re[i]); free(dh_im[i]); }
        free(dh_re); free(dh_im);
        return -1;
    }
    for (int j = 0; j < nf; j++) f_arr[j] = 0.0; /* caller's h_func should handle this */

    /* Compute numerical partial derivatives */
    double *params_plus  = (double*)malloc(nparam * sizeof(double));
    double *params_minus = (double*)malloc(nparam * sizeof(double));
    if (!params_plus || !params_minus) {
        for (int i = 0; i < nparam; i++) { free(dh_re[i]); free(dh_im[i]); }
        free(dh_re); free(dh_im); free(f_arr);
        free(params_plus); free(params_minus);
        return -1;
    }

    for (int i = 0; i < nparam; i++) {
        memcpy(params_plus,  params, nparam * sizeof(double));
        memcpy(params_minus, params, nparam * sizeof(double));
        params_plus[i]  += dparams[i];
        params_minus[i] -= dparams[i];

        h_func(work_re, work_im, params_plus,  f_arr, nf, df);
        h_func(work2_re, work2_im, params_minus, f_arr, nf, df);

        double two_dp = 2.0 * dparams[i];
        for (int j = 0; j < nf; j++) {
            dh_re[i][j] = (work_re[j] - work2_re[j]) / two_dp;
            dh_im[i][j] = (work_im[j] - work2_im[j]) / two_dp;
        }
    }

    /* Gamma_{ij} = (dh_i | dh_j) */
    for (int i = 0; i < nparam; i++) {
        for (int j = 0; j <= i; j++) {
            double val = gw_inner_product(dh_re[i], dh_im[i],
                                          dh_re[j], dh_im[j],
                                          psd, nf, df);
            Gamma[i*nparam + j] = val;
            Gamma[j*nparam + i] = val;  /* symmetric */
        }
    }

    for (int i = 0; i < nparam; i++) { free(dh_re[i]); free(dh_im[i]); }
    free(dh_re); free(dh_im);
    free(f_arr);
    free(params_plus);
    free(params_minus);
    return 0;
}

/* ================================================================
 * L3 — Fisher Uncertainties & Covariance
 * ================================================================ */

int gw_fisher_uncertainties(double *sigma, const double *Gamma, int nparam) {
    double *Gamma_inv = (double*)malloc(nparam * nparam * sizeof(double));
    if (!Gamma_inv) return -1;

    int ret = gw_matrix_invert_cholesky(Gamma_inv, Gamma, nparam);
    if (ret != 0) {
        free(Gamma_inv);
        return -1;
    }

    for (int i = 0; i < nparam; i++) {
        double var = Gamma_inv[i*nparam + i];
        sigma[i] = (var >= 0.0) ? sqrt(var) : 0.0;
    }

    free(Gamma_inv);
    return 0;
}

int gw_fisher_covariance(double *cov, const double *Gamma, int nparam) {
    return gw_matrix_invert_cholesky(cov, Gamma, nparam);
}

/* ================================================================
 * L3 — Cholesky Decomposition and Matrix Inversion
 *
 * For symmetric positive-definite matrix A, decompose A = L L^T,
 * then invert via back-substitution.
 * ================================================================ */

static int cholesky_decompose(double *L, const double *A, int n) {
    /* Copy A to L (will be overwritten by lower triangle) */
    for (int i = 0; i < n*n; i++) L[i] = A[i];

    for (int j = 0; j < n; j++) {
        /* Diagonal */
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += L[j*n + k] * L[j*n + k];
        }
        double d = L[j*n + j] - sum;
        if (d <= 0.0) return -1;  /* Not positive-definite */
        L[j*n + j] = sqrt(d);

        /* Off-diagonal below */
        for (int i = j+1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i*n + k] * L[j*n + k];
            }
            L[i*n + j] = (L[i*n + j] - sum) / L[j*n + j];
        }
    }

    /* Zero out upper triangle */
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            L[i*n + j] = 0.0;

    return 0;
}

int gw_matrix_invert_cholesky(double *A_inv, const double *A, int n) {
    if (n <= 0) return -1;

    double *L = (double*)malloc(n * n * sizeof(double));
    if (!L) return -1;

    int ret = cholesky_decompose(L, A, n);
    if (ret != 0) { free(L); return -1; }

    /* Invert by solving n linear systems L L^T x = e_j
     * First solve L y = e_j, then L^T x = y
     * Store x as columns of A_inv */

    for (int j = 0; j < n; j++) {
        double y[64]; /* stack-allocated for small n */
        double *y_dyn = NULL;
        double *y_ptr = y;
        if (n > 64) {
            y_dyn = (double*)malloc(n * sizeof(double));
            y_ptr = y_dyn;
        }

        /* Forward substitution: L y = e_j */
        for (int i = 0; i < n; i++) {
            double sum = (i == j) ? 1.0 : 0.0;
            for (int k = 0; k < i; k++) {
                sum -= L[i*n + k] * y_ptr[k];
            }
            y_ptr[i] = sum / L[i*n + i];
        }

        /* Back substitution: L^T x = y */
        for (int i = n-1; i >= 0; i--) {
            double sum = y_ptr[i];
            for (int k = i+1; k < n; k++) {
                sum -= L[k*n + i] * A_inv[k*n + j];
            }
            A_inv[i*n + j] = sum / L[i*n + i];
        }

        if (y_dyn) free(y_dyn);
    }

    free(L);
    return 0;
}

int gw_cholesky_solve(double *x, const double *A, const double *b, int n) {
    double *L = (double*)malloc(n * n * sizeof(double));
    if (!L) return -1;

    int ret = cholesky_decompose(L, A, n);
    if (ret != 0) { free(L); return -1; }

    /* Forward: L y = b */
    double *y = (double*)malloc(n * sizeof(double));
    if (!y) { free(L); return -1; }

    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int k = 0; k < i; k++) sum -= L[i*n + k] * y[k];
        y[i] = sum / L[i*n + i];
    }

    /* Back: L^T x = y */
    for (int i = n-1; i >= 0; i--) {
        double sum = y[i];
        for (int k = i+1; k < n; k++) sum -= L[k*n + i] * x[k];
        x[i] = sum / L[i*n + i];
    }

    free(y);
    free(L);
    return 0;
}

/* ================================================================
 * L5 — Time-Domain Matched Filter
 * ================================================================ */

void gw_matched_filter_td(double *z, const double *data,
                          const double *template_,
                          int n_data, int n_template) {
    for (int k = 0; k < n_data; k++) {
        z[k] = 0.0;
    }
    for (int k = 0; k < n_data - n_template + 1; k++) {
        double sum = 0.0;
        for (int j = 0; j < n_template; j++) {
            sum += data[k + j] * template_[j];
        }
        z[k] = sum;
    }
}

double gw_matched_filter_peak(const double *z, int n, int *peak_index) {
    double max_val = z[0];
    int    max_idx = 0;
    for (int i = 1; i < n; i++) {
        if (z[i] > max_val) {
            max_val = z[i];
            max_idx = i;
        }
    }
    if (peak_index) *peak_index = max_idx;
    return max_val;
}

/* ================================================================
 * L7 — LIGO/Virgo/LISA/ET Noise PSD Models
 * ================================================================ */

double gw_psd_aligo_design(double f) {
    /* Analytic fit to aLIGO design sensitivity.
     * S_n(f) = S_0 [ (f/f_0)^{-41} + ... ] simplified.
     * Using the LIGO Scientific Collaboration parameterization (approx). */

    if (f < 5.0) f = 5.0;

    double x = f / 215.0;
    double S0 = 1.0e-49;  /* baseline */

    double sn = S0 * (pow(x, -4.14)
               - 5.0 * pow(x, -2.0)
               + 111.0 * (1.0 - x*x + 0.5*x*x*x*x) / (1.0 + 0.5*x*x));
    /* Clamp negative values */
    if (sn < 1e-52) sn = 1e-52;
    return sn;
}

double gw_psd_aligo_o3(double f) {
    /* O3 sensitivity ~ 1.5-2x worse than design at some frequencies */
    return gw_psd_aligo_design(f) * 1.5;
}

double gw_psd_lisa(double f) {
    /* LISA power spectral density (approximate).
     * Based on Robson, Cornish & Liu (2019). */

    double L = 2.5e9;  /* arm length [m] */
    double S_acc = 9.0e-30 * (1.0 + pow(0.0004/f, 2.0))
                   * (1.0 + pow(f/0.008, 4.0))
                   / pow(2.0 * GW_PI * f, 4.0);  /* acceleration noise */
    double S_oms = 2.25e-22 * pow(2.0*GW_PI*f / GW_C, 2.0);  /* position noise */

    double S_n = (10.0/3.0) * (S_acc + S_oms) / (L*L)
                 * (1.0 + 0.6 * pow(f*L/GW_C, 2.0));
    return S_n;
}

double gw_psd_et(double f) {
    /* Einstein Telescope sensitivity (approximate, ET-D sum). */
    if (f < 1.0) f = 1.0;

    double x = f / 100.0;
    double sn = 1.0e-50 * (2.0 * pow(x, -4.0) + 1.0 + 0.5 * pow(x, 2.0));
    return sn;
}

/* ================================================================
 * L7 — Horizon Distance
 * ================================================================ */

double gw_horizon_distance(double Mc,
    double (*psd_func)(double f),
    double f_min, double f_max, int nf) {
    if (nf < 2) return -1.0;

    double df = (f_max - f_min) / (nf - 1.0);

    /* Build fiducial waveform at D = 1 Mpc, then scale */
    double D_fid = 1.0 * GW_MPC;
    double *h_re = (double*)calloc(nf, sizeof(double));
    double *h_im = (double*)calloc(nf, sizeof(double));
    double *psd = (double*)calloc(nf, sizeof(double));
    if (!h_re || !h_im || !psd) {
        free(h_re); free(h_im); free(psd);
        return -1.0;
    }

    for (int i = 0; i < nf; i++) {
        double fi = f_min + i * df;
        psd[i] = psd_func(fi);

        /* Newtonian inspiral amplitude |h(f)| ~ A_Mc * f^{-7/6} / D_L */
        double A = GW_G_OVER_C4 * GW_G * Mc / (D_fid * GW_C * GW_C)
                   * pow(GW_G * Mc * GW_PI * fi / (GW_C*GW_C*GW_C), 2.0/3.0);
        /* Simplified: leading-order amplitude scaling */
        double amp = (GW_C / D_fid)
                     * sqrt(5.0/24.0) * pow(GW_PI, -2.0/3.0)
                     * pow(GW_G * Mc / (GW_C*GW_C*GW_C), 5.0/6.0)
                     * pow(fi, -7.0/6.0);
        h_re[i] = amp;
        h_im[i] = 0.0;
    }

    double snr_fid = gw_optimal_snr(h_re, h_im, psd, nf, df);
    double snr_target = 8.0;

    free(h_re); free(h_im); free(psd);

    if (snr_fid <= 0.0) return 1e30;
    return D_fid * snr_fid / snr_target;
}

/* ================================================================
 * L7 — Detection Statistics
 * ================================================================ */

double gw_false_alarm_probability(double snr_threshold) {
    return exp(-0.5 * snr_threshold * snr_threshold);
}

double gw_false_alarm_rate(double fap, double N_trials, double T_obs) {
    if (T_obs <= 0.0) return 1e100;
    return fap * N_trials / T_obs;
}
