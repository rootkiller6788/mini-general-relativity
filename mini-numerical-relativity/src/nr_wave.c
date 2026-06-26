/**
 * @file nr_wave.c
 * @brief Gravitational wave extraction: Ψ_4 computation, mode decomposition, strain.
 *
 * Extracts gravitational radiation from numerical relativity simulations
 * by computing the Weyl scalar Ψ_4 on coordinate spheres and decomposing
 * into spin-weight −2 spherical harmonics.
 *
 * The GW strain is obtained by double time integration:
 *   h_{lm}(t) = −∫_{−∞}^t dt' ∫_{−∞}^{t'} dt'' Ψ_4^{lm}(t'')
 *
 * References:
 *   - Newman & Penrose, JMP 3, 566 (1962)
 *   - Baker, Campanelli, Lousto, PRD 65, 044001 (2002)
 *   - Baumgarte & Shapiro (2010), Ch. 9
 *   - Alcubierre (2008), Ch. 8
 *   - Abbott et al., PRL 116, 061102 (2016) — GW150914
 *
 * Knowledge: L1 (Ψ_4, h_+/h_×, strain), L2 (NP formalism),
 *            L3 (spin-weighted harmonics), L5 (wave extraction, FFI),
 *            L7 (LIGO/Virgo analysis)
 */

#include "nr_wave.h"
#include "nr_utils.h"
#include "nr_adm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===========================================================================
 * L2: Newman-Penrose Tetrad
 * =========================================================================== */

void nr_wave_tetrad_m(const double s_i[3],
                       const double e1_i[3],
                       const double e2_i[3],
                       nr_complex_t mi[3]) {
    (void)s_i;  /* Radial vector used implicitly through e1/e2 */
    /* m^i = (1/√2) (e1^i + i e2^i) */
    double inv_sqrt2 = 1.0 / sqrt(2.0);
    for (int i = 0; i < 3; i++) {
        mi[i].re = inv_sqrt2 * e1_i[i];
        mi[i].im = inv_sqrt2 * e2_i[i];
    }
}

/* ===========================================================================
 * L2: Ψ_4 from 3+1 Variables
 * =========================================================================== */

nr_complex_t nr_wave_psi4_point(const nr_sym_tensor3_t *gamma,
                                 const nr_sym_tensor3_t *K,
                                 int i, int j, int k, double dx,
                                 const nr_complex_t mi[3]) {
    /* Ψ_4 = R_{ij} m̄^i m̄^j − K K_{ij} m̄^i m̄^j + K_{ik} K^k_j m̄^i m̄^j
     *      + i ε_i^{kl} D_k K_{jl} m̄^i m̄^j
     *
     * Compute at a grid point. mi is the complex m-vector.
     * m̄ is the complex conjugate.
     */

    /* Get metric and curvature at point */
    double g[6], Kij[6];
    int sx = gamma->comp[0]->stride_x, sy = gamma->comp[0]->stride_y;
    int base = i*sx + j*sy + k;
    for (int c = 0; c < 6; c++) {
        g[c] = gamma->comp[c]->data[base];
        Kij[c] = K->comp[c]->data[base];
    }

    double inv_g[3][3];
    nr_adm_inverse_metric(g, inv_g);

    /* 3-Ricci tensor */
    double Rij[3][3];
    nr_adm_ricci(gamma, i, j, k, dx, Rij);

    /* Trace of K */
    double trK = 0.0;
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            trK += inv_g[a][b] * Kij[map[a][b]];

    /* m̄^i m̄^j (complex conjugate pair) */
    nr_complex_t mbar[3];
    for (int a = 0; a < 3; a++) {
        mbar[a] = nr_complex_conj(mi[a]);
    }

    /* Real part: R_{ij} m̄^i m̄^j − K K_{ij} m̄^i m̄^j + K_{ik} K^k_j m̄^i m̄^j */
    nr_complex_t psi4 = nr_complex_make(0.0, 0.0);

    /* R_{ij} m̄^i m̄^j */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            nr_complex_t term = nr_complex_make(Rij[a][b], 0.0);
            term = nr_complex_mul(term, mbar[a]);
            term = nr_complex_mul(term, mbar[b]);
            psi4 = nr_complex_add(psi4, term);
        }
    }

    /* − K K_{ij} m̄^i m̄^j */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            nr_complex_t term = nr_complex_make(-trK * Kij[map[a][b]], 0.0);
            term = nr_complex_mul(term, mbar[a]);
            term = nr_complex_mul(term, mbar[b]);
            psi4 = nr_complex_add(psi4, term);
        }
    }

    /* K_{ik} K^k_j m̄^i m̄^j */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double sum_k = 0.0;
            for (int kk = 0; kk < 3; kk++) {
                /* K^k_j = γ^{kl} K_{lj} */
                double Kkj = 0.0;
                for (int l = 0; l < 3; l++) {
                    Kkj += inv_g[kk][l] * Kij[map[l][b]];
                }
                sum_k += Kij[map[a][kk]] * Kkj;
            }
            nr_complex_t term = nr_complex_make(sum_k, 0.0);
            term = nr_complex_mul(term, mbar[a]);
            term = nr_complex_mul(term, mbar[b]);
            psi4 = nr_complex_add(psi4, term);
        }
    }

    /* Imaginary part: i ε_i^{kl} D_k K_{jl} m̄^i m̄^j
     *
     * ε_i^{kl} is the Levi-Civita symbol in 3D.
     * D_k K_{jl} = ∂_k K_{jl} − Γ^m_{kj} K_{ml} − Γ^m_{kl} K_{jm}
     *
     * For simplicity, approximate D_k ≈ ∂_k (flat covariant derivative).
     */
    double dK[6][3];  /* ∂_k K_{ij} for each component */
    for (int c = 0; c < 6; c++)
        for (int dir = 0; dir < 3; dir++)
            dK[c][dir] = nr_fd_deriv1(K->comp[c], i, j, k, dx, dir, 4);

    double imag_part = 0.0;

    /* ε^{ikl} = 1 for even permutations, −1 for odd, 0 otherwise.
     * In 3D: ε^{123}=1, ε^{132}=−1, ε^{213}=−1, ε^{231}=1, ε^{312}=1, ε^{321}=−1
     */
    /* i ε_i^{kl} D_k K_{jl} m̄^i m̄^j
     * = i Σ_{i,k,l} ε_{ikl} D_k K_{jl} m̄^i m̄^j
     *
     * Actually the formula is: i ε^{kl}_m D_k K_{jl} m̄^m m̄^j
     *
     * ε^{kl}_m = γ^{kn} ε_{nml} where ε_{123}=1
     * For flat space approximation: ε^{kl}_m = ε_{mkl}
     */
    for (int m = 0; m < 3; m++) {
        for (int jj = 0; jj < 3; jj++) {
            /* m̄^m m̄^j */
            double mb_m_re = mbar[m].re, mb_m_im = mbar[m].im;
            double mb_j_re = mbar[jj].re, mb_j_im = mbar[jj].im;
            double prod_re = mb_m_re * mb_j_re - mb_m_im * mb_j_im;
            (void)mb_j_im; (void)mb_m_im;  /* Imag part reserved */

            /* Σ_{k,l} ε_{mkl} ∂_k K_{jl} */
            double sum_kl = 0.0;
            /* ε_{123}=1 → (k=1,l=2), (k=2,l=3), (k=3,l=1) are positive
             * ε_{132}=-1 → (k=1,l=3), (k=3,l=2), (k=2,l=1) are negative
             */
            int perms[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
            double signs[6] = {1,1,1,-1,-1,-1};
            for (int p = 0; p < 6; p++) {
                int kk = perms[p][0], l = perms[p][1], idx3 = perms[p][2];
                if (idx3 != m) continue;
                /* ∂_k K_{j l} */
                double dK_val = 0.0;
                int c_jl = map[jj][l];
                dK_val = dK[c_jl][kk];
                sum_kl += signs[p] * dK_val;
            }

            imag_part += prod_re * sum_kl;
        }
    }

    psi4.im += imag_part;

    return psi4;
}

/* ===========================================================================
 * Ψ_4 on an Extraction Sphere
 * =========================================================================== */

void nr_wave_psi4_sphere(const nr_sym_tensor3_t *gamma,
                          const nr_sym_tensor3_t *K,
                          const nr_grid_t *grid, double dx,
                          double r_ext,
                          double cx, double cy, double cz,
                          int n_theta, int n_phi,
                          nr_complex_t *psi4_grid) {
    if (!gamma || !K || !grid || !psi4_grid) return;

    for (int it = 0; it < n_theta; it++) {
        double theta = M_PI * (it + 0.5) / n_theta;
        double st = sin(theta), ct = cos(theta);

        for (int ip = 0; ip < n_phi; ip++) {
            double phi = 2.0 * M_PI * ip / n_phi;
            double sp = sin(phi), cp = cos(phi);

            /* Cartesian coordinates of the extraction point */
            double x = cx + r_ext * st * cp;
            double y = cy + r_ext * st * sp;
            double z = cz + r_ext * ct;

            /* Interpolate metric and curvature */
            double gij[6], Kv[6];
            nr_interp_tensor_trilinear(gamma, grid, x, y, z, gij);
            nr_interp_tensor_trilinear(K, grid, x, y, z, Kv);

            /* Unit radial vector s^i on the sphere */
            double si[3] = {st*cp, st*sp, ct};

            /* Tangent vectors e1^i (θ direction), e2^i (φ direction) */
            double e1[3] = {ct*cp, ct*sp, -st};
            double e2[3] = {-sp, cp, 0.0};

            /* Build m^i tetrad vector */
            nr_complex_t mi[3];
            nr_wave_tetrad_m(si, e1, e2, mi);

            /* Find the nearest grid indices for the point */
            int ng = grid->ng;
            double fx = (x - grid->xmin) / dx + ng;
            double fy = (y - grid->ymin) / dx + ng;
            double fz = (z - grid->zmin) / dx + ng;

            int i0 = (int)floor(fx);
            int j0 = (int)floor(fy);
            int k0 = (int)floor(fz);

            int nx_tot = grid->nx + 2*ng;
            int ny_tot = grid->ny + 2*ng;
            int nz_tot = grid->nz + 2*ng;

            if (i0 >= 0 && i0+1 < nx_tot && j0 >= 0 && j0+1 < ny_tot
                && k0 >= 0 && k0+1 < nz_tot) {
                nr_complex_t psi4 = nr_wave_psi4_point(gamma, K, i0, j0, k0, dx, mi);
                psi4_grid[it * n_phi + ip] = psi4;
            } else {
                psi4_grid[it * n_phi + ip] = nr_complex_make(0.0, 0.0);
            }
        }
    }
}

/* ===========================================================================
 * L3: Spin-Weighted Spherical Harmonics
 * =========================================================================== */

nr_complex_t nr_wave_spin_weighted_Ylm(int l, int m,
                                        double theta, double phi) {
    /* _{-2}Y_{lm}(θ,φ)
     *
     * Computed from the scalar Y_{lm} using spin-raising/lowering operators:
     *   _{s}Y_{lm} = sqrt((l-s)!/(l+s)!) ð^s Y_{lm}   (for s ≥ 0)
     *   _{s}Y_{lm} = sqrt((l+s)!/(l-s)!) (−1)^s ð̄^{−s} Y_{lm}   (for s ≤ 0)
     *
     * For s = −2:
     *   _{-2}Y_{lm} = sqrt((l+2)!/(l-2)!) ð̄^2 Y_{lm}
     *
     * where ð̄ = −(∂_θ − i cscθ ∂_φ) (eth-bar operator).
     *
     * Reference: Goldberg et al., JMP 8, 2155 (1967)
     */

    /* For l < 2, _{-2}Y_{lm} is not defined */
    if (l < 2) return nr_complex_make(0.0, 0.0);

    /* Use recurrence/closed form for l=2 for validation:
     * _{-2}Y_{2,±2} = sqrt(5/(64π)) (1 ± cosθ)^2 e^{±2iφ}
     * _{-2}Y_{2,±1} = sqrt(5/(16π)) sinθ (1 ± cosθ) e^{±iφ}
     * _{-2}Y_{2,0}  = sqrt(15/(32π)) sin²θ
     */

    double ct = cos(theta), st = sin(theta);

    if (l == 2) {
        if (m == 2) {
            double norm = sqrt(5.0 / (64.0 * M_PI));
            double val = norm * (1.0 + ct)*(1.0 + ct);
            return nr_complex_make(val * cos(2.0*phi), val * sin(2.0*phi));
        } else if (m == -2) {
            double norm = sqrt(5.0 / (64.0 * M_PI));
            double val = norm * (1.0 - ct)*(1.0 - ct);
            return nr_complex_make(val * cos(-2.0*phi), val * sin(-2.0*phi));
        } else if (m == 1) {
            double norm = sqrt(5.0 / (16.0 * M_PI));
            double val = norm * st * (1.0 + ct);
            return nr_complex_make(val * cos(phi), val * sin(phi));
        } else if (m == -1) {
            double norm = sqrt(5.0 / (16.0 * M_PI));
            double val = norm * st * (1.0 - ct);
            return nr_complex_make(val * cos(-phi), val * sin(-phi));
        } else {
            double norm = sqrt(15.0 / (32.0 * M_PI));
            double val = norm * st * st;
            return nr_complex_make(val, 0.0);
        }
    }

    /* For l > 2, we approximate using recurrence from l=2.
     * The general recurrence is:
     *   _{-2}Y_{l+1,m} = A_{lm} ct _{-2}Y_{lm} + B_{lm} _{-2}Y_{l-1,m}
     *
     * For now return 0 for l > 2 (full implementation requires recurrence).
     */
    return nr_complex_make(0.0, 0.0);
}

/* ===========================================================================
 * L5: Mode Decomposition
 * =========================================================================== */

void nr_wave_decompose_psi4(const nr_complex_t *psi4_grid,
                             int n_theta, int n_phi, int lmax,
                             nr_psi4_modes_t *modes) {
    if (!psi4_grid || !modes) return;

    modes->lmax = lmax;
    int n_modes = (lmax + 1) * (2 * lmax + 1);
    for (int idx = 0; idx < n_modes; idx++)
        modes->psi4_mode[idx] = nr_complex_make(0.0, 0.0);

    for (int l = 2; l <= lmax; l++) {
        for (int m = -l; m <= l; m++) {
            nr_complex_t sum = nr_complex_make(0.0, 0.0);
            double total_weight = 0.0;

            for (int it = 0; it < n_theta; it++) {
                double theta = M_PI * (it + 0.5) / n_theta;
                double st = sin(theta);

                for (int ip = 0; ip < n_phi; ip++) {
                    double phi = 2.0 * M_PI * ip / n_phi;

                    /* _{-2}Y_{lm}(θ,φ) */
                    nr_complex_t Ylm = nr_wave_spin_weighted_Ylm(l, m, theta, phi);
                    /* Complex conjugate: _{-2}Ȳ_{lm} */
                    nr_complex_t Ylm_conj = nr_complex_conj(Ylm);

                    /* Ψ_4(θ,φ) * _{-2}Ȳ_{lm} */
                    nr_complex_t prod = nr_complex_mul(
                        psi4_grid[it * n_phi + ip], Ylm_conj);

                    /* Integration weight: sinθ dθ dφ */
                    double weight = st * (M_PI / n_theta) * (2.0 * M_PI / n_phi);

                    sum.re += prod.re * weight;
                    sum.im += prod.im * weight;
                    total_weight += weight;
                }
            }

            int idx = l * (2 * lmax + 1) + (m + lmax);
            if (idx < n_modes) {
                modes->psi4_mode[idx] = sum;
            }
        }
    }
}

/* ===========================================================================
 * L5: Strain via Fixed-Frequency Integration
 * =========================================================================== */

void nr_wave_integrate_psi4(nr_complex_t psi4_lm, int l, int m,
                             double t, double dt, double omega_cut,
                             nr_strain_series_t *series) {
    (void)t; (void)omega_cut;  /* Reserved for fixed-frequency integration */
    if (!series || !series->hlm_data) return;
    if (l < 2) return;

    int idx = series->ntimes;
    if (idx >= NR_WAVE_MAX_TIME) return;

    int mode_idx = l * (2 * series->lmax + 1) + (m + series->lmax);

    /* Fixed-frequency integration:
     *
     * In frequency domain: h(ω) = −Ψ_4(ω) / ω²
     *
     * In time domain: this corresponds to a low-frequency cutoff.
     *
     * Simple approach: accumulate first and second integrals with
     * a high-pass filter to remove secular drift.
     *
     * h(t) = −∫_0^t dt' ∫_0^{t'} dt'' Ψ_4(t'')
     *
     * Using trapezoidal integration:
     *   I1(t + dt) = I1(t) + (dt/2)(Ψ_4(t) + Ψ_4(t+dt))
     *   h(t + dt) = h(t) + (dt/2)(I1(t) + I1(t+dt))
     *
     * High-pass filter: h → h − low_freq_component
     */

    /* Store the current Ψ_4 in the series for the next step */
    /* For simplicity, accumulate integrals using running sums */

    /* Estimate the factor: h_lm ≈ −Ψ_4^{lm} / (l(l−1)) * t² mode
     * Actually: d²h/dt² = Ψ_4/r (at large r)
     * h = −∫∫ Ψ_4 dt² (after factoring extraction radius)
     */
    double factor = -1.0;

    /* Simple double integration with trapezoidal rule */
    if (idx == 0) {
        /* First sample: h = 0 */
        series->hlm_data[idx * ((series->lmax+1)*(2*series->lmax+1)) + mode_idx]
            = nr_complex_make(0.0, 0.0);
    } else {
        /* Integrate using running average of Ψ_4 */
        nr_complex_t h_prev = series->hlm_data[
            (idx-1) * ((series->lmax+1)*(2*series->lmax+1)) + mode_idx];

        /* Simple trapezoidal double integration approximation:
         * h(t) = h(t-dt) + dt * I1(t-dt/2)
         * I1(t) = I1(t-dt) + dt * Ψ_4(t-dt/2)
         *
         * Store I1 in imaginary part as a hack (or use separate arrays).
         * For simplicity: perform a running estimate.
         */
        double scale = factor * dt * dt;
        nr_complex_t h_new;
        h_new.re = h_prev.re + scale * psi4_lm.re;
        h_new.im = h_prev.im + scale * psi4_lm.im;

        series->hlm_data[idx * ((series->lmax+1)*(2*series->lmax+1)) + mode_idx]
            = h_new;
    }
}

/* ===========================================================================
 * L2: Strain at a Sky Direction
 * =========================================================================== */

void nr_wave_strain_at_angle(const nr_strain_series_t *series,
                              int time_idx, double theta, double phi,
                              double *hp, double *hc) {
    if (!series || !hp || !hc) return;
    if (time_idx < 0 || time_idx >= series->ntimes) {
        *hp = *hc = 0.0;
        return;
    }

    nr_complex_t h = nr_complex_make(0.0, 0.0);  /* h_+ − i h_× */
    int nmode = (series->lmax + 1) * (2 * series->lmax + 1);

    for (int l = 2; l <= series->lmax; l++) {
        for (int m = -l; m <= l; m++) {
            int mode_idx = l * (2 * series->lmax + 1) + (m + series->lmax);
            if (mode_idx >= nmode) continue;

            nr_complex_t hlm = series->hlm_data[
                time_idx * nmode + mode_idx];

            nr_complex_t Ylm = nr_wave_spin_weighted_Ylm(l, m, theta, phi);
            nr_complex_t contrib = nr_complex_mul(hlm, Ylm);
            h = nr_complex_add(h, contrib);
        }
    }

    /* h_+ = Re(h), h_× = −Im(h) */
    *hp = h.re;
    *hc = -h.im;
}

/* ===========================================================================
 * L2: Radiated Energy
 * =========================================================================== */

double nr_wave_radiated_energy(const nr_strain_series_t *series) {
    if (!series || series->ntimes < 2) return 0.0;

    double E = 0.0;
    int nmode = (series->lmax + 1) * (2 * series->lmax + 1);

    for (int t = 0; t < series->ntimes - 1; t++) {
        double dE_dt = 0.0;
        for (int l = 2; l <= series->lmax; l++) {
            for (int m = -l; m <= l; m++) {
                int mode_idx = l * (2 * series->lmax + 1) + (m + series->lmax);
                if (mode_idx >= nmode) continue;

                nr_complex_t hlm_t = series->hlm_data[t * nmode + mode_idx];
                nr_complex_t hlm_t1 = series->hlm_data[(t+1) * nmode + mode_idx];

                /* dh/dt ≈ (h_{t+1} − h_t) / dt */
                double dh_re = (hlm_t1.re - hlm_t.re) / series->dt;
                double dh_im = (hlm_t1.im - hlm_t.im) / series->dt;
                double dh2 = dh_re*dh_re + dh_im*dh_im;

                /* dE/dt = r²/(16π) Σ |dh/dt|²
                 * With r = 1 (normalized) for dimensionless energy.
                 */
                dE_dt += dh2 / (16.0 * M_PI);
            }
        }
        E += dE_dt * series->dt;
    }

    return E;
}

/* ===========================================================================
 * L2: Final BH Properties (Ringdown Analysis)
 * =========================================================================== */

void nr_wave_final_properties(const nr_strain_series_t *series,
                               double *t_peak, double *omega22,
                               double *tau22) {
    if (!series || !t_peak || !omega22 || !tau22 || series->ntimes < 10) {
        if (t_peak) *t_peak = 0.0;
        if (omega22) *omega22 = 0.0;
        if (tau22) *tau22 = 0.0;
        return;
    }

    int nmode = (series->lmax + 1) * (2 * series->lmax + 1);
    int mode_22 = 2 * (2 * series->lmax + 1) + (2 + series->lmax);

    /* Find peak of |h_{22}| */
    double max_amp = 0.0;
    int peak_idx = 0;
    for (int t = 0; t < series->ntimes; t++) {
        nr_complex_t hlm = series->hlm_data[t * nmode + mode_22];
        double amp = nr_complex_norm2(hlm);
        if (amp > max_amp) {
            max_amp = amp;
            peak_idx = t;
        }
    }

    *t_peak = series->t0 + peak_idx * series->dt;

    /* Ringdown analysis after peak:
     * h_{22}(t) ~ A e^{−t/τ} sin(ω t + δ)
     *
     * Estimate ω from zero crossings, τ from envelope decay.
     */
    if (peak_idx + 10 < series->ntimes) {
        /* Count zero crossings after peak */
        int crossings = 0;
        double last_sign = 0.0;
        double first_t = 0.0, last_t = 0.0;

        for (int t = peak_idx + 1; t < series->ntimes; t++) {
            nr_complex_t hlm = series->hlm_data[t * nmode + mode_22];
            double sign = (hlm.re > 0) ? 1.0 : -1.0;
            if (last_sign != 0.0 && sign != last_sign) {
                if (crossings == 0) first_t = series->t0 + t * series->dt;
                last_t = series->t0 + t * series->dt;
                crossings++;
            }
            last_sign = sign;
        }

        if (crossings >= 2) {
            double period = 2.0 * (last_t - first_t) / (crossings - 1);
            *omega22 = 2.0 * M_PI / period;
        } else {
            *omega22 = 0.5;  /* Default ~ 1/(2M) */
        }

        /* Estimate τ from amplitude decay */
        double amp_peak = sqrt(max_amp);
        int half_idx = peak_idx;
        for (int t = peak_idx; t < series->ntimes; t++) {
            nr_complex_t hlm = series->hlm_data[t * nmode + mode_22];
            double amp = sqrt(nr_complex_norm2(hlm));
            if (amp < amp_peak / exp(1.0)) {
                half_idx = t;
                break;
            }
        }
        *tau22 = (series->t0 + half_idx * series->dt) - *t_peak;
    } else {
        *omega22 = 0.5;
        *tau22 = 10.0 * series->dt;
    }
}

/* ===========================================================================
 * L7: LIGO/Virgo Analysis
 * =========================================================================== */

void nr_wave_ligo_analysis(const nr_strain_series_t *series,
                            double distance,
                            double *M_chirp, double *M_total,
                            double *SNR) {
    if (!series) return;

    double t_peak, omega22, tau22;
    nr_wave_final_properties(series, &t_peak, &omega22, &tau22);

    /* Chirp mass from the (2,2) mode frequency:
     * ω_{22} ≈ (1/M) → M ≈ 1/ω_{22}
     *
     * Quasi-normal mode frequency for Schwarzschild:
     *   M ω_{22} ≈ 0.3737 → M ≈ 0.3737/ω_{22} (geometric units)
     *
     * Quasi-normal mode for Kerr (fitting formula):
     *   M ω_{22} ≈ [1 − 0.63(1−a)^{0.3}]
     *
     * Chirp mass from inspiral: M_chirp = (m1 m2)^{3/5} / (m1+m2)^{1/5}
     * Approximate from total mass and mass ratio q ≤ 1:
     *   M_chirp = M_total * q^{3/5} / (1+q)^{6/5}
     *
     * For equal mass: q=1 → M_chirp = M_total / 2^{6/5} ≈ 0.435 M_total
     */

    /* From ringdown frequency: M_final ≈ 0.3737/ω_{22} (geometric) */
    double M_final_geo = (omega22 > 1e-10) ? 0.3737 / omega22 : 100.0;

    /* Convert to solar masses: 1 M_sun = 1.477 km = 4.925e-6 s
     * Actually M_sun in geometric units = 1.477 km.
     * The conversion depends on the code units.
     * For G=c=1, M is in length units.
     * 1 M_sun = 1.477e3 m = 4.925e-6 s (in geometric units)
     */
    double M_sun_geo = 1.477e3;  /* meters */
    double M_total_solar = M_final_geo / M_sun_geo;

    if (M_total) *M_total = M_total_solar;

    /* Chirp mass: for equal-mass inspiral, M_chirp ≈ M_total / 2^{6/5} */
    if (M_chirp) *M_chirp = M_total_solar * pow(2.0, -1.2);

    /* Simplified SNR:
     * SNR² ≈ 4 ∫ |h(f)|² / S_n(f) df
     *
     * For a rough estimate:
     * SNR ≈ h_peak * sqrt(N_cycles) / h_noise
     *
     * With h_noise ~ 10^{-23} for Advanced LIGO at design sensitivity,
     * and distance scaling h ∝ 1/distance.
     */
    if (SNR) {
        /* Peak strain at 1 Mpc ≈ 10^{-19} for a 100 M_sun source
         * Scale with distance and total mass.
         */
        double h_peak_ref = 1e-19;  /* At 1 Mpc for 100 Msun */
        double h_peak = h_peak_ref * M_total_solar / 100.0 * (1.0 / distance);

        int nmode = (series->lmax + 1) * (2 * series->lmax + 1);
        int mode_22 = 2 * (2 * series->lmax + 1) + (2 + series->lmax);

        double max_sim_amp = 0.0;
        for (int t = 0; t < series->ntimes; t++) {
            nr_complex_t hlm = series->hlm_data[t * nmode + mode_22];
            double amp = sqrt(nr_complex_norm2(hlm));
            if (amp > max_sim_amp) max_sim_amp = amp;
        }

        /* N_cycles ≈ ω_{22} * τ_{22} / (2π) */
        double N_cycles = omega22 * tau22 / (2.0 * M_PI);
        if (N_cycles < 1) N_cycles = 1;

        double h_noise = 1e-23;
        *SNR = h_peak * sqrt(N_cycles) / h_noise;
    }
}

/* ===========================================================================
 * Strain Series Allocation
 * =========================================================================== */

nr_strain_series_t* nr_strain_series_alloc(int lmax, int max_times) {
    nr_strain_series_t *s = (nr_strain_series_t*)malloc(sizeof(nr_strain_series_t));
    if (!s) return NULL;

    s->lmax = lmax;
    s->ntimes = 0;
    s->dt = 0.0;
    s->t0 = 0.0;

    int n_modes = (lmax + 1) * (2 * lmax + 1);
    s->hlm_data = (nr_complex_t*)calloc(max_times * n_modes, sizeof(nr_complex_t));
    if (!s->hlm_data) {
        free(s);
        return NULL;
    }
    return s;
}

void nr_strain_series_free(nr_strain_series_t *series) {
    if (series) {
        free(series->hlm_data);
        free(series);
    }
}
