/**
 * @file nr_adm.c
 * @brief ADM 3+1 decomposition: metric, curvature, constraints, evolution.
 *
 * Implements the core mathematical operations of the ADM formalism:
 * Christoffel symbols, Ricci tensor, Hamiltonian and momentum constraints,
 * and the right-hand side of the ADM evolution equations.
 *
 * References:
 *   - Arnowitt, Deser, Misner (1962)
 *   - Wald (1984), Appendix E
 *   - Baumgarte & Shapiro (2010), Ch. 2
 *   - Gourgoulhon (2012), "3+1 Formalism in General Relativity"
 *
 * Knowledge: L1 (ADM variables), L2 (foliation), L3 (3-Riemannian geometry),
 *            L4 (Einstein constraints, evolution equations)
 */

#include "nr_adm.h"
#include "nr_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===========================================================================
 * ADM State Management
 * =========================================================================== */

nr_adm_state_t* nr_adm_alloc(int nx, int ny, int nz, int ng) {
    nr_adm_state_t *state = (nr_adm_state_t*)malloc(sizeof(nr_adm_state_t));
    if (!state) return NULL;

    state->gamma = nr_sym_tensor3_alloc(nx, ny, nz, ng);
    state->K     = nr_sym_tensor3_alloc(nx, ny, nz, ng);
    state->alpha = nr_gf_alloc(nx, ny, nz, ng);
    state->beta  = nr_vector3_alloc(nx, ny, nz, ng);

    if (!state->gamma || !state->K || !state->alpha || !state->beta) {
        nr_adm_free(state);
        return NULL;
    }
    return state;
}

void nr_adm_free(nr_adm_state_t *state) {
    if (state) {
        nr_sym_tensor3_free(state->gamma);
        nr_sym_tensor3_free(state->K);
        nr_gf_free(state->alpha);
        nr_vector3_free(state->beta);
        free(state);
    }
}

void nr_adm_set_minkowski(nr_adm_state_t *state) {
    if (!state) return;
    /* γ_{ij} = δ_{ij} */
    nr_gf_set_all(state->gamma->comp[0], 1.0);  /* γ_xx */
    nr_gf_set_all(state->gamma->comp[1], 0.0);  /* γ_xy */
    nr_gf_set_all(state->gamma->comp[2], 0.0);  /* γ_xz */
    nr_gf_set_all(state->gamma->comp[3], 1.0);  /* γ_yy */
    nr_gf_set_all(state->gamma->comp[4], 0.0);  /* γ_yz */
    nr_gf_set_all(state->gamma->comp[5], 1.0);  /* γ_zz */

    /* K_{ij} = 0 */
    for (int c = 0; c < 6; c++) nr_gf_set_all(state->K->comp[c], 0.0);

    /* α = 1, β^i = 0 */
    nr_gf_set_all(state->alpha, 1.0);
    for (int c = 0; c < 3; c++) nr_gf_set_all(state->beta->comp[c], 0.0);
}

/* ===========================================================================
 * L3: 3-Metric Inversion
 * =========================================================================== */

double nr_adm_inverse_metric(const double gamma_ij[6],
                             double inv_gamma[3][3]) {
    /* Symmetric 3x3 matrix inversion.
     * gamma_ij = {xx, xy, xz, yy, yz, zz}
     *
     * Determinant via cofactor expansion:
     * det = gxx*gyy*gzz + 2*gxy*gyz*gxz - gxx*gyz^2 - gyy*gxz^2 - gzz*gxy^2
     */
    double gxx = gamma_ij[0], gxy = gamma_ij[1], gxz = gamma_ij[2];
    double gyy = gamma_ij[3], gyz = gamma_ij[4], gzz = gamma_ij[5];

    double det = gxx * gyy * gzz
               + 2.0 * gxy * gyz * gxz
               - gxx * gyz * gyz
               - gyy * gxz * gxz
               - gzz * gxy * gxy;

    if (fabs(det) < 1e-30) {
        /* Degenerate metric — return identity inverse */
        inv_gamma[0][0] = inv_gamma[1][1] = inv_gamma[2][2] = 1.0;
        inv_gamma[0][1] = inv_gamma[0][2] = 0.0;
        inv_gamma[1][0] = inv_gamma[1][2] = 0.0;
        inv_gamma[2][0] = inv_gamma[2][1] = 0.0;
        return det;
    }

    double inv_det = 1.0 / det;

    /* Cofactors of the symmetric matrix:
     * C^{xx} = gyy*gzz - gyz^2
     * C^{xy} = gyz*gxz - gxy*gzz
     * C^{xz} = gxy*gyz - gyy*gxz
     * C^{yy} = gxx*gzz - gxz^2
     * C^{yz} = gxy*gxz - gxx*gyz
     * C^{zz} = gxx*gyy - gxy^2
     */
    inv_gamma[0][0] = (gyy * gzz - gyz * gyz) * inv_det;
    inv_gamma[0][1] = (gyz * gxz - gxy * gzz) * inv_det;
    inv_gamma[0][2] = (gxy * gyz - gyy * gxz) * inv_det;
    inv_gamma[1][0] = inv_gamma[0][1];
    inv_gamma[1][1] = (gxx * gzz - gxz * gxz) * inv_det;
    inv_gamma[1][2] = (gxy * gxz - gxx * gyz) * inv_det;
    inv_gamma[2][0] = inv_gamma[0][2];
    inv_gamma[2][1] = inv_gamma[1][2];
    inv_gamma[2][2] = (gxx * gyy - gxy * gxy) * inv_det;

    return det;
}

/* ===========================================================================
 * L3: Christoffel Symbols & Ricci Tensor
 *
 * Γ^i_{jk} = (1/2) γ^{il} (∂_j γ_{lk} + ∂_k γ_{lj} − ∂_l γ_{jk})
 *
 * R_{ij} = ∂_k Γ^k_{ij} − ∂_i Γ^k_{kj} + Γ^k_{kl} Γ^l_{ij} − Γ^k_{il} Γ^l_{kj}
 * =========================================================================== */

/* Helper: get metric component from symmetric tensor at grid point */
static double sym_comp(const nr_sym_tensor3_t *t, int c,
                       int i, int j, int k) {
    return t->comp[c]->data[i * t->comp[c]->stride_x
                            + j * t->comp[c]->stride_y + k];
}

void nr_adm_christoffel(const nr_sym_tensor3_t *gamma,
                         int i, int j, int k, double dx,
                         double Gamma[3][3][3]) {
    /* Get the 3-metric at this point */
    double g[6];
    for (int c = 0; c < 6; c++) {
        g[c] = gamma->comp[c]->data[i * gamma->comp[c]->stride_x
                                   + j * gamma->comp[c]->stride_y + k];
    }

    /* Inverse metric */
    double inv[3][3];
    nr_adm_inverse_metric(g, inv);

    /* Compute ∂_a g_{bc} — 18 components (6 metric comps × 3 directions)
     * Use 4th-order FD
     */
    double d_g[6][3];  /* d_g[comp][dir] */
    int order = 4;
    for (int c = 0; c < 6; c++) {
        d_g[c][0] = nr_fd_deriv1(gamma->comp[c], i, j, k, dx, 0, order);
        d_g[c][1] = nr_fd_deriv1(gamma->comp[c], i, j, k, dx, 1, order);
        d_g[c][2] = nr_fd_deriv1(gamma->comp[c], i, j, k, dx, 2, order);
    }

    /* Map: metric component index → (a,b) pair
     * 0: xx (0,0), 1: xy (0,1), 2: xz (0,2),
     * 3: yy (1,1), 4: yz (1,2), 5: zz (2,2)
     */
    static const int map_a[6] = {0, 0, 0, 1, 1, 2};
    static const int map_b[6] = {0, 1, 2, 1, 2, 2};

    /* Initialize Gamma to zero */
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++)
                Gamma[a][b][c] = 0.0;

    /* Γ^i_{jk} = (1/2) Σ_l γ^{il} (∂_j g_{lk} + ∂_k g_{lj} − ∂_l g_{jk}) */
    for (int i_idx = 0; i_idx < 3; i_idx++) {
        for (int j_idx = 0; j_idx < 3; j_idx++) {
            for (int k_idx = 0; k_idx < 3; k_idx++) {
                double sum = 0.0;
                for (int l = 0; l < 3; l++) {
                    double dg_lkj = 0.0, dg_ljk = 0.0, dg_jkl = 0.0;
                    /* Find ∂_j g_{lk} */
                    for (int m = 0; m < 6; m++) {
                        if (map_a[m] == l && map_b[m] == k_idx) dg_lkj = d_g[m][j_idx];
                        else if (map_a[m] == k_idx && map_b[m] == l) dg_lkj = d_g[m][j_idx];
                        if (map_a[m] == l && map_b[m] == j_idx) dg_ljk = d_g[m][k_idx];
                        else if (map_a[m] == j_idx && map_b[m] == l) dg_ljk = d_g[m][k_idx];
                        if (map_a[m] == j_idx && map_b[m] == k_idx) dg_jkl = d_g[m][l];
                        else if (map_a[m] == k_idx && map_b[m] == j_idx) dg_jkl = d_g[m][l];
                    }
                    double term = dg_lkj + dg_ljk - dg_jkl;
                    sum += inv[i_idx][l] * term;
                }
                Gamma[i_idx][j_idx][k_idx] = 0.5 * sum;
            }
        }
    }
}

void nr_adm_ricci(const nr_sym_tensor3_t *gamma,
                  int i, int j, int k, double dx,
                  double Rij[3][3]) {
    int order = 4;
    (void)(gamma->comp[0]->ng);  /* ng accessed for consistency */

    /* Compute Gamma^i_{jk} at the center point */
    double G0[3][3][3];
    nr_adm_christoffel(gamma, i, j, k, dx, G0);

    /* We need derivatives of Gamma^i_{jk}.
     * For each (a,b,c), compute the 3-component "grid function" of Gamma^a_{bc}
     * evaluated at neighboring points, then FD.
     *
     * Since we don't have precomputed Gamma, we compute at stencil points.
     * For R_{ij} = ∂_k Γ^k_{ij} − ∂_i Γ^k_{kj} + Γ^k_{kl}Γ^l_{ij} − Γ^k_{il}Γ^l_{kj}
     *
     * We'll use a simplified approach: compute Gamma at stencil points
     * and then FD the terms.
     *
     * For a 3D grid with N ~ 10-100, computing Christoffel symbols at
     * stencil points is O(stencil_points) = O(1).
     */

    /* We need ∂_k Γ^k_{ij} and ∂_i Γ^k_{kj}.
     *
     * ∂_k Γ^k_{ij} = Σ_{k=0}^{2} ∂_k (Gamma[k][i][j])
     *
     * Each term ∂_k Gamma[k][i][j] is a directional derivative.
     * We compute it via FD on the "Gamma function".
     */

    /* Precompute Gamma at all stencil points (±1 in each of 3 directions) */
    (void)order;  /* Use default order from FD calls */
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            Rij[a][b] = 0.0;

    double G_p[3][3][3][3];  /* G_p[dir][a][b][c] */
    double G_n[3][3][3][3];
    int di_dirs[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

    for (int d = 0; d < 3; d++) {
        int di = di_dirs[d][0], dj = di_dirs[d][1], dk = di_dirs[d][2];
        nr_adm_christoffel(gamma, i+di, j+dj, k+dk, dx, G_p[d]);
        nr_adm_christoffel(gamma, i-di, j-dj, k-dk, dx, G_n[d]);
    }

    /* Term 1: ∂_k Γ^k_{ij} — derivative of Gamma in its upper index direction
     * ∂_0 Γ^0_{ij} + ∂_1 Γ^1_{ij} + ∂_2 Γ^2_{ij}
     */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int d = 0; d < 3; d++) {
                double fp = G_p[d][d][a][b];
                double fm = G_n[d][d][a][b];
                Rij[a][b] += (fp - fm) / (2.0 * dx);  /* ∂_d Γ^d_{ab} */
            }
        }
    }

    /* Term 2: −∂_i Γ^k_{kj} = −∂_a Γ^k_{kb}
     * For each (a,b): compute Σ_k ∂_a Γ^k_{kb}
     */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            /* Σ_k ∂_a Γ^k_{kb} */
            for (int kk = 0; kk < 3; kk++) {
                double fp = G_p[a][kk][kk][b];
                double fm = G_n[a][kk][kk][b];
                Rij[a][b] -= (fp - fm) / (2.0 * dx);
            }
        }
    }

    /* Term 3: Γ^k_{kl} Γ^l_{ij} */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            /* Σ_{k,l} Γ^k_{kl} · Γ^l_{ab} */
            double Gk_kl = 0.0;
            for (int kk = 0; kk < 3; kk++) {
                for (int l = 0; l < 3; l++) {
                    Gk_kl += G0[kk][kk][l];
                }
            }
            /* Actually: Σ_{k,l} Γ^k_{kl} Γ^l_{ij} */
            for (int kk = 0; kk < 3; kk++) {
                for (int l = 0; l < 3; l++) {
                    Rij[a][b] += G0[kk][kk][l] * G0[l][a][b];
                }
            }
        }
    }

    /* Term 4: −Γ^k_{il} Γ^l_{kj} = −Σ_{k,l} Γ^k_{ia}·Γ^l_{kb} must use correct indices
     * Correct: −Γ^k_{i l} Γ^l_{k j}
     */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int kk = 0; kk < 3; kk++) {
                for (int l = 0; l < 3; l++) {
                    Rij[a][b] -= G0[kk][a][l] * G0[l][kk][b];
                }
            }
        }
    }
}

double nr_adm_ricci_scalar(const double Rij[3][3],
                           const double inv_gamma[3][3]) {
    double R = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R += inv_gamma[i][j] * Rij[i][j];
        }
    }
    return R;
}

/* ===========================================================================
 * L4: Constraints
 * =========================================================================== */

double nr_adm_hamiltonian_constraint(const nr_sym_tensor3_t *gamma,
                                      const nr_sym_tensor3_t *K,
                                      int i, int j, int k, double dx) {
    /* H = R + (tr K)^2 − K_{ij} K^{ij} */

    /* 1. 3-Ricci scalar R */
    double Rij[3][3];
    nr_adm_ricci(gamma, i, j, k, dx, Rij);

    /* 2. Inverse metric */
    double g_ij[6];
    for (int c = 0; c < 6; c++) {
        g_ij[c] = sym_comp(gamma, c, i, j, k);
    }
    double inv_g[3][3];
    nr_adm_inverse_metric(g_ij, inv_g);

    double R = nr_adm_ricci_scalar(Rij, inv_g);

    /* 3. Trace of extrinsic curvature: K = γ^{ij} K_{ij} */
    double trK = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            /* K_{ab} — get from symmetric tensor */
            double Kab = 0.0;
            /* Map (a,b) → component index */
            int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
            int c = map[a][b];
            Kab = sym_comp(K, c, i, j, k);
            trK += inv_g[a][b] * Kab;
        }
    }

    /* 4. K_{ij} K^{ij} = K_{ab} K_{cd} γ^{ac} γ^{bd} */
    double K2 = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                for (int d = 0; d < 3; d++) {
                    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
                    int ci = map[a][b], cj = map[c][d];
                    double Kab = sym_comp(K, ci, i, j, k);
                    double Kcd = sym_comp(K, cj, i, j, k);
                    K2 += inv_g[a][c] * inv_g[b][d] * Kab * Kcd;
                }
            }
        }
    }

    return R + trK * trK - K2;
}

void nr_adm_momentum_constraint(const nr_sym_tensor3_t *gamma,
                                 const nr_sym_tensor3_t *K,
                                 int i, int j, int k, double dx,
                                 double Mi[3]) {
    /* M_i = D_j K^j_i − D_i K
     *
     * D_j K^j_i = γ^{jl} D_l K_{ji}
     *           = γ^{jl} (∂_l K_{ji} − Γ^m_{lj} K_{mi} − Γ^m_{li} K_{jm})
     *
     * D_i K = ∂_i K
     */

    int order = 4;

    /* Get Gamma at this point */
    double G[3][3][3];
    nr_adm_christoffel(gamma, i, j, k, dx, G);

    /* Inverse metric */
    double g_ij[6];
    for (int c = 0; c < 6; c++) {
        g_ij[c] = sym_comp(gamma, c, i, j, k);
    }
    double inv_g[3][3];
    nr_adm_inverse_metric(g_ij, inv_g);

    /* Compute K components at this point */
    double Kij[3][3];
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            Kij[a][b] = sym_comp(K, map[a][b], i, j, k);

    /* Trace of K: K = γ^{ij} K_{ij} */
    double trK = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            trK += inv_g[a][b] * Kij[a][b];

    /* Compute partial derivatives of K and K_{ij} */
    double d_K[3];  /* ∂_i K */
    double d_Kij[6][3];  /* ∂_k K_{ij} for each component */

    for (int dir = 0; dir < 3; dir++) {
        /* ∂_i K: Σ_{a,b} ∂_i (γ^{ab} K_{ab})
         * Simplified: compute from FG stencil of K field
         * We compute ∂_i K = ∂_i(γ^{ab} K_{ab})
         *
         * For simplicity: use an anonymous gf to hold K at neighboring
         * points and compute FD derivative.
         */
        /* Approximate ∂_i K via finite diff on trK at stencil:
         * K(i+1) computed explicitly */
        double trK_p = 0.0, trK_m = 0.0;
        /* K at i+1 in dir */
        {
            double g_p[6], inv_p[3][3], K_p[3][3];
            int di=0,dj=0,dk=0;
            if(dir==0) di=1; else if(dir==1) dj=1; else dk=1;
            for (int c=0;c<6;c++) g_p[c]=sym_comp(gamma,c,i+di,j+dj,k+dk);
            nr_adm_inverse_metric(g_p, inv_p);
            for (int a=0;a<3;a++) for(int b=0;b<3;b++)
                K_p[a][b]=sym_comp(K,map[a][b],i+di,j+dj,k+dk);
            for (int a=0;a<3;a++) for(int b=0;b<3;b++)
                trK_p += inv_p[a][b] * K_p[a][b];
        }
        {
            double g_m[6], inv_m[3][3], K_m[3][3];
            int di=0,dj=0,dk=0;
            if(dir==0) di=-1; else if(dir==1) dj=-1; else dk=-1;
            for (int c=0;c<6;c++) g_m[c]=sym_comp(gamma,c,i+di,j+dj,k+dk);
            nr_adm_inverse_metric(g_m, inv_m);
            for (int a=0;a<3;a++) for(int b=0;b<3;b++)
                K_m[a][b]=sym_comp(K,map[a][b],i+di,j+dj,k+dk);
            for (int a=0;a<3;a++) for(int b=0;b<3;b++)
                trK_m += inv_m[a][b] * K_m[a][b];
        }
        d_K[dir] = (trK_p - trK_m) / (2.0 * dx);

        /* ∂_k K_{ab} for each component */
        for (int c = 0; c < 6; c++) {
            d_Kij[c][dir] = nr_fd_deriv1(K->comp[c], i, j, k, dx, dir, order);
        }
    }

    /* Compute M_i = D_j K^j_i − ∂_i K */
    for (int mi = 0; mi < 3; mi++) {
        double Dj_Kji = 0.0;

        /* D_j K^j_i = γ^{jl} (∂_l K_{ji} − Γ^m_{lj} K_{mi} − Γ^m_{li} K_{jm}) */
        for (int jj = 0; jj < 3; jj++) {
            for (int l = 0; l < 3; l++) {
                double inv_jl = inv_g[jj][l];

                /* ∂_l K_{ji} */
                double dl_Kji = 0.0;
                int c_ji = map[jj][mi];
                dl_Kji = d_Kij[c_ji][l];

                /* − Γ^m_{lj} K_{mi} */
                double Gamma_K1 = 0.0;
                for (int m = 0; m < 3; m++) {
                    Gamma_K1 += G[m][l][jj] * Kij[m][mi];
                }

                /* − Γ^m_{li} K_{jm} */
                double Gamma_K2 = 0.0;
                for (int m = 0; m < 3; m++) {
                    Gamma_K2 += G[m][l][mi] * Kij[jj][m];
                }

                Dj_Kji += inv_jl * (dl_Kji - Gamma_K1 - Gamma_K2);
            }
        }

        Mi[mi] = Dj_Kji - d_K[mi];
    }
}

void nr_adm_constraint_norms(const nr_adm_state_t *state, double dx,
                              double *norm_H, double *norm_M) {
    if (!state) { *norm_H = *norm_M = 0.0; return; }
    int ng = state->gamma->comp[0]->ng;
    int nx = state->gamma->comp[0]->nx;
    int ny = state->gamma->comp[0]->ny;
    int nz = state->gamma->comp[0]->nz;

    double sum_H = 0.0, sum_M = 0.0;
    int count = 0;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                double H = nr_adm_hamiltonian_constraint(
                    state->gamma, state->K, i, j, k, dx);
                sum_H += H * H;

                double Mi[3];
                nr_adm_momentum_constraint(
                    state->gamma, state->K, i, j, k, dx, Mi);
                sum_M += Mi[0]*Mi[0] + Mi[1]*Mi[1] + Mi[2]*Mi[2];

                count++;
            }
        }
    }

    *norm_H = (count > 0) ? sqrt(sum_H / count) : 0.0;
    *norm_M = (count > 0) ? sqrt(sum_M / count) : 0.0;
}

/* ===========================================================================
 * L4: ADM Evolution RHS
 * =========================================================================== */

void nr_adm_rhs_point(const nr_adm_state_t *state,
                       int i, int j, int k, double dx,
                       double dt_gamma[6], double dt_K[6]) {
    if (!state) return;

    int order = 4;

    /* Get metric and curvature at this point */
    double g[6], Kij[6];
    for (int c = 0; c < 6; c++) {
        g[c] = sym_comp(state->gamma, c, i, j, k);
        Kij[c] = sym_comp(state->K, c, i, j, k);
    }
    double alpha = state->alpha->data[i * state->alpha->stride_x
                                     + j * state->alpha->stride_y + k];

    /* Inverse metric */
    double inv_g[3][3];
    nr_adm_inverse_metric(g, inv_g);

    /* 3-Ricci tensor */
    double Rij[3][3];
    nr_adm_ricci(state->gamma, i, j, k, dx, Rij);

    /* Christoffel symbols */
    double G[3][3][3];
    nr_adm_christoffel(state->gamma, i, j, k, dx, G);

    /* Trace of K */
    double trK = 0.0;
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            trK += inv_g[a][b] * Kij[map[a][b]];

    /* ∂_i α, ∂_i ∂_j α */
    double da[3];
    for (int dir = 0; dir < 3; dir++)
        da[dir] = nr_fd_deriv1(state->alpha, i, j, k, dx, dir, order);

    /* ∂_k g_{ij} */
    double dg[6][3];
    for (int c = 0; c < 6; c++)
        for (int dir = 0; dir < 3; dir++)
            dg[c][dir] = nr_fd_deriv1(state->gamma->comp[c], i, j, k, dx, dir, order);

    /* ∂_k β^i */
    double db[3][3];
    for (int comp = 0; comp < 3; comp++)
        for (int dir = 0; dir < 3; dir++)
            db[comp][dir] = nr_fd_deriv1(state->beta->comp[comp], i, j, k, dx, dir, order);

    /* β^i */
    double beta[3];
    for (int c = 0; c < 3; c++)
        beta[c] = state->beta->comp[c]->data[i * state->beta->comp[c]->stride_x
                                             + j * state->beta->comp[c]->stride_y + k];

    /* ========= ∂_t γ_{ij} = −2α K_{ij} + D_i β_j + D_j β_i ========= */
    /* D_i β_j = ∂_i β_j − Γ^k_{ij} β_k */
    /* where β_j = γ_{jk} β^k */
    double beta_lower[3];
    for (int a = 0; a < 3; a++) {
        beta_lower[a] = 0.0;
        for (int b = 0; b < 3; b++) {
            int c = map[a][b];
            beta_lower[a] += g[c] * beta[b];
        }
    }

    /* D_a β_b */
    double Dbeta[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double d_a_beta_b = 0.0;
            /* ∂_a γ_{bc} β^c + γ_{bc} ∂_a β^c */
            for (int c = 0; c < 3; c++) {
                int ci = map[b][c];
                d_a_beta_b += dg[ci][a] * beta[c] + g[ci] * db[c][a];
            }
            double Gamma_term = 0.0;
            for (int kk = 0; kk < 3; kk++)
                Gamma_term += G[kk][a][b] * beta_lower[kk];
            Dbeta[a][b] = d_a_beta_b - Gamma_term;
        }
    }

    for (int c = 0; c < 6; c++) {
        int aa[6] = {0,0,0,1,1,2};
        int bb[6] = {0,1,2,1,2,2};
        int a = aa[c], b = bb[c];
        dt_gamma[c] = -2.0 * alpha * Kij[c] + Dbeta[a][b] + Dbeta[b][a];
    }

    /* ========= ∂_t K_{ij} = ... ========= */
    /* ∂_i ∂_j α */
    double dda[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            if (a == b) {
                dda[a][b] = nr_fd_deriv2(state->alpha, i, j, k, dx, a, order);
            } else {
                dda[a][b] = nr_fd_deriv2_mixed(state->alpha, i, j, k, dx, dx, a, b, order);
            }
        }
    }

    /* Covariant derivative D_i D_j α = ∂_i ∂_j α − Γ^k_{ij} ∂_k α */
    double DDalpha[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            DDalpha[a][b] = dda[a][b];
            for (int kk = 0; kk < 3; kk++)
                DDalpha[a][b] -= G[kk][a][b] * da[kk];
        }
    }

    /* β^k ∂_k K_{ij} */
    double beta_dK[6];
    for (int c = 0; c < 6; c++) {
        beta_dK[c] = 0.0;
        for (int dir = 0; dir < 3; dir++) {
            beta_dK[c] += beta[dir] * nr_fd_deriv1(state->K->comp[c], i, j, k, dx, dir, order);
        }
    }

    for (int c = 0; c < 6; c++) {
        int aa[6] = {0,0,0,1,1,2};
        int bb[6] = {0,1,2,1,2,2};
        int a = aa[c], b = bb[c];

        /* − D_a D_b α */
        double term1 = -DDalpha[a][b];

        /* α (R_{ab} + K K_{ab} − 2 Σ_{k} K_{ak} K^k_b) */
        double KaiKkb = 0.0;
        for (int kk = 0; kk < 3; kk++) {
            int c_ak = map[a][kk];
            /* K^k_b = γ^{kl} K_{lb} */
            double Kkb = 0.0;
            for (int l = 0; l < 3; l++) {
                int c_lb = map[l][b];
                Kkb += inv_g[kk][l] * Kij[c_lb];
            }
            KaiKkb += Kij[c_ak] * Kkb;
        }
        double term2 = alpha * (Rij[a][b] + trK * Kij[c] - 2.0 * KaiKkb);

        /* β^k ∂_k K_{ij} */
        double term3 = beta_dK[c];

        /* K_{ik} ∂_j β^k + K_{kj} ∂_i β^k */
        double term4 = 0.0;
        for (int kk = 0; kk < 3; kk++) {
            int c_ak = map[a][kk], c_kb = map[kk][b];
            term4 += Kij[c_ak] * db[kk][b] + Kij[c_kb] * db[kk][a];
        }

        dt_K[c] = term1 + term2 + term3 + term4;
    }
}

void nr_adm_rhs(const nr_adm_state_t *state, double dx,
                nr_adm_state_t *dt_state) {
    if (!state || !dt_state) return;
    int ng = state->gamma->comp[0]->ng;
    int nx = state->gamma->comp[0]->nx;
    int ny = state->gamma->comp[0]->ny;
    int nz = state->gamma->comp[0]->nz;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                double dt_g[6], dt_K[6];
                nr_adm_rhs_point(state, i, j, k, dx, dt_g, dt_K);
                for (int c = 0; c < 6; c++) {
                    dt_state->gamma->comp[c]->data[
                        i*dt_state->gamma->comp[c]->stride_x
                        + j*dt_state->gamma->comp[c]->stride_y + k] = dt_g[c];
                    dt_state->K->comp[c]->data[
                        i*dt_state->K->comp[c]->stride_x
                        + j*dt_state->K->comp[c]->stride_y + k] = dt_K[c];
                }
                /* Lapse and shift RHS set to 0 for now (handled by gauge) */
            }
        }
    }
}

/* ===========================================================================
 * L2: ADM Mass
 * =========================================================================== */

double nr_adm_mass(const nr_gf_t *psi, const nr_grid_t *grid, double r_ext) {
    /* M_ADM = −1/(2π) ∮_S_∞ ∂_r ψ dA
     *
     * For isotropic Schwarzschild: ψ = 1 + M/(2r) → M = 2r(ψ−1)
     *
     * We estimate M from the far-field average of 2r(ψ−1),
     * weighted by distance from the origin. This is more robust
     * on coarse grids than the surface integral approach.
     */
    if (!psi || !grid) return 0.0;

    int ng = grid->ng;
    double sum = 0.0;
    double weight_sum = 0.0;

    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];
                double r = sqrt(x*x + y*y + z*z);

                /* Use points near the extraction radius */
                double dr = fabs(r - r_ext);
                double weight = exp(-dr*dr/(grid->dx*grid->dx));
                if (r > 1.0 && weight > 1e-3) {
                    double psi_val = psi->data[i*psi->stride_x
                                               + j*psi->stride_y + k];
                    double M_est = 2.0 * r * (psi_val - 1.0);
                    sum += M_est * weight;
                    weight_sum += weight;
                }
            }
        }
    }

    if (weight_sum < 1e-15) return 0.0;
    return sum / weight_sum;
}

double nr_adm_komar_mass(const nr_gf_t *alpha,
                          const nr_sym_tensor3_t *gamma,
                          const nr_grid_t *grid, double r_ext) {
    /* M_Komar = (1/(4π)) ∮_S D^i α dS_i
     *
     * For Schwarzschild in isotropic coords: α = (1−M/(2r))/(1+M/(2r))
     * M_Komar = M.
     */
    if (!alpha || !gamma || !grid) return 0.0;

    int ng = grid->ng;
    double sum = 0.0;
    int count = 0;

    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];
                double r = sqrt(x*x + y*y + z*z);

                if (fabs(r - r_ext) < grid->dx * 0.75 && r > 0) {
                    /* D^i α = γ^{ij} ∂_j α */
                    double da_dx = nr_fd_deriv1(alpha, i, j, k, grid->dx, 0, 2);
                    double da_dy = nr_fd_deriv1(alpha, i, j, k, grid->dx, 1, 2);
                    double da_dz = nr_fd_deriv1(alpha, i, j, k, grid->dx, 2, 2);

                    /* Get metric for γ^{ij} (use identity for conformally flat) */
                    double g[6];
                    for (int c = 0; c < 6; c++)
                        g[c] = gamma->comp[c]->data[
                            i*gamma->comp[c]->stride_x
                            + j*gamma->comp[c]->stride_y + k];
                    double inv_g[3][3];
                    nr_adm_inverse_metric(g, inv_g);

                    double D_alpha[3];
                    D_alpha[0] = inv_g[0][0]*da_dx + inv_g[0][1]*da_dy + inv_g[0][2]*da_dz;
                    D_alpha[1] = inv_g[1][0]*da_dx + inv_g[1][1]*da_dy + inv_g[1][2]*da_dz;
                    D_alpha[2] = inv_g[2][0]*da_dx + inv_g[2][1]*da_dy + inv_g[2][2]*da_dz;

                    /* Radial unit vector: n_i = (x/r, y/r, z/r) */
                    double nx = x/r, ny = y/r, nz = z/r;

                    /* D^i α n_i * dA,  dA = r^2 dΩ */
                    double Ddotn = D_alpha[0]*nx + D_alpha[1]*ny + D_alpha[2]*nz;
                    double dOmega = (grid->dx * grid->dy) / (r * r);
                    sum += Ddotn * r * r * dOmega;
                    count++;
                }
            }
        }
    }

    if (count == 0) return 0.0;
    return sum / (4.0 * M_PI);
}
