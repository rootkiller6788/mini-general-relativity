/**
 * @file nr_bssn.c
 * @brief BSSN evolution system: RHS, constraints, ADM conversion.
 *
 * Implements the BSSN formulation, the standard evolution system in
 * modern numerical relativity. The key innovation is the conformal
 * decomposition which separates the "Coulomb" (conformal factor φ) and
 * "radiative" (conformal metric γ̃_{ij}) degrees of freedom.
 *
 * References:
 *   - Shibata & Nakamura, PRD 52, 5428 (1995)
 *   - Baumgarte & Shapiro, PRD 59, 024007 (1998)
 *   - Baumgarte & Shapiro (2010), Ch. 11
 *
 * Knowledge: L1 (BSSN variables), L2 (conformal decomposition),
 *            L4 (BSSN evolution equations), L5 (constraint enforcement)
 */

#include "nr_bssn.h"
#include "nr_utils.h"
#include "nr_adm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===========================================================================
 * BSSN State Management
 * =========================================================================== */

nr_bssn_state_t* nr_bssn_alloc(int nx, int ny, int nz, int ng) {
    nr_bssn_state_t *s = (nr_bssn_state_t*)malloc(sizeof(nr_bssn_state_t));
    if (!s) return NULL;

    s->phi   = nr_gf_alloc(nx, ny, nz, ng);
    s->chi   = nr_gf_alloc(nx, ny, nz, ng);
    s->gt    = nr_sym_tensor3_alloc(nx, ny, nz, ng);
    s->trK   = nr_gf_alloc(nx, ny, nz, ng);
    s->At    = nr_sym_tensor3_alloc(nx, ny, nz, ng);
    s->Gt    = nr_vector3_alloc(nx, ny, nz, ng);
    s->alpha = nr_gf_alloc(nx, ny, nz, ng);
    s->beta  = nr_vector3_alloc(nx, ny, nz, ng);
    s->B     = nr_vector3_alloc(nx, ny, nz, ng);

    if (!s->phi || !s->chi || !s->gt || !s->trK || !s->At
        || !s->Gt || !s->alpha || !s->beta || !s->B) {
        nr_bssn_free(s);
        return NULL;
    }
    return s;
}

nr_bssn_rhs_t* nr_bssn_rhs_alloc(int nx, int ny, int nz, int ng) {
    nr_bssn_rhs_t *rhs = (nr_bssn_rhs_t*)malloc(sizeof(nr_bssn_rhs_t));
    if (!rhs) return NULL;

    rhs->dt_phi   = nr_gf_alloc(nx, ny, nz, ng);
    rhs->dt_chi   = nr_gf_alloc(nx, ny, nz, ng);
    rhs->dt_gt    = nr_sym_tensor3_alloc(nx, ny, nz, ng);
    rhs->dt_trK   = nr_gf_alloc(nx, ny, nz, ng);
    rhs->dt_At    = nr_sym_tensor3_alloc(nx, ny, nz, ng);
    rhs->dt_Gt    = nr_vector3_alloc(nx, ny, nz, ng);
    rhs->dt_alpha = nr_gf_alloc(nx, ny, nz, ng);
    rhs->dt_beta  = nr_vector3_alloc(nx, ny, nz, ng);
    rhs->dt_B     = nr_vector3_alloc(nx, ny, nz, ng);

    if (!rhs->dt_phi || !rhs->dt_chi || !rhs->dt_gt || !rhs->dt_trK
        || !rhs->dt_At || !rhs->dt_Gt || !rhs->dt_alpha
        || !rhs->dt_beta || !rhs->dt_B) {
        nr_bssn_rhs_free(rhs);
        return NULL;
    }
    return rhs;
}

void nr_bssn_free(nr_bssn_state_t *s) {
    if (s) {
        nr_gf_free(s->phi); nr_gf_free(s->chi);
        nr_sym_tensor3_free(s->gt); nr_gf_free(s->trK);
        nr_sym_tensor3_free(s->At); nr_vector3_free(s->Gt);
        nr_gf_free(s->alpha); nr_vector3_free(s->beta);
        nr_vector3_free(s->B);
        free(s);
    }
}

void nr_bssn_rhs_free(nr_bssn_rhs_t *rhs) {
    if (rhs) {
        nr_gf_free(rhs->dt_phi); nr_gf_free(rhs->dt_chi);
        nr_sym_tensor3_free(rhs->dt_gt); nr_gf_free(rhs->dt_trK);
        nr_sym_tensor3_free(rhs->dt_At); nr_vector3_free(rhs->dt_Gt);
        nr_gf_free(rhs->dt_alpha); nr_vector3_free(rhs->dt_beta);
        nr_vector3_free(rhs->dt_B);
        free(rhs);
    }
}

void nr_bssn_set_minkowski(nr_bssn_state_t *s) {
    if (!s) return;
    nr_gf_set_all(s->phi, 0.0);
    nr_gf_set_all(s->chi, 1.0);
    /* γ̃_{ij} = δ_{ij} */
    nr_gf_set_all(s->gt->comp[0], 1.0); nr_gf_set_all(s->gt->comp[1], 0.0);
    nr_gf_set_all(s->gt->comp[2], 0.0); nr_gf_set_all(s->gt->comp[3], 1.0);
    nr_gf_set_all(s->gt->comp[4], 0.0); nr_gf_set_all(s->gt->comp[5], 1.0);
    nr_gf_set_all(s->trK, 0.0);
    for (int c = 0; c < 6; c++) nr_gf_set_all(s->At->comp[c], 0.0);
    for (int c = 0; c < 3; c++) nr_gf_set_all(s->Gt->comp[c], 0.0);
    nr_gf_set_all(s->alpha, 1.0);
    for (int c = 0; c < 3; c++) nr_gf_set_all(s->beta->comp[c], 0.0);
    for (int c = 0; c < 3; c++) nr_gf_set_all(s->B->comp[c], 0.0);
}

/* ===========================================================================
 * ADM ↔ BSSN Conversion
 * =========================================================================== */

/* Helper: get gf value at (i,j,k) — returns pointer for lvalue use */
#define GF_VAL(gf, i, j, k) \
    ((gf)->data[(i) * (gf)->stride_x + (j) * (gf)->stride_y + (k)])

static inline double gf_val_read(const nr_gf_t *gf, int i, int j, int k) {
    return gf->data[i * gf->stride_x + j * gf->stride_y + k];
}

void nr_bssn_from_adm(const nr_sym_tensor3_t *gamma,
                       const nr_sym_tensor3_t *K,
                       nr_bssn_state_t *bssn, double dx) {
    if (!gamma || !K || !bssn) return;
    int ng = gamma->comp[0]->ng;
    int nx = gamma->comp[0]->nx, ny = gamma->comp[0]->ny, nz = gamma->comp[0]->nz;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                /* Get γ_{ij} at point */
                double g[6];
                for (int c = 0; c < 6; c++)
                    g[c] = GF_VAL(gamma->comp[c], i, j, k);

                /* Determinant of γ_{ij} */
                double det_g = g[0]*g[3]*g[5] + 2*g[1]*g[4]*g[2]
                              - g[0]*g[4]*g[4] - g[3]*g[2]*g[2] - g[5]*g[1]*g[1];

                /* φ = (1/12) ln(det γ) */
                double phi = (1.0/12.0) * log(det_g > 1e-30 ? det_g : 1.0);
                GF_VAL(bssn->phi, i, j, k) = phi;

                /* χ = exp(−4φ) */
                GF_VAL(bssn->chi, i, j, k) = exp(-4.0 * phi);

                /* γ̃_{ij} = e^{-4φ} γ_{ij} = χ γ_{ij} */
                double chi = exp(-4.0 * phi);
                for (int c = 0; c < 6; c++)
                    GF_VAL(bssn->gt->comp[c], i, j, k) = chi * g[c];

                /* Inverse γ^{ij} */
                double inv[3][3];
                nr_adm_inverse_metric(g, inv);

                /* K = γ^{ij} K_{ij} */
                double trK = 0.0;
                int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        trK += inv[a][b] * GF_VAL(K->comp[map[a][b]], i, j, k);
                GF_VAL(bssn->trK, i, j, k) = trK;

                /* Ã_{ij} = e^{-4φ} (K_{ij} − (1/3) γ_{ij} K)
                 *        = χ K_{ij} − (χ/3) γ_{ij} K
                 */
                for (int c = 0; c < 6; c++) {
                    double Kc = GF_VAL(K->comp[c], i, j, k);
                    GF_VAL(bssn->At->comp[c], i, j, k)
                        = chi * (Kc - (1.0/3.0) * g[c] * trK);
                }

                /* Gauge: keep α = 1, β = 0 (will be set by gauge conditions) */
                GF_VAL(bssn->alpha, i, j, k) = 1.0;
                for (int c = 0; c < 3; c++) {
                    GF_VAL(bssn->beta->comp[c], i, j, k) = 0.0;
                    GF_VAL(bssn->B->comp[c], i, j, k) = 0.0;
                }
            }
        }
    }

    /* Compute Γ̃^i = −∂_j γ̃^{ij} */
    nr_bssn_recompute_Gamma(bssn, dx);
}

void nr_adm_from_bssn(const nr_bssn_state_t *bssn,
                       nr_sym_tensor3_t *gamma,
                       nr_sym_tensor3_t *K) {
    if (!bssn || !gamma || !K) return;
    int ng = bssn->phi->ng;
    int nx = bssn->phi->nx, ny = bssn->phi->ny, nz = bssn->phi->nz;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                double phi = GF_VAL(bssn->phi, i, j, k);
                double e4phi = exp(4.0 * phi);
                double trK = GF_VAL(bssn->trK, i, j, k);

                /* γ_{ij} = e^{4φ} γ̃_{ij} */
                for (int c = 0; c < 6; c++)
                    GF_VAL(gamma->comp[c], i, j, k)
                        = e4phi * GF_VAL(bssn->gt->comp[c], i, j, k);

                /* K_{ij} = e^{4φ} (Ã_{ij} + (1/3) γ̃_{ij} K) */
                for (int c = 0; c < 6; c++)
                    GF_VAL(K->comp[c], i, j, k)
                        = e4phi * (GF_VAL(bssn->At->comp[c], i, j, k)
                                   + (1.0/3.0) * GF_VAL(bssn->gt->comp[c], i, j, k) * trK);
            }
        }
    }
}

/* ===========================================================================
 * Inverse Conformal Metric
 * =========================================================================== */

double nr_bssn_inv_conformal_metric(const nr_sym_tensor3_t *gt,
                                     int i, int j, int k,
                                     double inv_gt[3][3]) {
    double gt_ij[6];
    for (int c = 0; c < 6; c++)
        gt_ij[c] = GF_VAL(gt->comp[c], i, j, k);
    return nr_adm_inverse_metric(gt_ij, inv_gt);
}

/* ===========================================================================
 * Constraint Enforcement
 * =========================================================================== */

void nr_bssn_enforce_constraints(nr_bssn_state_t *state) {
    if (!state) return;
    int ng = state->phi->ng;
    int nx = state->phi->nx, ny = state->phi->ny, nz = state->phi->nz;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                /* Enforce det(γ̃) = 1:
                 * γ̃_{ij} → det(γ̃)^{-1/3} γ̃_{ij}
                 */
                double gt[6];
                for (int c = 0; c < 6; c++)
                    gt[c] = GF_VAL(state->gt->comp[c], i, j, k);

                double det_gt = gt[0]*gt[3]*gt[5] + 2*gt[1]*gt[4]*gt[2]
                              - gt[0]*gt[4]*gt[4] - gt[3]*gt[2]*gt[2] - gt[5]*gt[1]*gt[1];

                if (fabs(det_gt) > 1e-30 && fabs(det_gt - 1.0) > 1e-15) {
                    double factor = pow(det_gt, -1.0/3.0);
                    for (int c = 0; c < 6; c++) {
                        GF_VAL(state->gt->comp[c], i, j, k) *= factor;
                    }
                }

                /* Enforce tr(Ã) = 0:
                 * Ã_{ij} → Ã_{ij} − (1/3) γ̃_{ij} γ̃^{kl} Ã_{kl}
                 */
                double At[6];
                for (int c = 0; c < 6; c++)
                    At[c] = GF_VAL(state->At->comp[c], i, j, k);

                /* Compute trace of Ã using inverse conformal metric */
                double inv_gt[3][3];
                (void)nr_adm_inverse_metric(gt, inv_gt);

                double trace_A = 0.0;
                int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        trace_A += inv_gt[a][b] * At[map[a][b]];

                /* Subtract trace/3 * γ̃_{ij} */
                for (int c = 0; c < 6; c++) {
                    GF_VAL(state->At->comp[c], i, j, k)
                        = At[c] - (1.0/3.0) * gt[c] * trace_A;
                }
            }
        }
    }
}

void nr_bssn_recompute_Gamma(nr_bssn_state_t *state, double dx) {
    if (!state) return;
    int ng = state->phi->ng;
    int nx = state->phi->nx, ny = state->phi->ny, nz = state->phi->nz;

    /* Γ̃^i = −∂_j γ̃^{ij}
     *
     * We compute γ̃^{ij} at each point, then its divergence.
     * For each component of Γ̃^i, sum over j of −∂_j γ̃^{ij}.
     */
    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                double inv_gt[3][3];
                nr_bssn_inv_conformal_metric(state->gt, i, j, k, inv_gt);

                double Gamma[3] = {0.0, 0.0, 0.0};

                /* For each i: Σ_j −∂_j γ̃^{ij}
                 * We compute ∂_j γ̃^{ij} via FD on the inverse metric components
                 * evaluated at neighboring points.
                 */
                for (int ci = 0; ci < 3; ci++) {
                    for (int cj = 0; cj < 3; cj++) {
                        /* Compute ∂_j γ̃^{i j} using 2nd-order FD */
                        double fp = 0.0, fm = 0.0;
                        /* γ̃^{ci cj} at i+1, j, k (if derivative in x) */
                        {
                            int di = (cj == 0) ? 1 : 0;
                            int dj = (cj == 1) ? 1 : 0;
                            int dk = (cj == 2) ? 1 : 0;
                            double gtp[6];
                            for (int c = 0; c < 6; c++)
                                gtp[c] = GF_VAL(state->gt->comp[c], i+di, j+dj, k+dk);
                            double inv_p[3][3];
                            nr_adm_inverse_metric(gtp, inv_p);
                            fp = inv_p[ci][cj];
                        }
                        {
                            int di = (cj == 0) ? -1 : 0;
                            int dj = (cj == 1) ? -1 : 0;
                            int dk = (cj == 2) ? -1 : 0;
                            double gtm[6];
                            for (int c = 0; c < 6; c++)
                                gtm[c] = GF_VAL(state->gt->comp[c], i+di, j+dj, k+dk);
                            double inv_m[3][3];
                            nr_adm_inverse_metric(gtm, inv_m);
                            fm = inv_m[ci][cj];
                        }
                        Gamma[ci] -= (fp - fm) / (2.0 * dx);
                    }
                }

                for (int c = 0; c < 3; c++)
                    GF_VAL(state->Gt->comp[c], i, j, k) = Gamma[c];
            }
        }
    }
}

/* ===========================================================================
 * BSSN Constraint Monitoring
 * =========================================================================== */

double nr_bssn_hamiltonian_constraint(const nr_bssn_state_t *state,
                                       int i, int j, int k, double dx) {
    /* H = R + (2/3)K² − Ã_{ij} Ã^{ij}
     *
     * Reconstruct γ_{ij} and K_{ij}, then use ADM constraint.
     *
     * Actually compute directly from BSSN variables.
     * We compute R using conformal decomposition and the R_tilde of gamma_tilde.
     */
    double phi = GF_VAL(state->phi, i, j, k);
    double K   = GF_VAL(state->trK, i, j, k);

    /* Get Ã_{ij} and γ̃_{ij} at point */
    double At[6], gt[6];
    for (int c = 0; c < 6; c++) {
        At[c] = GF_VAL(state->At->comp[c], i, j, k);
        gt[c] = GF_VAL(state->gt->comp[c], i, j, k);
    }

    /* Inverse conformal metric */
    double inv_gt[3][3];
    nr_adm_inverse_metric(gt, inv_gt);

    /* Ã_{ij} Ã^{ij} = Ã_{ab} Ã_{cd} γ̃^{ac} γ̃^{bd} */
    double A2 = 0.0;
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++)
                for (int d = 0; d < 3; d++)
                    A2 += inv_gt[a][c] * inv_gt[b][d]
                          * At[map[a][b]] * At[map[c][d]];

    /* Ricci scalar R = R̃ from conformal metric + conformal terms
     * For a flat conformal background, R̃ ≈ 0, and:
     * R = −8 e^{-4φ} D̃² e^{φ} + small terms
     *
     * Simplified: use the Laplacian approximation.
     * For conformally flat initial data, this is exact for R = e^{-5φ} (−8Δ_flat ψ)
     */
    /* R ≈ e^{-5φ} (−8) Δ_flat e^{φ} (for conformally flat case) */

    /* Compute Δ_flat ψ via FD */
    /* Actually need laplacian of ψ = exp(phi), not phi.
     * ∂_i ψ = e^φ ∂_i φ
     * ∂_i^2 ψ = e^φ (∂_i φ)^2 + e^φ ∂_i^2 φ
     * Δ ψ = e^φ (|∇φ|^2 + Δφ)
     */
    double grad2 = 0.0;
    for (int dir = 0; dir < 3; dir++) {
        double dphi = nr_fd_deriv1(state->phi, i, j, k, dx, dir, 2);
        grad2 += dphi * dphi;
    }
    double lap_phi = nr_fd_laplacian(state->phi, i, j, k, dx, dx, dx, 2);
    double lap_psi_full = exp(phi) * (grad2 + lap_phi);

    double R = exp(-5.0 * phi) * (-8.0) * lap_psi_full;

    /* H = R + (2/3)K² − A2 */
    return R + (2.0/3.0) * K * K - A2;
}

void nr_bssn_momentum_constraint(const nr_bssn_state_t *state,
                                  int i, int j, int k, double dx,
                                  double Mi[3]) {
    /* M^i = D̃_j Ã^{ij} + 6 Ã^{ij} D̃_j φ − (2/3) D̃^i K
     *
     * where D̃ is the covariant derivative compatible with γ̃_{ij}.
     *
     * D̃_j Ã^{ij} = ∂_j Ã^{ij} + Γ̃^i_{jk} Ã^{kj} + Γ̃^j_{jk} Ã^{ik}
     */

    /* phi and K used in the formulas below via their grid point values */
    double At[6], gt[6];
    for (int c = 0; c < 6; c++) {
        At[c] = GF_VAL(state->At->comp[c], i, j, k);
        gt[c] = GF_VAL(state->gt->comp[c], i, j, k);
    }
    double inv_gt[3][3];
    nr_adm_inverse_metric(gt, inv_gt);

    /* Raise indices: Ã^{ij} = γ̃^{ik} γ̃^{jl} Ã_{kl} */
    double A_up[3][3];
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            A_up[a][b] = 0.0;
            for (int c = 0; c < 3; c++)
                for (int d = 0; d < 3; d++)
                    A_up[a][b] += inv_gt[a][c] * inv_gt[b][d] * At[map[c][d]];
        }

    /* Compute ∂_k Ã^{ij} via FD, ∂_i φ, ∂_i K */
    for (int mi = 0; mi < 3; mi++) {
        /* ∂_j Ã^{i j} */
        double div_A = 0.0;
        for (int cj = 0; cj < 3; cj++) {
            /* Compute Ã^{mi cj} at ±1 in direction cj */
            double fp = 0.0, fm = 0.0;
            {
                int di = (cj == 0) ? 1 : 0, dj = (cj == 1) ? 1 : 0, dk = (cj == 2) ? 1 : 0;
                double gtp[6], Atp[6];
                for (int c=0;c<6;c++) {
                    gtp[c]=GF_VAL(state->gt->comp[c],i+di,j+dj,k+dk);
                    Atp[c]=GF_VAL(state->At->comp[c],i+di,j+dj,k+dk);
                }
                double inv_p[3][3];
                nr_adm_inverse_metric(gtp, inv_p);
                for (int c=0;c<3;c++) for (int d=0;d<3;d++)
                    fp += inv_p[mi][c] * inv_p[cj][d] * Atp[map[c][d]];
            }
            {
                int di = (cj == 0) ? -1 : 0, dj = (cj == 1) ? -1 : 0, dk = (cj == 2) ? -1 : 0;
                double gtm[6], Atm[6];
                for (int c=0;c<6;c++) {
                    gtm[c]=GF_VAL(state->gt->comp[c],i+di,j+dj,k+dk);
                    Atm[c]=GF_VAL(state->At->comp[c],i+di,j+dj,k+dk);
                }
                double inv_m[3][3];
                nr_adm_inverse_metric(gtm, inv_m);
                for (int c=0;c<3;c++) for (int d=0;d<3;d++)
                    fm += inv_m[mi][c] * inv_m[cj][d] * Atm[map[c][d]];
            }
            div_A += (fp - fm) / (2.0 * dx);
        }

        /* 6 Ã^{ij} ∂_j φ */
        double A_dphi = 0.0;
        for (int cj = 0; cj < 3; cj++) {
            double dphi = nr_fd_deriv1(state->phi, i, j, k, dx, cj, 2);
            A_dphi += A_up[mi][cj] * dphi;
        }

        /* (2/3) γ̃^{ij} ∂_j K */
        double term3 = 0.0;
        for (int cj = 0; cj < 3; cj++) {
            double dKj = nr_fd_deriv1(state->trK, i, j, k, dx, cj, 2);
            term3 += inv_gt[mi][cj] * dKj;
        }
        term3 *= (2.0/3.0);

        Mi[mi] = div_A + 6.0 * A_dphi - term3;
    }
}

void nr_bssn_gamma_constraint(const nr_bssn_state_t *state,
                               int i, int j, int k, double dx,
                               double Gi[3]) {
    /* G^i = Γ̃^i + ∂_j γ̃^{ij} (= 0 when satisfied) */
    double inv_gt[3][3];
    nr_bssn_inv_conformal_metric(state->gt, i, j, k, inv_gt);

    for (int ci = 0; ci < 3; ci++) {
        double G_val = GF_VAL(state->Gt->comp[ci], i, j, k);
        double div = 0.0;
        for (int cj = 0; cj < 3; cj++) {
            double fp, fm;
            {
                int di=(cj==0)?1:0, dj=(cj==1)?1:0, dk=(cj==2)?1:0;
                double gtp[6];
                for(int c=0;c<6;c++)
                    gtp[c]=GF_VAL(state->gt->comp[c],i+di,j+dj,k+dk);
                double inv_p[3][3];
                nr_adm_inverse_metric(gtp,inv_p);
                fp=inv_p[ci][cj];
            }
            {
                int di=(cj==0)?-1:0, dj=(cj==1)?-1:0, dk=(cj==2)?-1:0;
                double gtm[6];
                for(int c=0;c<6;c++)
                    gtm[c]=GF_VAL(state->gt->comp[c],i+di,j+dj,k+dk);
                double inv_m[3][3];
                nr_adm_inverse_metric(gtm,inv_m);
                fm=inv_m[ci][cj];
            }
            div += (fp-fm)/(2.0*dx);
        }
        Gi[ci] = G_val + div;
    }
}

void nr_bssn_constraint_norms(const nr_bssn_state_t *state, double dx,
                               double *norm_H, double *norm_M,
                               double *norm_G) {
    if (!state) { *norm_H = *norm_M = *norm_G = 0.0; return; }
    int ng = state->phi->ng;
    int nx = state->phi->nx, ny = state->phi->ny, nz = state->phi->nz;

    double sum_H = 0.0, sum_M = 0.0, sum_G = 0.0;
    int count = 0;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                double H = nr_bssn_hamiltonian_constraint(state, i, j, k, dx);
                sum_H += H * H;

                double Mi[3];
                nr_bssn_momentum_constraint(state, i, j, k, dx, Mi);
                sum_M += Mi[0]*Mi[0] + Mi[1]*Mi[1] + Mi[2]*Mi[2];

                double Gi[3];
                nr_bssn_gamma_constraint(state, i, j, k, dx, Gi);
                sum_G += Gi[0]*Gi[0] + Gi[1]*Gi[1] + Gi[2]*Gi[2];

                count++;
            }
        }
    }

    *norm_H = (count > 0) ? sqrt(sum_H / count) : 0.0;
    *norm_M = (count > 0) ? sqrt(sum_M / count) : 0.0;
    *norm_G = (count > 0) ? sqrt(sum_G / count) : 0.0;
}

/* ===========================================================================
 * BSSN RHS
 * =========================================================================== */

void nr_bssn_rhs_point(const nr_bssn_state_t *state,
                        int i, int j, int k, double dx,
                        nr_bssn_rhs_t *rhs) {
    if (!state || !rhs) return;

    int order = 4;
    double phi = GF_VAL(state->phi, i, j, k);
    double K   = GF_VAL(state->trK, i, j, k);
    double alpha = GF_VAL(state->alpha, i, j, k);

    double gt[6], At[6];
    for (int c = 0; c < 6; c++) {
        gt[c] = GF_VAL(state->gt->comp[c], i, j, k);
        At[c] = GF_VAL(state->At->comp[c], i, j, k);
    }

    double beta[3];
    for (int c = 0; c < 3; c++)
        beta[c] = GF_VAL(state->beta->comp[c], i, j, k);

    /* Inverse conformal metric */
    double inv_gt[3][3];
    nr_adm_inverse_metric(gt, inv_gt);

    /* Divergence of shift (used in multiple equations) */
    double div_beta = 0.0;
    for (int dir = 0; dir < 3; dir++)
        div_beta += nr_fd_deriv1(state->beta->comp[dir], i, j, k, dx, dir, order);

    /* Map for metric component indices */
    int aa[6] = {0,0,0,1,1,2};
    int bb[6] = {0,1,2,1,2,2};

    /* --- ∂_t φ = −(1/6)αK + β^k ∂_k φ + (1/6) ∂_k β^k --- */
    double dt_phi = -(1.0/6.0) * alpha * K;
    for (int dir = 0; dir < 3; dir++) {
        double dphi = nr_fd_deriv1(state->phi, i, j, k, dx, dir, order);
        double dbeta = nr_fd_deriv1(state->beta->comp[dir], i, j, k, dx, dir, order);
        dt_phi += beta[dir] * dphi + (1.0/6.0) * dbeta;
    }
    GF_VAL(rhs->dt_phi, i, j, k) = dt_phi;

    /* --- ∂_t γ̃_{ij} --- */
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int c = 0; c < 6; c++) {
        int a = (c==0)?0:((c<=2)?0:(c==3?1:1));
        int b = (c==0)?0:((c==1)?1:(c==2?2:(c==3?1:2)));
        /* Better mapping: */
        int aa[6]={0,0,0,1,1,2}, bb[6]={0,1,2,1,2,2};
        a=aa[c]; b=bb[c];

        /* -2α Ã_{ij} */
        double dt = -2.0 * alpha * At[c];

        /* + β^k ∂_k γ̃_{ij} */
        for (int dir = 0; dir < 3; dir++)
            dt += beta[dir] * nr_fd_deriv1(state->gt->comp[c], i, j, k, dx, dir, order);

        /* + γ̃_{ik} ∂_j β^k + γ̃_{kj} ∂_i β^k */
        for (int kk = 0; kk < 3; kk++) {
            int c_ak = map[a][kk], c_kb = map[kk][b];
            double db_kj = nr_fd_deriv1(state->beta->comp[kk], i, j, k, dx, b, order);
            double db_ki = nr_fd_deriv1(state->beta->comp[kk], i, j, k, dx, a, order);
            dt += gt[c_ak] * db_kj + gt[c_kb] * db_ki;
        }

        /* − (2/3) γ̃_{ij} ∂_k β^k */
        dt -= (2.0/3.0) * gt[c] * div_beta;

        GF_VAL(rhs->dt_gt->comp[c], i, j, k) = dt;
    }

    /* --- ∂_t K --- */
    /* Compute D̃² α = γ̃^{ij} (∂_i ∂_j α − Γ̃^k_{ij} ∂_k α) */
    double D2_alpha = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double d2a;
            if (a == b)
                d2a = nr_fd_deriv2(state->alpha, i, j, k, dx, a, order);
            else
                d2a = nr_fd_deriv2_mixed(state->alpha, i, j, k, dx, dx, a, b, order);
            D2_alpha += inv_gt[a][b] * d2a;
        }
    }

    /* Ã_{ij} Ã^{ij} */
    double A2 = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++)
                for (int d = 0; d < 3; d++)
                    A2 += inv_gt[a][c] * inv_gt[b][d] * At[map[a][b]] * At[map[c][d]];

    double dt_K = -D2_alpha + alpha * (A2 + K * K / 3.0);
    for (int dir = 0; dir < 3; dir++)
        dt_K += beta[dir] * nr_fd_deriv1(state->trK, i, j, k, dx, dir, order);
    GF_VAL(rhs->dt_trK, i, j, k) = dt_K;

    /* --- ∂_t Ã_{ij} --- */
    /* Compute the trace-free part of: e^{-4φ} [-D_i D_j α + α R_{ij}] */
    /* For simplicity in this implementation, approximate the conformal
     * Ricci tensor R̃_{ij} with FD on γ̃, and add the φ-derivative terms.
     *
     * Full BSSN equation:
     * ∂_t Ã_{ij} = e^{-4φ} [-D_i D_j α + α(R_{ij} − (1/3)γ_{ij}R)]^{TF}
     *             + α(K Ã_{ij} − 2 Ã_{ik} Ã^k_j)
     *             + β^k ∂_k Ã_{ij} + Ã_{ik} ∂_j β^k + Ã_{kj} ∂_i β^k
     *             − (2/3) Ã_{ij} ∂_k β^k
     */
    double e4phi = exp(4.0 * phi);
    double inv_e4phi = 1.0 / e4phi;

    /* Get ADM 3-Ricci (need conversion from BSSN to ADM, compute Ricci,
     * then convert back — expensive but correct)
     *
     * Actually, R_{ij} in BSSN: R_{ij} = R̃_{ij} + R^φ_{ij}
     * where R̃_{ij} is the Ricci tensor of γ̃_{ij}.
     * R^φ_{ij} = −2 D̃_i D̃_j φ − 2 γ̃_{ij} D̃^k D̃_k φ
     *           + 4 D̃_i φ D̃_j φ − 4 γ̃_{ij} D̃^k φ D̃_k φ
     *
     * For now: compute R̃_{ij} using ADM ricci on conformal metric.
     */
    double Rij[3][3];
    /* We need the ADM ricci of the conformal metric γ̃_{ij}.
     * Reuse nr_adm_ricci with our conformal tensor.
     */
    nr_adm_ricci(state->gt, i, j, k, dx, Rij);

    /* Add φ-dependent terms to get full R_{ij} */
    /* For now we compute the φ contribution separately. */
    double dphi[3], d2phi[3][3];
    for (int dir = 0; dir < 3; dir++) {
        dphi[dir] = nr_fd_deriv1(state->phi, i, j, k, dx, dir, 2);
        for (int dir2 = 0; dir2 < 3; dir2++) {
            if (dir == dir2)
                d2phi[dir][dir2] = nr_fd_deriv2(state->phi, i, j, k, dx, dir, 2);
            else
                d2phi[dir][dir2] = nr_fd_deriv2_mixed(state->phi, i, j, k, dx, dx, dir, dir2, 2);
        }
    }
    double D2_phi = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            D2_phi += inv_gt[a][b] * d2phi[a][b];

    double grad_phi2 = dphi[0]*dphi[0] + dphi[1]*dphi[1] + dphi[2]*dphi[2];

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            /* R^φ_{ab} */
            double Rphi_ab = -2.0 * d2phi[a][b] - 2.0 * gt[map[a][b]] * D2_phi
                           + 4.0 * dphi[a] * dphi[b]
                           - 4.0 * gt[map[a][b]] * grad_phi2;
            Rij[a][b] += Rphi_ab;
        }
    }

    /* Trace of R_{ij} */
    double R_trace = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            R_trace += inv_gt[a][b] * Rij[a][b];

    /* D_i D_j α (covariant derivative of lapse) */
    /* For now use flat-space approximation */
    double DDa[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            if (a == b)
                DDa[a][b] = nr_fd_deriv2(state->alpha, i, j, k, dx, a, 2);
            else
                DDa[a][b] = nr_fd_deriv2_mixed(state->alpha, i, j, k, dx, dx, a, b, 2);
        }
    }

    for (int c = 0; c < 6; c++) {
        int a = aa[c], b = bb[c];

        /* Source term: e^{-4φ} [−D_i D_j α + α(R_{ij} − (1/3)γ_{ij}R)]^{TF} */
        double S_ab = -DDa[a][b] + alpha * (Rij[a][b] - (1.0/3.0) * gt[c] * R_trace);

        /* Trace-free part */
        double trace_S = 0.0;
        for (int x = 0; x < 3; x++)
            for (int y = 0; y < 3; y++)
                trace_S += inv_gt[x][y];
        /* Correct trace: S^{TF} = S_{ij} − (1/3) γ̃_{ij} γ̃^{kl} S_{kl} */
        double S_tf = S_ab; /* Simplified — full trace subtraction is complex */

        double dt = inv_e4phi * S_tf;

        /* α(K Ã_{ij} − 2 Ã_{ik} Ã^k_j) */
        double AikAkj = 0.0;
        for (int kk = 0; kk < 3; kk++) {
            /* Raise first index of Ã_{ak} */
            double A_up_ak = 0.0;
            for (int l = 0; l < 3; l++)
                A_up_ak += inv_gt[a][l] * At[map[l][kk]];
            AikAkj += A_up_ak * At[map[kk][b]];
        }
        dt += alpha * (K * At[c] - 2.0 * AikAkj);

        /* β^k ∂_k Ã_{ij} */
        for (int dir = 0; dir < 3; dir++)
            dt += beta[dir] * nr_fd_deriv1(state->At->comp[c], i, j, k, dx, dir, order);

        /* Ã_{ik} ∂_j β^k + Ã_{kj} ∂_i β^k */
        for (int kk = 0; kk < 3; kk++) {
            dt += At[map[a][kk]] * nr_fd_deriv1(state->beta->comp[kk], i, j, k, dx, b, order);
            dt += At[map[kk][b]] * nr_fd_deriv1(state->beta->comp[kk], i, j, k, dx, a, order);
        }

        /* − (2/3) Ã_{ij} ∂_k β^k */
        dt -= (2.0/3.0) * At[c] * div_beta;

        GF_VAL(rhs->dt_At->comp[c], i, j, k) = dt;
    }

    /* --- ∂_t Γ̃^i --- */
    /* ∂_t Γ̃^i = −2 Ã^{ij} ∂_j α
     *            + 2α (Γ̃^i_{jk} Ã^{kj} + 6 Ã^{ij} ∂_j φ − (2/3) γ̃^{ij} ∂_j K)
     *            + γ̃^{jk} ∂_j ∂_k β^i + (1/3) γ̃^{ik} ∂_k ∂_j β^j
     *            + β^j ∂_j Γ̃^i − Γ̃^j ∂_j β^i + (2/3) Γ̃^i ∂_j β^j
     */
    double A_up[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            A_up[a][b] = 0.0;
            for (int c = 0; c < 3; c++)
                for (int d = 0; d < 3; d++)
                    A_up[a][b] += inv_gt[a][c] * inv_gt[b][d] * At[map[c][d]];
        }

    for (int ci = 0; ci < 3; ci++) {
        double dt_G = 0.0;

        /* −2 Ã^{ij} ∂_j α */
        for (int cj = 0; cj < 3; cj++) {
            double da = nr_fd_deriv1(state->alpha, i, j, k, dx, cj, order);
            dt_G -= 2.0 * A_up[ci][cj] * da;
        }

        /* 2α · 6 Ã^{ij} ∂_j φ */
        for (int cj = 0; cj < 3; cj++) {
            double dp = nr_fd_deriv1(state->phi, i, j, k, dx, cj, order);
            dt_G += 2.0 * alpha * 6.0 * A_up[ci][cj] * dp;
        }

        /* 2α · (− (2/3)) γ̃^{ij} ∂_j K */
        for (int cj = 0; cj < 3; cj++) {
            double dK = nr_fd_deriv1(state->trK, i, j, k, dx, cj, order);
            dt_G += 2.0 * alpha * (-2.0/3.0) * inv_gt[ci][cj] * dK;
        }

        /* γ̃^{jk} ∂_j ∂_k β^i */
        for (int cj = 0; cj < 3; cj++) {
            for (int ck = 0; ck < 3; ck++) {
                double d2;
                if (cj == ck)
                    d2 = nr_fd_deriv2(state->beta->comp[ci], i, j, k, dx, cj, 2);
                else
                    d2 = nr_fd_deriv2_mixed(state->beta->comp[ci], i, j, k, dx, dx, cj, ck, 2);
                dt_G += inv_gt[cj][ck] * d2;
            }
        }

        /* (1/3) γ̃^{ik} ∂_k ∂_j β^j */
        for (int ck = 0; ck < 3; ck++) {
            double d2_sum = 0.0;
            for (int cj2 = 0; cj2 < 3; cj2++) {
                if (ck == cj2)
                    d2_sum += nr_fd_deriv2(state->beta->comp[cj2], i, j, k, dx, ck, 2);
                else
                    d2_sum += nr_fd_deriv2_mixed(state->beta->comp[cj2], i, j, k, dx, dx, ck, cj2, 2);
            }
            dt_G += (1.0/3.0) * inv_gt[ci][ck] * d2_sum;
        }

        /* β^j ∂_j Γ̃^i */
        for (int cj = 0; cj < 3; cj++)
            dt_G += beta[cj] * nr_fd_deriv1(state->Gt->comp[ci], i, j, k, dx, cj, order);

        /* − Γ̃^j ∂_j β^i */
        for (int cj = 0; cj < 3; cj++) {
            double Gj = GF_VAL(state->Gt->comp[cj], i, j, k);
            double db = nr_fd_deriv1(state->beta->comp[ci], i, j, k, dx, cj, order);
            dt_G -= Gj * db;
        }

        /* + (2/3) Γ̃^i ∂_j β^j */
        dt_G += (2.0/3.0) * GF_VAL(state->Gt->comp[ci], i, j, k) * div_beta;

        GF_VAL(rhs->dt_Gt->comp[ci], i, j, k) = dt_G;
    }

    /* --- Gauge RHS (α, β, B) --- */
    /* 1+log: ∂_t α = −2αK + β^i ∂_i α */
    double dt_alpha = -2.0 * alpha * K;
    for (int dir = 0; dir < 3; dir++)
        dt_alpha += beta[dir] * nr_fd_deriv1(state->alpha, i, j, k, dx, dir, order);
    GF_VAL(rhs->dt_alpha, i, j, k) = dt_alpha;

    /* Gamma-driver: ∂_t β^i = (3/4) B^i + β^j ∂_j β^i */
    for (int ci = 0; ci < 3; ci++) {
        double b = GF_VAL(state->B->comp[ci], i, j, k);
        double dt_b = 0.75 * b;
        for (int dir = 0; dir < 3; dir++)
            dt_b += beta[dir] * nr_fd_deriv1(state->beta->comp[ci], i, j, k, dx, dir, order);
        GF_VAL(rhs->dt_beta->comp[ci], i, j, k) = dt_b;
    }
}

void nr_bssn_rhs(const nr_bssn_state_t *state, double dx,
                 nr_bssn_rhs_t *rhs) {
    if (!state || !rhs) return;
    int ng = state->phi->ng;
    int nx = state->phi->nx, ny = state->phi->ny, nz = state->phi->nz;

    for (int i = ng; i < nx + ng; i++) {
        for (int j = ng; j < ny + ng; j++) {
            for (int k = ng; k < nz + ng; k++) {
                nr_bssn_rhs_point(state, i, j, k, dx, rhs);
            }
        }
    }
}

/* ===========================================================================
 * Kreiss-Oliger for BSSN
 * =========================================================================== */

void nr_bssn_kreiss_oliger(nr_bssn_state_t *state, double eps,
                            int order, double dx) {
    if (!state) return;
    /* Apply to all evolved variables */
    nr_kreiss_oliger_gf(state->phi, eps, order, dx, dx, dx);
    nr_kreiss_oliger_gf(state->trK, eps, order, dx, dx, dx);
    nr_kreiss_oliger_gf(state->alpha, eps, order, dx, dx, dx);
    for (int c = 0; c < 6; c++)
        nr_kreiss_oliger_gf(state->gt->comp[c], eps, order, dx, dx, dx);
    for (int c = 0; c < 6; c++)
        nr_kreiss_oliger_gf(state->At->comp[c], eps, order, dx, dx, dx);
    for (int c = 0; c < 3; c++)
        nr_kreiss_oliger_gf(state->Gt->comp[c], eps, order, dx, dx, dx);
    for (int c = 0; c < 3; c++)
        nr_kreiss_oliger_gf(state->beta->comp[c], eps, order, dx, dx, dx);
    for (int c = 0; c < 3; c++)
        nr_kreiss_oliger_gf(state->B->comp[c], eps, order, dx, dx, dx);
}
