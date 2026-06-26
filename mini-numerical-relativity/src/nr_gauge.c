/**
 * @file nr_gauge.c
 * @brief Gauge conditions for numerical relativity: slicing and shift.
 *
 * Gauge conditions determine how the spatial hypersurfaces foliate
 * spacetime (slicing) and how spatial coordinates propagate (shift).
 * These are freely specifiable in GR — the physics is independent
 * of the gauge choice.
 *
 * Standard gauge conditions in NR:
 *   1+log slicing:    ∂_t α = −2αK + β^i ∂_i α
 *   Gamma-driver:     ∂_t β^i = (3/4) B^i + β^j ∂_j β^i
 *                     ∂_t B^i = ∂_t Γ̃^i − η B^i + β^j ∂_j B^i − β^j ∂_j Γ̃^i
 *
 * References:
 *   - Bona, Masso, Seidel, Stela, PRD 56, 3405 (1997) — 1+log family
 *   - Alcubierre et al., PRD 67, 084023 (2003) — Gamma-driver
 *   - Baumgarte & Shapiro (2010), Ch. 4 and Sec. 11.5
 *   - Campanelli et al., PRL 96, 111101 (2006)
 *
 * Knowledge Coverage:
 *   L1: Lapse α, shift β^i, 1+log slicing, Gamma-driver
 *   L2: Gauge freedom, slicing conditions, spatial coordinate conditions
 *   L4: Gauge evolution equations
 *   L5: Implementation of gauge RHS
 *   L6: Trumpet slice (stationary 1+log slicing of Schwarzschild)
 */

#include "nr_bssn.h"
#include "nr_utils.h"
#include <math.h>

/* ===========================================================================
 * L4: 1+Log Slicing
 *
 * The 1+log family: ∂_t α = −α² f(α) K + β^i ∂_i α
 *
 * Standard choice f(α) = 2/α gives:
 *   ∂_t α = −2α K + β^i ∂_i α
 *
 * This is a "singularity-avoiding" slicing: as the lapse collapses near
 * the singularity, proper time evolution slows, effectively freezing the
 * interior evolution while the exterior proceeds normally.
 *
 * Reference: Bona et al., PRD 56, 3405 (1997), Eq. (18)
 * =========================================================================== */

/**
 * Compute ∂_t α for the standard 1+log slicing condition.
 *
 * @param alpha     Lapse function on grid
 * @param trK       Trace of extrinsic curvature
 * @param beta[3]   Shift vector
 * @param i,j,k     Grid indices
 * @param dx        Grid spacing
 * @return          ∂_t α at the grid point
 *
 * Complexity: O(1).
 */
double nr_gauge_1pluslog_point(const nr_gf_t *alpha, const nr_gf_t *trK,
                                const nr_gf_t *beta[3],
                                int i, int j, int k, double dx) {
    double a = alpha->data[i * alpha->stride_x + j * alpha->stride_y + k];
    double K = trK->data[i * trK->stride_x + j * trK->stride_y + k];

    /* β^i ∂_i α */
    double beta_d_alpha = 0.0;
    for (int dir = 0; dir < 3; dir++) {
        double b = beta[dir]->data[i * beta[dir]->stride_x
                                  + j * beta[dir]->stride_y + k];
        double da = nr_fd_deriv1(alpha, i, j, k, dx, dir, 4);
        beta_d_alpha += b * da;
    }

    return -2.0 * a * K + beta_d_alpha;
}

/**
 * Compute ∂_t α for the 1+log slicing over the entire interior grid.
 *
 * Complexity: O(N^3).
 */
void nr_gauge_1pluslog_rhs(const nr_gf_t *alpha, const nr_gf_t *trK,
                            const nr_gf_t *beta[3],
                            nr_gf_t *dt_alpha, double dx) {
    if (!alpha || !trK || !dt_alpha) return;
    int ng = alpha->ng;
    for (int i = ng; i < alpha->nx + ng; i++) {
        for (int j = ng; j < alpha->ny + ng; j++) {
            for (int k = ng; k < alpha->nz + ng; k++) {
                dt_alpha->data[i * dt_alpha->stride_x + j * dt_alpha->stride_y + k]
                    = nr_gauge_1pluslog_point(alpha, trK, beta, i, j, k, dx);
            }
        }
    }
}

/* ===========================================================================
 * L4: Gamma-Driver Shift Condition
 *
 * The Gamma-driver condition is a hyperbolic shift condition that
 * mimics the harmonic shift in the BSSN system:
 *
 *   ∂_t β^i  = (3/4) B^i + β^j ∂_j β^i
 *   ∂_t B^i  = ∂_t Γ̃^i − η B^i + β^j ∂_j B^i − β^j ∂_j Γ̃^i
 *
 * The parameter η controls the damping timescale of the shift.
 * Typical values: η = 1/(2M) for single BH, η = 2/M for binary.
 *
 * Reference: Alcubierre et al., PRD 67, 084023 (2003), Eq. (36)-(37)
 * =========================================================================== */

/**
 * Compute ∂_t β^i and ∂_t B^i for the Gamma-driver at a grid point.
 *
 * @param beta[3]   Shift vector
 * @param B[3]      Auxiliary B^i field
 * @param dt_Gt[3]  ∂_t Γ̃^i (from BSSN evolution)
 * @param Gt[3]     Conformal connection Γ̃^i
 * @param eta       Damping parameter
 * @param i,j,k     Grid indices
 * @param dx        Grid spacing
 * @param dt_beta[3]  Output: ∂_t β^i
 * @param dt_B[3]     Output: ∂_t B^i
 *
 * Complexity: O(1).
 */
void nr_gauge_gamma_driver_point(const nr_gf_t *beta[3],
                                  const nr_gf_t *B[3],
                                  const double dt_Gt[3],
                                  const nr_gf_t *Gt[3],
                                  double eta,
                                  int i, int j, int k, double dx,
                                  double dt_beta[3], double dt_B[3]) {
    /* Get β^i and B^i at this point */
    double b[3], Bval[3];
    for (int c = 0; c < 3; c++) {
        b[c]    = beta[c]->data[i * beta[c]->stride_x + j * beta[c]->stride_y + k];
        Bval[c] = B[c]->data[i * B[c]->stride_x + j * B[c]->stride_y + k];
    }

    /* β^j ∂_j β^i */
    double b_db[3] = {0.0, 0.0, 0.0};
    /* β^j ∂_j B^i */
    double b_dB[3] = {0.0, 0.0, 0.0};
    /* β^j ∂_j Γ̃^i */
    double b_dG[3] = {0.0, 0.0, 0.0};

    for (int ci = 0; ci < 3; ci++) {
        for (int dir = 0; dir < 3; dir++) {
            double bj = b[dir];
            b_db[ci] += bj * nr_fd_deriv1(beta[ci], i, j, k, dx, dir, 4);
            b_dB[ci] += bj * nr_fd_deriv1(B[ci], i, j, k, dx, dir, 4);
            b_dG[ci] += bj * nr_fd_deriv1(Gt[ci], i, j, k, dx, dir, 4);
        }
    }

    /* ∂_t β^i = (3/4) B^i + β^j ∂_j β^i */
    for (int ci = 0; ci < 3; ci++) {
        dt_beta[ci] = 0.75 * Bval[ci] + b_db[ci];
    }

    /* ∂_t B^i = ∂_t Γ̃^i − η B^i + β^j ∂_j B^i − β^j ∂_j Γ̃^i */
    for (int ci = 0; ci < 3; ci++) {
        dt_B[ci] = dt_Gt[ci] - eta * Bval[ci] + b_dB[ci] - b_dG[ci];
    }
}

/**
 * Compute Gamma-driver RHS over the entire interior grid.
 *
 * Complexity: O(N^3).
 */
void nr_gauge_gamma_driver_rhs(const nr_bssn_state_t *state,
                                const nr_bssn_rhs_t *rhs,
                                double eta, double dx) {
    if (!state || !rhs) return;
    int ng = state->alpha->ng;

    /* Wrap beta and B as gf* arrays for the point function */
    const nr_gf_t *beta_arr[3] = {
        state->beta->comp[0], state->beta->comp[1], state->beta->comp[2]
    };
    const nr_gf_t *B_arr[3] = {
        state->B->comp[0], state->B->comp[1], state->B->comp[2]
    };
    const nr_gf_t *Gt_arr[3] = {
        state->Gt->comp[0], state->Gt->comp[1], state->Gt->comp[2]
    };

    nr_gf_t *dt_beta_arr[3] = {
        rhs->dt_beta->comp[0], rhs->dt_beta->comp[1], rhs->dt_beta->comp[2]
    };
    nr_gf_t *dt_B_arr[3] = {
        rhs->dt_B->comp[0], rhs->dt_B->comp[1], rhs->dt_B->comp[2]
    };

    for (int i = ng; i < state->alpha->nx + ng; i++) {
        for (int j = ng; j < state->alpha->ny + ng; j++) {
            for (int k = ng; k < state->alpha->nz + ng; k++) {
                /* Get dt_Gt from rhs */
                double dt_Gt_local[3];
                for (int c = 0; c < 3; c++) {
                    dt_Gt_local[c] = rhs->dt_Gt->comp[c]->data[
                        i * rhs->dt_Gt->comp[c]->stride_x
                        + j * rhs->dt_Gt->comp[c]->stride_y + k];
                }

                double dt_b[3], dt_B[3];
                nr_gauge_gamma_driver_point(
                    beta_arr, B_arr, dt_Gt_local, Gt_arr, eta,
                    i, j, k, dx, dt_b, dt_B);

                for (int c = 0; c < 3; c++) {
                    dt_beta_arr[c]->data[i * dt_beta_arr[c]->stride_x
                                         + j * dt_beta_arr[c]->stride_y + k]
                        = dt_b[c];
                    dt_B_arr[c]->data[i * dt_B_arr[c]->stride_x
                                     + j * dt_B_arr[c]->stride_y + k]
                        = dt_B[c];
                }
            }
        }
    }
}

/* ===========================================================================
 * L6: Trumpet Slice — Stationary 1+Log of Schwarzschild
 *
 * The stationary end-state of 1+log slicing of Schwarzschild is the
 * "trumpet" slice. Key features:
 *   - The slice reaches a finite areal radius R_0 as r → 0
 *   - R_0 ≈ 1.312 M for the standard 1+log condition
 *   - The lapse α → 0 at r → 0 (puncture field freezing)
 *
 * Radial profile in isotropic radius r:
 *   α(r) = sqrt(1 − 2M/R(r))
 *   where R(r) = r(1 + M/(2r))^2
 *
 * The conformal factor:
 *   ψ(r) = 1 + M/(2r) + C/r  (approximately, with C ≈ 0)
 *
 * Reference: Hannam et al., PRL 99, 241102 (2007)
 * =========================================================================== */

/**
 * Compute the trumpet lapse α(r) for Schwarzschild.
 *
 * @param r     Isotropic coordinate radius (must be > 0)
 * @param M     ADM mass
 * @return      Lapse α(r)
 *
 * Complexity: O(1).
 */
double nr_gauge_trumpet_alpha(double r, double M) {
    if (r < 1e-15) return 0.0;
    double R = r * (1.0 + M / (2.0 * r)) * (1.0 + M / (2.0 * r));
    double interior = 1.0 - 2.0 * M / R;
    if (interior < 0.0) return 0.0;
    return sqrt(interior);
}

/**
 * Compute the trumpet conformal factor ψ(r).
 *
 * In the exact trumpet: ψ(r) = 1 + M/(2r).
 * This is the same as isotropic Schwarzschild.
 *
 * @param r     Isotropic radius
 * @param M     ADM mass
 * @return      Conformal factor ψ(r)
 *
 * Complexity: O(1).
 */
double nr_gauge_trumpet_psi(double r, double M) {
    if (r < 1e-15) return 1.0 + M / (2.0 * 1e-15);  /* large but finite */
    return 1.0 + M / (2.0 * r);
}

/**
 * Compute the trumpet areal radius R(r) from isotropic radius.
 *
 * R(r) = r ψ(r)^2 = r (1 + M/(2r))^2
 *
 * Complexity: O(1).
 */
double nr_gauge_trumpet_areal_radius(double r, double M) {
    double psi = nr_gauge_trumpet_psi(r, M);
    return r * psi * psi;
}

/**
 * Set up the trumpet initial data on a BSSN grid.
 *
 * Fills the BSSN state with Schwarzschild in trumpet coordinates.
 *
 * @param state    BSSN state to fill
 * @param grid     Grid definition
 * @param M        ADM mass
 *
 * Complexity: O(N^3).
 */
void nr_gauge_set_trumpet(nr_bssn_state_t *state,
                           const nr_grid_t *grid, double M) {
    if (!state || !grid) return;

    int ng = grid->ng;
    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];
                double r = sqrt(x*x + y*y + z*z);

                /* Avoid division by zero */
                if (r < 1e-15) {
                    /* At the puncture: lapse → 0, psi large */
                    int idx_a = i*state->alpha->stride_x + j*state->alpha->stride_y + k;
                    state->alpha->data[idx_a] = 0.0;
                    state->phi->data[idx_a] = 4.0;  /* arbitrary large value */
                    state->trK->data[idx_a] = 0.0;
                    for (int c = 0; c < 6; c++) {
                        int cidx = i*state->gt->comp[c]->stride_x
                                  + j*state->gt->comp[c]->stride_y + k;
                        state->gt->comp[c]->data[cidx] = (c==0||c==3||c==5) ? 1.0 : 0.0;
                        state->At->comp[c]->data[cidx] = 0.0;
                    }
                    for (int c = 0; c < 3; c++) {
                        int cidx = i*state->beta->comp[c]->stride_x
                                  + j*state->beta->comp[c]->stride_y + k;
                        state->beta->comp[c]->data[cidx] = 0.0;
                        state->B->comp[c]->data[cidx] = 0.0;
                        state->Gt->comp[c]->data[cidx] = 0.0;
                    }
                    continue;
                }

                double psi = nr_gauge_trumpet_psi(r, M);
                double alpha = nr_gauge_trumpet_alpha(r, M);

                /* BSSN conformal factor: φ = ln(ψ) */
                state->phi->data[i*state->phi->stride_x
                                 + j*state->phi->stride_y + k] = log(psi);

                /* Lapse */
                state->alpha->data[i*state->alpha->stride_x
                                   + j*state->alpha->stride_y + k] = alpha;

                /* Conformal metric: γ̃_{ij} = ψ^{-4} γ_{ij} = ψ^{-4} · ψ^4 δ_{ij} = δ_{ij} */
                double *gt[6];
                for (int c = 0; c < 6; c++) gt[c] = &state->gt->comp[c]->data[
                    i*state->gt->comp[c]->stride_x + j*state->gt->comp[c]->stride_y + k];
                *gt[0] = 1.0; *gt[1] = 0.0; *gt[2] = 0.0;
                *gt[3] = 1.0; *gt[4] = 0.0; *gt[5] = 1.0;

                /* Extrinsic curvature: K_{ij} = 0 for time-symmetric data */
                for (int c = 0; c < 6; c++)
                    state->At->comp[c]->data[i*state->At->comp[c]->stride_x
                                             + j*state->At->comp[c]->stride_y + k] = 0.0;
                state->trK->data[i*state->trK->stride_x
                                 + j*state->trK->stride_y + k] = 0.0;

                /* Shift and auxiliary fields */
                for (int c = 0; c < 3; c++) {
                    state->beta->comp[c]->data[i*state->beta->comp[c]->stride_x
                                              + j*state->beta->comp[c]->stride_y + k] = 0.0;
                    state->B->comp[c]->data[i*state->B->comp[c]->stride_x
                                          + j*state->B->comp[c]->stride_y + k] = 0.0;
                    state->Gt->comp[c]->data[i*state->Gt->comp[c]->stride_x
                                           + j*state->Gt->comp[c]->stride_y + k] = 0.0;
                }
            }
        }
    }
}

/* ===========================================================================
 * L2: Maximal Slicing
 *
 * Maximal slicing sets K = 0 and ∂_t K = 0, leading to an elliptic
 * equation for the lapse:
 *
 *   D_i D^i α = α K_{ij} K^{ij}
 *
 * This is computationally expensive (3D elliptic solve at each time step)
 * but provides ideal singularity avoidance. 1+log slicing approximates
 * maximal slicing with a local (parabolic) evolution equation.
 *
 * Reference: Baumgarte & Shapiro (2010), Sec. 4.1; Estabrook et al., PRD 7, 2814 (1973)
 * =========================================================================== */

/**
 * Check if the current slice is approximately maximal (K ≈ 0).
 *
 * @param trK    Trace of extrinsic curvature on the grid
 * @param tol    Tolerance for |K|_∞
 * @return       true if max|K| < tol
 *
 * Complexity: O(N^3).
 */
bool nr_gauge_is_maximal(const nr_gf_t *trK, double tol) {
    return nr_gf_maxabs(trK) < tol;
}

/* ===========================================================================
 * L6: Harmonic Slicing
 *
 * Harmonic slicing: □ t = 0, which gives:
 *   ∂_t α = −α^2 K + β^i ∂_i α
 *
 * This is the f(α)=1/α case of the Bona-Masso family.
 * Unlike 1+log, harmonic slicing does not avoid singularities
 * (it can reach them in finite proper time). Used primarily for
 * analytic studies and as a simple test case.
 *
 * Reference: Bona et al., PRD 56, 3405 (1997), Eq. (16)
 * =========================================================================== */

/**
 * Compute ∂_t α for harmonic slicing.
 *
 * @param alpha     Lapse grid function
 * @param trK       Trace of K
 * @param beta[3]   Shift vector
 * @param i,j,k     Grid indices
 * @param dx        Grid spacing
 * @return          ∂_t α
 *
 * Complexity: O(1).
 */
double nr_gauge_harmonic_point(const nr_gf_t *alpha, const nr_gf_t *trK,
                                const nr_gf_t *beta[3],
                                int i, int j, int k, double dx) {
    double a = alpha->data[i * alpha->stride_x + j * alpha->stride_y + k];
    double K = trK->data[i * trK->stride_x + j * trK->stride_y + k];

    double beta_d_alpha = 0.0;
    for (int dir = 0; dir < 3; dir++) {
        double b = beta[dir]->data[i * beta[dir]->stride_x
                                  + j * beta[dir]->stride_y + k];
        double da = nr_fd_deriv1(alpha, i, j, k, dx, dir, 4);
        beta_d_alpha += b * da;
    }

    return -a * a * K + beta_d_alpha;
}

/**
 * Compute the Bona-Masso slicing parameter f(α).
 *
 * The general Bona-Masso slicing condition is:
 *   ∂_t α = −α² f(α) K + β^i ∂_i α
 *
 * Standard choices:
 *   f(α) = 1       → harmonic slicing
 *   f(α) = 2/α     → 1+log slicing
 *   f(α) = 1 + c/α² → shock-avoiding slicing
 *
 * @param alpha     Lapse
 * @param slicing_type  0=harmonic, 1=1+log, 2=shock-avoiding
 * @param c         Parameter for shock-avoiding
 * @return          f(α)
 *
 * Complexity: O(1).
 * Reference: Bona et al., PRD 56, 3405 (1997)
 */
double nr_gauge_bona_masso_f(double alpha, int slicing_type, double c) {
    if (slicing_type == 0) {
        /* Harmonic: f(α) = 1 */
        return 1.0;
    } else if (slicing_type == 1) {
        /* 1+log: f(α) = 2/α */
        if (alpha < 1e-15) return 2.0 / 1e-15;
        return 2.0 / alpha;
    } else {
        /* Shock-avoiding: f(α) = 1 + c/α² */
        if (alpha < 1e-15) return 1.0 + c / (1e-30);
        return 1.0 + c / (alpha * alpha);
    }
}
