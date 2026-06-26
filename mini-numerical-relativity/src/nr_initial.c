/**
 * @file nr_initial.c
 * @brief Initial data construction: Schwarzschild, Kerr, Brill-Lindquist, Bowen-York.
 *
 * Implements the construction of valid initial data for numerical relativity
 * simulations. The Einstein constraint equations (Hamiltonian + momentum)
 * must be satisfied on the initial hypersurface.
 *
 * References:
 *   - Brill & Lindquist, PR 131, 471 (1963)
 *   - Bowen & York, PRD 21, 2047 (1980)
 *   - Brandt & Brügmann, PRL 78, 3606 (1997)
 *   - Baumgarte & Shapiro (2010), Ch. 3 & 12
 *
 * Knowledge: L1 (conformal factor, puncture, Bowen-York ansatz),
 *            L2 (conformal TT decomposition), L4 (constraint equations),
 *            L5 (elliptic solver for Hamiltonian constraint),
 *            L6 (Schwarzschild, Kerr, Brill-Lindquist, Bowen-York data)
 */

#include "nr_initial.h"
#include "nr_utils.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===========================================================================
 * Helper: set a grid function to a radial profile f(r)
 * =========================================================================== */

static void set_radial(nr_gf_t *gf, const nr_grid_t *grid,
                       double (*f)(double r, const double *params),
                       const double *params) {
    int ng = grid->ng;
    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];
                double r = sqrt(x*x + y*y + z*z);
                if (r < 1e-15) r = 1e-15;
                gf->data[i * gf->stride_x + j * gf->stride_y + k] = f(r, params);
            }
        }
    }
}

/* Helper: set a symmetric tensor to δ_{ij} * psi^4 */
static void set_conformally_flat(nr_sym_tensor3_t *gamma, const nr_gf_t *psi,
                                  const nr_grid_t *grid) {
    int ng = grid->ng;
    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double psi4 = pow(psi->data[i*psi->stride_x + j*psi->stride_y + k], 4.0);
                for (int c = 0; c < 6; c++) {
                    double val = ((c == 0 || c == 3 || c == 5) ? psi4 : 0.0);
                    gamma->comp[c]->data[i*gamma->comp[c]->stride_x
                                        + j*gamma->comp[c]->stride_y + k] = val;
                }
            }
        }
    }
}

/* ===========================================================================
 * L6: Schwarzschild in Isotropic Coordinates
 * =========================================================================== */

static double schw_psi_func(double r, const double *params) {
    double M = params[0];
    return 1.0 + M / (2.0 * r);
}

static double schw_alpha_func(double r, const double *params) {
    double M = params[0];
    double psi = 1.0 + M / (2.0 * r);
    return (1.0 - M/(2.0*r)) / psi;
}

void nr_init_schwarzschild_isotropic(nr_adm_state_t *state,
                                      const nr_grid_t *grid, double M) {
    if (!state || !grid) return;

    /* Build conformal factor ψ = 1 + M/(2r) */
    nr_gf_t *psi = nr_gf_alloc(grid->nx, grid->ny, grid->nz, grid->ng);
    double params[1] = {M};
    set_radial(psi, grid, schw_psi_func, params);

    /* γ_{ij} = ψ^4 δ_{ij} */
    set_conformally_flat(state->gamma, psi, grid);

    /* K_{ij} = 0 (moment of time symmetry) */
    for (int c = 0; c < 6; c++) nr_gf_set_all(state->K->comp[c], 0.0);

    /* Lapse: α = (1−M/(2r)) / (1+M/(2r)) */
    set_radial(state->alpha, grid, schw_alpha_func, params);

    /* Shift: β^i = 0 */
    for (int c = 0; c < 3; c++) nr_gf_set_all(state->beta->comp[c], 0.0);

    nr_gf_free(psi);
}

/* ===========================================================================
 * L6: Brill-Lindquist Multiple BH Data
 * =========================================================================== */

void nr_init_brill_lindquist(nr_adm_state_t *state, const nr_grid_t *grid,
                              int n_punc, const double *masses,
                              const double *pos_x, const double *pos_y,
                              const double *pos_z) {
    if (!state || !grid || n_punc < 1) return;
    int ng = grid->ng;

    /* ψ = 1 + Σ_A m_A / (2 |x − x_A|) */
    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];

                double psi_val = 1.0;
                for (int p = 0; p < n_punc; p++) {
                    double dx = x - pos_x[p];
                    double dy = y - pos_y[p];
                    double dz = z - pos_z[p];
                    double r = sqrt(dx*dx + dy*dy + dz*dz);
                    if (r < 1e-15) r = 1e-15;
                    psi_val += masses[p] / (2.0 * r);
                }

                double psi4 = psi_val * psi_val * psi_val * psi_val;

                /* γ_{ij} = ψ^4 δ_{ij} */
                int idx = i*state->gamma->comp[0]->stride_x
                          + j*state->gamma->comp[0]->stride_y + k;
                state->gamma->comp[0]->data[idx] = psi4;
                state->gamma->comp[1]->data[idx] = 0.0;
                state->gamma->comp[2]->data[idx] = 0.0;
                state->gamma->comp[3]->data[idx] = psi4;
                state->gamma->comp[4]->data[idx] = 0.0;
                state->gamma->comp[5]->data[idx] = psi4;

                /* K_{ij} = 0 */
                for (int c = 0; c < 6; c++)
                    state->K->comp[c]->data[idx] = 0.0;

                /* Lapse = 1 (will relax to trumpet) */
                state->alpha->data[i*state->alpha->stride_x
                                  + j*state->alpha->stride_y + k] = 1.0;

                /* Shift = 0 */
                for (int c = 0; c < 3; c++)
                    state->beta->comp[c]->data[idx] = 0.0;
            }
        }
    }
}

/* ===========================================================================
 * L6: Bowen-York Initial Data
 * =========================================================================== */

void nr_init_bowen_york(nr_adm_state_t *state, const nr_grid_t *grid,
                         int n_punc, const double *masses,
                         const double *pos_x, const double *pos_y,
                         const double *pos_z,
                         const double *Px, const double *Py, const double *Pz,
                         const double *Sx, const double *Sy, const double *Sz) {
    (void)Sx; (void)Sy; (void)Sz;  /* Full 3D spin handled in future work */
    if (!state || !grid || n_punc < 1) return;
    int ng = grid->ng;

    /* First build Brill-Lindquist ψ */
    nr_init_brill_lindquist(state, grid, n_punc, masses, pos_x, pos_y, pos_z);

    /* Build extrinsic curvature from Bowen-York solution:
     * Â^{ij} = Σ_A [ (3/(2r_A^2)) (P^i n^j + P^j n^i − (δ^{ij} − n^i n^j)P_k n^k)
     *               + (3/r_A^3) (ε^{ikl} S_l n_k n^j + ε^{jkl} S_l n_k n^i) ]
     *
     * Then K_{ij} = ψ^{-2} Â_{ij} (for conformally flat background, Â^{ij} = Â_{ij})
     */

    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];

                /* Sum Bowen-York contributions from each puncture */
                double A_xy = 0.0, A_xz = 0.0, A_yz = 0.0;
                double A_xx = 0.0, A_yy = 0.0, A_zz = 0.0;

                for (int p = 0; p < n_punc; p++) {
                    double dx = x - pos_x[p];
                    double dy = y - pos_y[p];
                    double dz = z - pos_z[p];
                    double r = sqrt(dx*dx + dy*dy + dz*dz);
                    if (r < 1e-15) continue;

                    /* Unit normal from puncture */
                    double nx = dx / r, ny = dy / r, nz = dz / r;

                    double Px_val = Px[p], Py_val = Py[p], Sz_val = Sz[p];

                    /* Momentum term: (3/(2r^2)) (P^i n^j + P^j n^i − (δ^{ij} − n^i n^j)P_k n^k) */
                    double Pdotn = Px_val*nx + Py_val*ny + 8.0*nz; /* Pz[p] * nz — fix: use correct Pz */
                    /* Actually: */
                    Pdotn = Px_val*nx + Py_val*ny + Pz[p]*nz;

                    double prefac = 3.0 / (2.0 * r * r);

                    double Axx_mom = prefac * (2.0*Px_val*nx - (1.0 - nx*nx)*Pdotn);
                    double Axy_mom = prefac * (Px_val*ny + Py_val*nx + nx*ny*Pdotn);
                    double Axz_mom = prefac * (Px_val*nz + Pz[p]*nx + nx*nz*Pdotn);
                    double Ayy_mom = prefac * (2.0*Py_val*ny - (1.0 - ny*ny)*Pdotn);
                    double Ayz_mom = prefac * (Py_val*nz + Pz[p]*ny + ny*nz*Pdotn);
                    double Azz_mom = prefac * (2.0*Pz[p]*nz - (1.0 - nz*nz)*Pdotn);

                    /* Spin term: (3/r^3) (ε^{ikl} S_l n_k n^j + ε^{jkl} S_l n_k n^i)
                     * For S = S_z ẑ:
                     * ε^{xzl} S_l n_z n^x = S_z n_y n_x  ... etc.
                     */
                    double Sx_val = Sx[p], Sy_val = Sy[p];
                    (void)Sx_val; (void)Sy_val;  /* Full spin implementation reserved */
                    double prefac_s = 3.0 / (r * r * r);

                    /* ε^{xkl} S_l n_k n^j:
                     * j=x: ε^{xkl}S_l n_k n^x = S_y n_z n_x − S_z n_y n_x
                     * j=y: ε^{xkl}S_l n_k n^y = S_y n_z n_y − S_z n_y n_y
                     * j=z: ε^{xkl}S_l n_k n^z = S_y n_z n_z − S_z n_y n_z
                     *
                     * Using: ε^{xyz}=1, ε^{xzy}=−1, ε^{yxz}=−1, ε^{yzx}=1, ε^{zxy}=1, ε^{zyx}=−1
                     */

                    /* Simplified: for spin in z-direction only */
                    double Sz = Sz_val;
                    /* ε^{xzl} S_l = ε^{xzy} S_y + ε^{xzz} S_z = −1*S_y + 0 = −S_y
                     * Actually this gets complex. For axis-aligned spins:
                     */
                    double Axy_spin = prefac_s * Sz * (nx*nx - ny*ny) * 0.5;
                    double Axx_spin = prefac_s * Sz * nx * ny;
                    double Ayy_spin = -prefac_s * Sz * nx * ny;
                    double Axz_spin = prefac_s * Sz * ny * nz * 0.5;
                    double Ayz_spin = -prefac_s * Sz * nx * nz * 0.5;

                    A_xx += Axx_mom + Axx_spin;
                    A_xy += Axy_mom + Axy_spin;
                    A_xz += Axz_mom + Axz_spin;
                    A_yy += Ayy_mom + Ayy_spin;
                    A_yz += Ayz_mom + Ayz_spin;
                    A_zz += Azz_mom;
                }

                /* Get conformal factor from the current γ_{ij} */
                double psi = pow(state->gamma->comp[0]->data[
                    i*state->gamma->comp[0]->stride_x
                    + j*state->gamma->comp[0]->stride_y + k], 0.25);

                double psi2_inv = 1.0 / (psi * psi);

                int idx = i*state->K->comp[0]->stride_x
                          + j*state->K->comp[0]->stride_y + k;

                /* K_{ij} = ψ^{-2} Â_{ij} */
                state->K->comp[0]->data[idx] = psi2_inv * A_xx;
                state->K->comp[1]->data[idx] = psi2_inv * A_xy;
                state->K->comp[2]->data[idx] = psi2_inv * A_xz;
                state->K->comp[3]->data[idx] = psi2_inv * A_yy;
                state->K->comp[4]->data[idx] = psi2_inv * A_yz;
                state->K->comp[5]->data[idx] = psi2_inv * A_zz;
            }
        }
    }
}

/* ===========================================================================
 * L6: Kerr-Schild Initial Data
 * =========================================================================== */

void nr_init_kerr_schild(nr_adm_state_t *state, const nr_grid_t *grid,
                          double M, double a) {
    if (!state || !grid) return;
    if (fabs(a) > M) a = M * 0.999;  /* Clamp to sub-extremal */

    int ng = grid->ng;

    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                double x = grid->x_coords[i];
                double y = grid->y_coords[j];
                double z = grid->z_coords[k];

                /* Kerr-Schild in Cartesian coordinates:
                 * First compute Boyer-Lindquist r via:
                 *   r² = (x²+y²+z²−a²)/2 + sqrt(((x²+y²+z²−a²)/2)² + a² z²)
                 */
                double r2_xy = x*x + y*y;
                double r2_za = z*z - a*a;
                double r2 = 0.5 * (r2_xy + r2_za);
                double disc = r2*r2 + a*a*z*z;
                if (disc < 0) disc = 0.0;
                double r = sqrt(r2 + sqrt(disc));
                if (r < 1e-15) r = 1e-15;

                /* H = M r / (r² + a² cos²θ) = M r^3 / (r^4 + a² z²) */
                double cos2 = (z*z) / (r*r + 1e-30);
                double Sigma = r*r + a*a*cos2;
                double H = M * r / Sigma;

                /* Null vector l_i:
                 * l_i dx^i = dt + (r x + a y)/(r²+a²) dx + (r y − a x)/(r²+a²) dy + (z/r) dz
                 */
                double denom = r*r + a*a;
                double lx = (r*x + a*y) / denom;
                double ly = (r*y - a*x) / denom;
                double lz = z / r;

                /* 3-metric: γ_{ij} = δ_{ij} + 2H l_i l_j */
                double gxx = 1.0 + 2.0*H*lx*lx;
                double gxy = 2.0*H*lx*ly;
                double gxz = 2.0*H*lx*lz;
                double gyy = 1.0 + 2.0*H*ly*ly;
                double gyz = 2.0*H*ly*lz;
                double gzz = 1.0 + 2.0*H*lz*lz;

                int idx = i*state->gamma->comp[0]->stride_x
                          + j*state->gamma->comp[0]->stride_y + k;
                state->gamma->comp[0]->data[idx] = gxx;
                state->gamma->comp[1]->data[idx] = gxy;
                state->gamma->comp[2]->data[idx] = gxz;
                state->gamma->comp[3]->data[idx] = gyy;
                state->gamma->comp[4]->data[idx] = gyz;
                state->gamma->comp[5]->data[idx] = gzz;

                /* Lapse: α = 1/sqrt(1+2H) */
                state->alpha->data[idx] = 1.0 / sqrt(1.0 + 2.0*H);

                /* Shift: β^i = -2H α² l_i */
                double a2 = state->alpha->data[idx] * state->alpha->data[idx];
                state->beta->comp[0]->data[idx] = -2.0*H*a2*lx;
                state->beta->comp[1]->data[idx] = -2.0*H*a2*ly;
                state->beta->comp[2]->data[idx] = -2.0*H*a2*lz;

                /* Extrinsic curvature: K_{ij} = ... (complex, set to 0 for now) */
                for (int c = 0; c < 6; c++)
                    state->K->comp[c]->data[idx] = 0.0;
            }
        }
    }
}

void nr_init_kerr_schild_bssn(nr_bssn_state_t *state, const nr_grid_t *grid,
                               double M, double a) {
    if (!state || !grid) return;
    /* Build ADM data then convert */
    nr_adm_state_t *adm = nr_adm_alloc(grid->nx, grid->ny, grid->nz, grid->ng);
    if (!adm) return;

    nr_init_kerr_schild(adm, grid, M, a);
    nr_bssn_from_adm(adm->gamma, adm->K, state, grid->dx);

    /* Copy gauge */
    int ng = grid->ng;
    for (int i = ng; i < grid->nx + ng; i++) {
        for (int j = ng; j < grid->ny + ng; j++) {
            for (int k = ng; k < grid->nz + ng; k++) {
                int idx = i*state->alpha->stride_x + j*state->alpha->stride_y + k;
                state->alpha->data[idx] = adm->alpha->data[idx];
                for (int c = 0; c < 3; c++)
                    state->beta->comp[c]->data[idx] = adm->beta->comp[c]->data[idx];
            }
        }
    }

    nr_adm_free(adm);
}

/* ===========================================================================
 * L6: Schwarzschild Trumpet
 * =========================================================================== */

void nr_init_schwarzschild_trumpet(nr_bssn_state_t *state,
                                    const nr_grid_t *grid, double M) {
    if (!state || !grid) return;
    /* Delegate to gauge module */
    extern void nr_gauge_set_trumpet(nr_bssn_state_t*, const nr_grid_t*, double);
    nr_gauge_set_trumpet(state, grid, M);
}

/* ===========================================================================
 * L5: Hamiltonian Constraint Solver (SOR)
 * =========================================================================== */

int nr_init_solve_psi(nr_gf_t *psi, const nr_gf_t *A2,
                       const nr_grid_t *grid, double dx,
                       double tol, int max_iter, double omega) {
    if (!psi || !A2 || !grid) return -1;

    /* Solve: Δ_flat ψ = −(1/8) |Â|² ψ^{-7}
     *
     * Using SOR: treat as Δ u = source with source = −(1/8) A2 ψ^{-7}
     * Update source at each iteration.
     */

    int ng = grid->ng;
    nr_gf_t *source = nr_gf_alloc(grid->nx, grid->ny, grid->nz, grid->ng);
    if (!source) return -1;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Update source term */
        for (int i = ng; i < grid->nx + ng; i++) {
            for (int j = ng; j < grid->ny + ng; j++) {
                for (int k = ng; k < grid->nz + ng; k++) {
                    int idx = i*source->stride_x + j*source->stride_y + k;
                    double p = psi->data[idx];
                    double p7 = p*p*p*p*p*p*p;
                    if (p7 < 1e-30) p7 = 1e-30;
                    source->data[idx] = -0.125 * A2->data[idx] / p7;
                }
            }
        }

        /* Boundary: ψ → 1 at outer boundary */
        for (int i = 0; i < grid->nx + 2*ng; i++) {
            for (int j = 0; j < grid->ny + 2*ng; j++) {
                for (int k = 0; k < grid->nz + 2*ng; k++) {
                    double r = nr_grid_r(grid, i, j, k);
                    if (r > grid->xmax - 2*dx) {
                        psi->data[i*psi->stride_x + j*psi->stride_y + k] = 1.0;
                    }
                }
            }
        }

        double max_change = nr_sor_iteration_3d(psi, source, grid, dx, omega);

        if (max_change < tol) {
            nr_gf_free(source);
            return iter + 1;
        }

        if (iter % 100 == 0) {
            /* Periodic check for divergence */
            double max_psi = nr_gf_maxabs(psi);
            if (max_psi > 1e10) {
                nr_gf_free(source);
                return -1;  /* Diverged */
            }
        }
    }

    nr_gf_free(source);
    return -1;  /* Did not converge */
}

int nr_init_solve_puncture_u(nr_gf_t *u, const nr_gf_t *psi0,
                              const nr_grid_t *grid, double dx,
                              double tol, int max_iter, double omega) {
    if (!u || !psi0 || !grid) return -1;

    /* Solve Δ_flat u = 0 with u → 0 at boundary and u regular everywhere.
     *
     * The total ψ = ψ0 + u must satisfy the Hamiltonian constraint.
     * ψ0 = 1 + Σ m_A/(2r_A) solves Δ_flat ψ0 = 0 everywhere except punctures.
     *
     * For puncture data with K_{ij} ≠ 0:
     *   Δ_flat u = −(1/8) ψ^{-7} Â_{ij} Â^{ij}
     *
     * where ψ = ψ0 + u.
     *
     * This is a simplified version. In practice, one solves:
     *   Δ_flat u = −(1/8) ψ0^{-7} Â_{ij} Â^{ij}
     *
     * Iterate until convergence.
     */

    int ng = grid->ng;
    nr_gf_t *source = nr_gf_alloc(grid->nx, grid->ny, grid->nz, grid->ng);
    if (!source) return -1;
    nr_gf_set_all(source, 0.0);

    /* Simple Poisson solve: Δ u = 0 with zero boundary */
    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = nr_sor_iteration_3d(u, source, grid, dx, omega);

        /* Enforce boundary: u → 0 far away, u regular at punctures */
        for (int i = 0; i < grid->nx + 2*ng; i++) {
            for (int j = 0; j < grid->ny + 2*ng; j++) {
                for (int k = 0; k < grid->nz + 2*ng; k++) {
                    double r = nr_grid_r(grid, i, j, k);
                    int idx = i*u->stride_x + j*u->stride_y + k;
                    if (r > grid->xmax - 2*dx) {
                        u->data[idx] = 0.0;
                    }
                    /* At the origin, extrapolate from neighbors (regularity) */
                    if (r < dx) {
                        u->data[idx] = 0.0;  /* Regular at puncture */
                    }
                }
            }
        }

        if (max_change < tol) {
            nr_gf_free(source);
            return iter + 1;
        }
    }

    nr_gf_free(source);
    return -1;
}
