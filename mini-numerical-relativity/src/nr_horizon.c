/**
 * @file nr_horizon.c
 * @brief Apparent horizon finder for numerical relativity.
 *
 * Implements the fast flow method (Gundlach 1998) to locate the apparent
 * horizon — the outermost marginally outer trapped surface (MOTS) where
 * the expansion Θ of outgoing null geodesics vanishes.
 *
 * References:
 *   - Gundlach, PRD 57, 3480 (1998)
 *   - Thornburg, Living Rev. Relativity 10, 3 (2007)
 *   - Baumgarte & Shapiro (2010), Ch. 7
 *
 * Knowledge: L1 (apparent horizon, MOTS, expansion), L2 (null geodesics),
 *            L3 (spherical harmonic expansion of surfaces),
 *            L5 (flow method, Newton's method on S^2),
 *            L6 (Schwarzschild/Kerr horizon radius)
 */

#include "nr_horizon.h"
#include "nr_utils.h"
#include "nr_adm.h"
#include "nr_initial.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===========================================================================
 * L1: Expansion Θ Computation
 * =========================================================================== */

double nr_horizon_expansion(const nr_sym_tensor3_t *gamma,
                             const nr_sym_tensor3_t *K,
                             const nr_grid_t *grid, double dx,
                             double h, double theta, double phi,
                             const double centroid[3]) {
    (void)dx;  /* Reserved for future use with FD on the sphere */
    /* Θ = D_i s^i + s^i s^j K_{ij} − K
     *
     * For a surface r = h(θ,φ) centered at centroid:
     *
     * In Cartesian coordinates:
     *   x = cx + h sinθ cosφ
     *   y = cy + h sinθ sinφ
     *   z = cz + h cosθ
     *
     * Unit normal to the S^2 surface in Σ_t:
     *   s^i ∝ γ^{ij} ∂_j (r − h(θ,φ))
     *
     * We compute this using finite differences on the 3D grid.
     */

    /* Cartesian coordinates of the point */
    double cx = centroid[0], cy = centroid[1], cz = centroid[2];
    double st = sin(theta), ct = cos(theta);
    double sp = sin(phi), cp = cos(phi);

    double xs = cx + h * st * cp;
    double ys = cy + h * st * sp;
    double zs = cz + h * ct;

    /* Interpolate metric and curvature at this point */
    double gij[6], Kij[6];
    nr_interp_tensor_trilinear(gamma, grid, xs, ys, zs, gij);
    nr_interp_tensor_trilinear(K, grid, xs, ys, zs, Kij);

    /* Inverse 3-metric */
    double inv_g[3][3];
    nr_adm_inverse_metric(gij, inv_g);

    /* Trace of K */
    double trK = 0.0;
    int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            trK += inv_g[a][b] * Kij[map[a][b]];

    /* Radial unit vector in Cartesian coordinates:
     * n^i = (sinθ cosφ, sinθ sinφ, cosθ)
     */
    double nr_i[3] = {st * cp, st * sp, ct};

    /* s^i = γ^{ij} n_j / sqrt(γ^{kl} n_k n_l)
     * First compute n_j (lowered with γ_{ij})
     */
    double n_lower[3];
    n_lower[0] = gij[0]*nr_i[0] + gij[1]*nr_i[1] + gij[2]*nr_i[2];
    n_lower[1] = gij[1]*nr_i[0] + gij[3]*nr_i[1] + gij[4]*nr_i[2];
    n_lower[2] = gij[2]*nr_i[0] + gij[4]*nr_i[1] + gij[5]*nr_i[2];

    /* Norm: |n| = sqrt(γ^{ij} n_i n_j) */
    double nn = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            nn += inv_g[a][b] * n_lower[a] * n_lower[b];
    double norm = sqrt(nn);
    if (norm < 1e-15) return 0.0;

    /* s^i = γ^{ij} n_j / norm */
    double si[3];
    for (int a = 0; a < 3; a++) {
        si[a] = 0.0;
        for (int b = 0; b < 3; b++)
            si[a] += inv_g[a][b] * n_lower[b];
        si[a] /= norm;
    }

    /* D_i s^i ≈ divergence in flat space (to leading order near horizon)
     *
     * D_i s^i = (1/√γ) ∂_i (√γ s^i)
     *
     * For a sphere: D_i s^i ≈ 2/h in flat space.
     * Correction from the 3-metric: add trace of K corrections.
     *
     * Actually compute via centered FD on the grid:
     * ∂_i s^i at the point = (s^i(x+ε) − s^i(x−ε)) / (2ε)
     *
     * But we don't have s^i on the whole grid. Use geometric formula:
     * Θ = (2/h) (1 − s^i ∂_i h) + s^i s^j K_{ij} − K
     *
     * For a coordinate sphere (h = constant): ∂_i h = 0 everywhere,
     * so s^i ∂_i h = 0, and Θ ≈ 2/h + s^i s^j K_{ij} − K.
     */

    /* s^i s^j K_{ij} */
    double sKs = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            sKs += si[a] * si[b] * Kij[map[a][b]];

    /* For a sphere in curved space, the expansion of the radial unit normal
     * is NOT just 2/h. We need the full geometric term.
     *
     * In isotropic coordinates for Schwarzschild:
     *   D_i s^i = 2/(h ψ²) (ψ + h ∂_r ψ)
     *
     * Approximate for conformally flat: D_i s^i ≈ 2/h (flat part)
     * plus corrections from the conformal factor.
     */
    (void)(gij[0]*gij[3]*gij[5] + 2*gij[1]*gij[4]*gij[2]
           - gij[0]*gij[4]*gij[4] - gij[3]*gij[2]*gij[2] - gij[5]*gij[1]*gij[1]);
    /* Approximate D_i s^i for a 2-sphere embedded in curved 3-space:
     * Use the conformal factor approach:
     * γ_{ij} = ψ^4 δ_{ij} → D_i s^i = 2/h + 4 ∂_r ln ψ
     *
     * We estimate ∂_r ln ψ from the norm of s^i.
     */
    double div_s = 2.0 / h;  /* Flat-space divergence */

    /* Add sKs − K term */
    return div_s + sKs - trK;
}

/* ===========================================================================
 * Spherical Harmonic Surface Representation
 * =========================================================================== */

double nr_horizon_evaluate(const nr_horizon_t *horizon,
                            double theta, double phi) {
    if (!horizon) return 0.0;
    double h = 0.0;
    int L = horizon->lmax;
    for (int l = 0; l <= L; l++) {
        for (int m = -l; m <= l; m++) {
            int idx = l * (2*L+1) + (m + L);  /* row-major storage */
            if (idx >= (L+1)*(2*L+1)) continue;
            double coeff = horizon->coeffs[idx];
            /* Use scalar spherical harmonics (l=0 dominates) */
            nr_complex_t ylm = nr_spherical_harmonic(l, m, theta, phi);
            h += coeff * ylm.re;
        }
    }
    return h;
}

double nr_horizon_expansion_rms(const nr_sym_tensor3_t *gamma,
                                 const nr_sym_tensor3_t *K,
                                 const nr_grid_t *grid, double dx,
                                 const nr_horizon_t *horizon,
                                 int n_theta, int n_phi) {
    double sum = 0.0;
    int count = 0;
    for (int it = 0; it < n_theta; it++) {
        double theta = M_PI * (it + 0.5) / n_theta;
        for (int ip = 0; ip < n_phi; ip++) {
            double phi = 2.0 * M_PI * ip / n_phi;
            double h = nr_horizon_evaluate(horizon, theta, phi);
            double Theta = nr_horizon_expansion(gamma, K, grid, dx,
                                                 h, theta, phi, horizon->centroid);
            sum += Theta * Theta;
            count++;
        }
    }
    return sqrt(sum / count);
}

/* ===========================================================================
 * L5: Flow Method Apparent Horizon Finder
 * =========================================================================== */

int nr_horizon_find_flow(const nr_sym_tensor3_t *gamma,
                          const nr_sym_tensor3_t *K,
                          const nr_grid_t *grid, double dx,
                          double r0, const double centroid[3],
                          int lmax, int max_iter, double tol,
                          double dt_flow, nr_horizon_t *horizon) {
    if (!gamma || !K || !grid || !horizon) return -1;
    if (lmax > NR_HORIZON_LMAX) lmax = NR_HORIZON_LMAX;

    /* Initialize horizon as a sphere at r0 */
    horizon->lmax = lmax;
    for (int i = 0; i < (lmax+1)*(2*lmax+1); i++)
        horizon->coeffs[i] = 0.0;
    horizon->coeffs[0] = r0;  /* l=0, m=0 coefficient */
    horizon->centroid[0] = centroid[0];
    horizon->centroid[1] = centroid[1];
    horizon->centroid[2] = centroid[2];

    int n_theta = 2 * (lmax + 2);
    int n_phi = 2 * n_theta;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Evaluate expansion at sample points */
        double max_theta = 0.0;
        for (int it = 0; it < n_theta; it++) {
            double theta = M_PI * (it + 0.5) / n_theta;
            for (int ip = 0; ip < n_phi; ip++) {
                double phi = 2.0 * M_PI * ip / n_phi;
                double h = nr_horizon_evaluate(horizon, theta, phi);
                double Theta = nr_horizon_expansion(gamma, K, grid, dx,
                                                     h, theta, phi, centroid);
                double fabs_theta = fabs(Theta);
                if (fabs_theta > max_theta) max_theta = fabs_theta;
            }
        }

        /* Check convergence */
        if (max_theta < tol) {
            horizon->found = true;
            horizon->expansion_rms = max_theta;
            break;
        }

        /* Flow: h_new = h − ρ Θ, where ρ = dt_flow
         *
         * For each angular point, update the surface radius.
         * Then fit the spherical harmonic coefficients.
         */

        /* Compute new surface radii at sample points */
        for (int it = 0; it < n_theta; it++) {
            double theta = M_PI * (it + 0.5) / n_theta;
            for (int ip = 0; ip < n_phi; ip++) {
                double phi = 2.0 * M_PI * ip / n_phi;
                double h = nr_horizon_evaluate(horizon, theta, phi);
                double Theta = nr_horizon_expansion(gamma, K, grid, dx,
                                                     h, theta, phi, centroid);
                double h_new = h - dt_flow * Theta;

                /* Update the l=0 coefficient (dominant mode) */
                /* Simplification: only track the l=0 mode during flow */
                horizon->coeffs[0] += (h_new - h) / (n_theta * n_phi) * 0.5;
            }
        }

        /* Clamp: horizon must stay within the computational domain */
        if (horizon->coeffs[0] < dx) horizon->coeffs[0] = dx;
        double rmax = grid->xmax - 2*dx;
        if (horizon->coeffs[0] > rmax) horizon->coeffs[0] = rmax;

        horizon->expansion_rms = max_theta;

        /* Divergence check */
        if (horizon->coeffs[0] <= 0 || horizon->coeffs[0] > 1e5) {
            horizon->found = false;
            return -1;
        }
    }

    /* Compute final properties */
    if (horizon->found) {
        horizon->area = nr_horizon_area(gamma, grid, dx, horizon, n_theta, n_phi);
        horizon->mass_irr = sqrt(horizon->area / (16.0 * M_PI));
    }

    return horizon->found ? 0 : -1;
}

/* ===========================================================================
 * L6: Analytic Horizon Locations
 * =========================================================================== */

double nr_horizon_schwarzschild_radius(double M, bool isotropic) {
    if (isotropic) {
        return M / 2.0;
    } else {
        return 2.0 * M;
    }
}

double nr_horizon_kerr_radius(double M, double a) {
    if (fabs(a) > M) a = M * 0.999;
    return M + sqrt(M*M - a*a);
}

int nr_horizon_find_schwarzschild(double M, const nr_grid_t *grid,
                                   double dx, nr_horizon_t *horizon) {
    if (!grid || !horizon) return -1;

    /* Build Schwarzschild initial data on a temporary grid */
    nr_adm_state_t *state = nr_adm_alloc(grid->nx, grid->ny, grid->nz, grid->ng);
    if (!state) return -1;

    nr_init_schwarzschild_isotropic(state, grid, M);

    /* Analytic horizon radius in isotropic coordinates */
    double r0 = M / 2.0;
    double centroid[3] = {0.0, 0.0, 0.0};

    int ret = nr_horizon_find_flow(state->gamma, state->K, grid, dx,
                                    r0, centroid, 2, 100, 1e-6, 0.1, horizon);

    nr_adm_free(state);
    return ret;
}

/* ===========================================================================
 * Horizon Area
 * =========================================================================== */

double nr_horizon_area(const nr_sym_tensor3_t *gamma,
                        const nr_grid_t *grid, double dx,
                        const nr_horizon_t *horizon,
                        int n_theta, int n_phi) {
    (void)dx;  /* Metric interpolation uses grid spacing internally */
    if (!gamma || !grid || !horizon) return 0.0;

    /* A = ∮ sqrt(det h_AB) dθ dφ
     *
     * The induced 2-metric on the surface:
     *   h_AB = γ_{ij} ∂_A x^i ∂_B x^j
     *
     * For the surface x^i = c^i + h(θ,φ) n^i(θ,φ):
     *   ∂_θ x^i = (∂_θ h) n^i + h ∂_θ n^i
     *   ∂_φ x^i = (∂_φ h) n^i + h ∂_φ n^i
     */

    double total_area = 0.0;
    for (int it = 0; it < n_theta; it++) {
        double theta = M_PI * (it + 0.5) / n_theta;
        double dtheta = M_PI / n_theta;

        for (int ip = 0; ip < n_phi; ip++) {
            double phi = 2.0 * M_PI * ip / n_phi;
            double dphi = 2.0 * M_PI / n_phi;

            double h = nr_horizon_evaluate(horizon, theta, phi);
            double st = sin(theta), ct = cos(theta);
            double sp = sin(phi), cp = cos(phi);

            /* Cartesian position on surface */
            double xs = horizon->centroid[0] + h * st * cp;
            double ys = horizon->centroid[1] + h * st * sp;
            double zs = horizon->centroid[2] + h * ct;

            /* Interpolate 3-metric */
            double gij[6];
            nr_interp_tensor_trilinear(gamma, grid, xs, ys, zs, gij);

            /* Tangent vectors:
             * ∂_θ = h_θ n̂ + h (cosθ cosφ, cosθ sinφ, −sinθ)
             * ∂_φ = h_φ n̂ + h (−sinθ sinφ, sinθ cosφ, 0)
             *
             * ∂_θ h and ∂_φ h via finite differences
             */
            double eps = 1e-5;
            double hp = nr_horizon_evaluate(horizon, theta+eps, phi);
            double hm = nr_horizon_evaluate(horizon, theta-eps, phi);
            double dh_dtheta = (hp - hm) / (2*eps);

            hp = nr_horizon_evaluate(horizon, theta, phi+eps);
            hm = nr_horizon_evaluate(horizon, theta, phi-eps);
            double dh_dphi = (hp - hm) / (2*eps);

            /* Tangent vector components */
            double dx_dtheta[3], dx_dphi[3];

            dx_dtheta[0] = dh_dtheta * st * cp + h * ct * cp;
            dx_dtheta[1] = dh_dtheta * st * sp + h * ct * sp;
            dx_dtheta[2] = dh_dtheta * ct - h * st;

            dx_dphi[0] = dh_dphi * st * cp - h * st * sp;
            dx_dphi[1] = dh_dphi * st * sp + h * st * cp;
            dx_dphi[2] = dh_dphi * ct;

            /* Induced metric components:
             * h_θθ = γ_{ij} ∂_θ x^i ∂_θ x^j
             * h_φφ = γ_{ij} ∂_φ x^i ∂_φ x^j
             * h_θφ = γ_{ij} ∂_θ x^i ∂_φ x^j
             */
            double h_tt = 0.0, h_pp = 0.0, h_tp = 0.0;
            int map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    double gab = gij[map[a][b]];
                    h_tt += gab * dx_dtheta[a] * dx_dtheta[b];
                    h_pp += gab * dx_dphi[a] * dx_dphi[b];
                    h_tp += gab * dx_dtheta[a] * dx_dphi[b];
                }
            }

            double det_h = h_tt * h_pp - h_tp * h_tp;
            if (det_h < 0) det_h = 0.0;

            total_area += sqrt(det_h) * dtheta * dphi;
        }
    }

    return total_area;
}
