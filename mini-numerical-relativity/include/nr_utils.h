/**
 * @file nr_utils.h
 * @brief Numerical utilities for finite differences, integration, and interpolation
 *
 * Core numerical methods used throughout numerical relativity:
 *   - Finite difference stencils (2nd, 4th, 6th, 8th order)
 *   - Runge-Kutta time integration (RK2, RK4)
 *   - Trilinear interpolation on 3D grids
 *   - Legendre polynomials and spherical harmonics
 *   - SOR elliptic solver utilities
 *
 * References:
 *   - Fornberg, "A Practical Guide to Pseudospectral Methods" (1998)
 *   - LeVeque, "Finite Difference Methods for Ordinary and Partial Differential Equations" (2007)
 *   - Baumgarte & Shapiro (2010), Appendix B
 *   - Alcubierre (2008), Ch. 4
 *
 * Knowledge Coverage:
 *   L5: Finite difference stencils, RK4, interpolation, elliptic solvers
 *   L3: Discrete derivative operators on 3D meshes
 */

#ifndef NR_UTILS_H
#define NR_UTILS_H

#include "nr_grid.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
 * L5: Finite Difference Stencils
 *
 * Centered finite difference approximations for first and second derivatives.
 *
 * For a uniform grid with spacing h:
 *
 * 2nd-order:
 *   ∂_x f_i = (f_{i+1} − f_{i-1}) / (2h) + O(h^2)
 *   ∂_xx f_i = (f_{i+1} − 2f_i + f_{i-1}) / h^2 + O(h^2)
 *
 * 4th-order:
 *   ∂_x f_i = (−f_{i+2} + 8f_{i+1} − 8f_{i-1} + f_{i-2}) / (12h) + O(h^4)
 *   ∂_xx f_i = (−f_{i+2} + 16f_{i+1} − 30f_i + 16f_{i-1} − f_{i-2}) / (12h^2) + O(h^4)
 *
 * 6th-order:
 *   ∂_x f_i = (f_{i+3} − 9f_{i+2} + 45f_{i+1} − 45f_{i-1} + 9f_{i-2} − f_{i-3}) / (60h) + O(h^6)
 *   ∂_xx f_i = (2f_{i+3} − 27f_{i+2} + 270f_{i+1} − 490f_i + 270f_{i-1} − 27f_{i-2} + 2f_{i-3}) / (180h^2)
 *
 * 8th-order:
 *   ∂_x f_i = (3f_{i+4} − 32f_{i+3} + 168f_{i+2} − 672f_{i+1}
 *              + 672f_{i-1} − 168f_{i-2} + 32f_{i-3} − 3f_{i-4}) / (840h) + O(h^8)
 *
 * =========================================================================== */

/**
 * L5: Compute first derivative ∂_x f at a grid point using centered FD.
 *
 * @param gf     Grid function
 * @param i,j,k  Grid indices
 * @param dx     Grid spacing in x
 * @param dir    0=x derivative, 1=y derivative, 2=z derivative
 * @param order  FD order (2, 4, 6, or 8)
 * @return       ∂_i f at the grid point
 *
 * Uses one-sided stencils near boundaries when ghost zones are insufficient.
 * Complexity: O(1).
 */
double nr_fd_deriv1(const nr_gf_t *gf, int i, int j, int k,
                     double dx, int dir, int order);

/**
 * L5: Compute second derivative ∂²_xx f at a grid point using centered FD.
 *
 * @param gf     Grid function
 * @param i,j,k  Grid indices
 * @param dx     Grid spacing in x
 * @param dir    0=∂²/∂x², 1=∂²/∂y², 2=∂²/∂z²
 * @param order  FD order (2, 4, 6, or 8)
 * @return       ∂²_ii f at the grid point
 *
 * Complexity: O(1).
 */
double nr_fd_deriv2(const nr_gf_t *gf, int i, int j, int k,
                     double dx, int dir, int order);

/**
 * L5: Compute mixed second derivative ∂²_{xy} f at a grid point.
 *
 * Uses successive application of first-derivative stencils:
 *   ∂²_{xy} f = ∂_x (∂_y f)
 *
 * Complexity: O(1).
 */
double nr_fd_deriv2_mixed(const nr_gf_t *gf, int i, int j, int k,
                           double dx, double dy,
                           int dir1, int dir2, int order);

/**
 * L5: Compute the flat-space Laplacian Δ f = ∂²_x f + ∂²_y f + ∂²_z f.
 *
 * Complexity: O(1).
 */
double nr_fd_laplacian(const nr_gf_t *gf, int i, int j, int k,
                        double dx, double dy, double dz, int order);

/**
 * L5: Apply upwind (one-sided) derivative for advection terms.
 *
 * Upwind differencing is needed for the β^k ∂_k terms in BSSN:
 *   β^x ∂_x f ≈ β^x (f_{i+1} − f_{i-1})/(2dx) − |β^x| (f_{i+1} − 2f_i + f_{i-1})/(2dx)
 *
 * This adds numerical viscosity in the upwind direction.
 *
 * @param gf     Grid function
 * @param i,j,k  Grid indices
 * @param dx     Grid spacing in x
 * @param beta   Advection speed β
 * @param dir    Direction (0=x, 1=y, 2=z)
 * @param order  FD order for the centered part (upwind correction is always 2nd-order)
 * @return       Upwind-biased first derivative
 *
 * Complexity: O(1).
 */
double nr_fd_upwind(const nr_gf_t *gf, int i, int j, int k,
                     double dx, double beta, int dir, int order);

/* ===========================================================================
 * L5: Runge-Kutta Time Integration
 * =========================================================================== */

/**
 * L5: RK4 integration step for a single grid function.
 *
 * Advances f by one time step dt using the classic 4th-order Runge-Kutta:
 *   k1 = dt * RHS(f)
 *   k2 = dt * RHS(f + k1/2)
 *   k3 = dt * RHS(f + k2/2)
 *   k4 = dt * RHS(f + k3)
 *   f_new = f + (k1 + 2k2 + 2k3 + k4) / 6
 *
 * @param f         Current values (in)
 * @param rhs       RHS function: computes ∂_t f given current f
 * @param dt        Time step
 * @param f_new     Updated values (out)
 * @param ctx       User context passed to rhs
 * @param N         Number of grid points per function
 *
 * RHS callback signature:
 *   void rhs_func(const double *f, double *dfdt, void *ctx);
 *
 * Complexity: O(N) per substep, total 4 × O(N).
 * Reference: Butcher, "Numerical Methods for Ordinary Differential Equations" (2008)
 */
typedef void (*nr_rhs_func_t)(const double *f, double *dfdt, void *ctx);

void nr_rk4_step(const double *f, nr_rhs_func_t rhs, double dt,
                  double *f_new, void *ctx, int N);

/**
 * L5: RK2 (midpoint) integration for simpler problems.
 *
 * k1 = dt * RHS(f)
 * k2 = dt * RHS(f + k1/2)
 * f_new = f + k2
 *
 * Complexity: O(N).
 */
void nr_rk2_step(const double *f, nr_rhs_func_t rhs, double dt,
                  double *f_new, void *ctx, int N);

/**
 * L5: Low-storage 3rd-order Runge-Kutta (Williamson 1980).
 *
 * Uses only 2 storage registers instead of 4 for RK4, advantageous
 * for large grids in numerical relativity.
 *
 * Reference: Williamson, J. Comp. Phys. 35, 48 (1980)
 * Complexity: O(N).
 */
void nr_rk3_low_storage_step(const double *f, nr_rhs_func_t rhs,
                              double dt, double *f_new,
                              void *ctx, int N);

/* ===========================================================================
 * L5: Adaptive Time Stepping (CFL Condition)
 * =========================================================================== */

/**
 * Compute the CFL-limited time step for explicit integration.
 *
 * For the BSSN system, the Courant condition is:
 *   Δt ≤ CFL * min(dx, dy, dz) / c_max
 *
 * where c_max is the maximum characteristic speed (lapse + |shift| + gauge speed).
 * Typical CFL = 0.25–0.5 for RK4.
 *
 * @param alpha_max   Maximum lapse on grid
 * @param beta_max    Maximum shift magnitude on grid
 * @param dx, dy, dz  Grid spacings
 * @param cfl_factor  CFL safety factor (e.g., 0.25)
 * @return            Maximum stable time step
 *
 * Complexity: O(1).
 * Reference: Baumgarte & Shapiro (2010), Sec.9.3
 */
double nr_cfl_timestep(double alpha_max, double beta_max,
                        double dx, double dy, double dz,
                        double cfl_factor);

/* ===========================================================================
 * L5: Trilinear Interpolation
 *
 * Interpolate grid function values to arbitrary positions within the
 * computational domain. Essential for wave extraction on spheres,
 * apparent horizon evaluation, etc.
 * =========================================================================== */

/**
 * Trilinear interpolation of a grid function at (x, y, z).
 *
 * Finds the enclosing cell and interpolates using the 8 corner values.
 *
 * @param gf     Grid function
 * @param grid   Grid structure (for coordinate mapping)
 * @param x,y,z  Physical coordinates to interpolate at
 * @return       Interpolated value, or 0.0 if outside domain
 *
 * Complexity: O(1).
 */
double nr_interp_trilinear(const nr_gf_t *gf, const nr_grid_t *grid,
                            double x, double y, double z);

/**
 * Trilinear interpolation of all 6 symmetric tensor components at (x, y, z).
 *
 * @param tensor  Symmetric 3-tensor
 * @param grid    Grid structure
 * @param x,y,z   Physical coordinates
 * @param result[6]  Output tensor components {xx, xy, xz, yy, yz, zz}
 *
 * Complexity: O(1) per component.
 */
void nr_interp_tensor_trilinear(const nr_sym_tensor3_t *tensor,
                                 const nr_grid_t *grid,
                                 double x, double y, double z,
                                 double result[6]);

/**
 * Trilinear interpolation of a 3-vector at (x, y, z).
 *
 * @param vec     Spatial vector
 * @param grid    Grid structure
 * @param x,y,z   Physical coordinates
 * @param result[3]  Output vector components
 *
 * Complexity: O(1) per component.
 */
void nr_interp_vector_trilinear(const nr_vector3_t *vec,
                                 const nr_grid_t *grid,
                                 double x, double y, double z,
                                 double result[3]);

/* ===========================================================================
 * L3: Special Functions
 * =========================================================================== */

/**
 * Associated Legendre polynomial P_l^m(x) for −1 ≤ x ≤ 1.
 *
 * Computed via stable recurrence:
 *   (l−m) P_l^m = x(2l−1) P_{l−1}^m − (l+m−1) P_{l−2}^m
 *
 * @param l  Degree (l ≥ 0)
 * @param m  Order (0 ≤ m ≤ l)
 * @param x  Argument in [−1, 1]
 * @return   P_l^m(x)
 *
 * Reference: Press et al., "Numerical Recipes" (2007), Sec. 6.8
 * Complexity: O(l).
 */
double nr_legendre_P(int l, int m, double x);

/**
 * Spherical harmonic Y_{lm}(θ, φ) — scalar (spin-weight 0).
 *
 * Y_{lm}(θ, φ) = N_{lm} P_l^m(cos θ) e^{i m φ}
 * N_{lm} = sqrt((2l+1)/(4π) * (l−m)!/(l+m)!)
 *
 * Complexity: O(l).
 */
nr_complex_t nr_spherical_harmonic(int l, int m, double theta, double phi);

/**
 * Factorial function n! for small n (up to 20).
 * Complexity: O(n).
 */
double nr_factorial(int n);

/**
 * Double factorial n!! = n·(n−2)·(n−4)···.
 * Complexity: O(n).
 */
double nr_double_factorial(int n);

/* ===========================================================================
 * L5: SOR (Successive Over-Relaxation) Solver
 *
 * Used for solving the elliptic Hamiltonian constraint equation.
 * =========================================================================== */

/**
 * One SOR iteration for the 7-point Laplacian on a 3D grid.
 *
 * Updates u to better satisfy:
 *   Δ u = source  (in 3D with uniform spacing h)
 *
 * The 7-point stencil (in 3D):
 *   u_{i,j,k}^{new} = (1−ω) u_{i,j,k}^{old} + (ω/6) (u_{i+1} + u_{i−1}
 *                     + u_{j+1} + u_{j−1} + u_{k+1} + u_{k−1} − h² source_{i,j,k})
 *
 * @param u        Current solution (modified in place)
 * @param source   Source term (on all grid points)
 * @param grid     Grid structure
 * @param h        Grid spacing (assumed uniform)
 * @param omega    Over-relaxation parameter (1.0 = Gauss-Seidel, 1.0 < ω < 2.0)
 * @return         Maximum change |u_new − u_old|
 *
 * Complexity: O(N^3).
 * Reference: Press et al. (2007), Sec. 20.5
 */
double nr_sor_iteration_3d(nr_gf_t *u, const nr_gf_t *source,
                            const nr_grid_t *grid, double h,
                            double omega);

/**
 * Solve Δ u = source using SOR iteration to tolerance.
 *
 * @param u        Initial guess + output solution
 * @param source   Source term
 * @param grid     Grid
 * @param h        Spacing
 * @param omega    SOR parameter
 * @param tol      Convergence tolerance (max change)
 * @param max_iter Maximum iterations
 * @return         Number of iterations, or -1 if not converged
 *
 * Complexity: O(N^3 · iter).
 */
int nr_sor_solve_3d(nr_gf_t *u, const nr_gf_t *source,
                     const nr_grid_t *grid, double h, double omega,
                     double tol, int max_iter);

/* ===========================================================================
 * L5: Kreiss-Oliger Dissipation (Generic)
 * =========================================================================== */

/**
 * Apply Kreiss-Oliger dissipation to a single grid function.
 *
 * For order 2p:
 *   f_new = f + ε (−1)^p h^{2p−1} D_+^p D_-^p f / 2^{2p}
 *
 * @param gf     Grid function (modified in place)
 * @param eps    Dissipation strength (typ. 0.1–0.5)
 * @param order  Dissipation order (even ≥ 2, typ. 6)
 * @param dx     Grid spacing in x
 * @param dy     Grid spacing in y (same as dx for uniform grids)
 * @param dz     Grid spacing in z
 *
 * Complexity: O(N^3).
 * Reference: Kreiss & Oliger, GARP Pub. Ser. 10 (1973)
 */
void nr_kreiss_oliger_gf(nr_gf_t *gf, double eps, int order,
                          double dx, double dy, double dz);

#ifdef __cplusplus
}
#endif

#endif /* NR_UTILS_H */
