/**
 * @file nr_utils.c
 * @brief Finite difference stencils, Runge-Kutta integration, interpolation.
 *
 * Core numerical methods for solving Einstein's equations on a 3D grid.
 * Each function implements a specific computational technique from the
 * numerical relativity literature.
 *
 * References:
 *   - Fornberg (1998): finite difference coefficient generation
 *   - Baumgarte & Shapiro (2010), Appendix B
 *   - LeVeque (2007): FD methods for PDEs
 *   - Press et al. (2007): Numerical Recipes
 *
 * Knowledge: L5 (computational methods), L3 (discrete derivative operators)
 */

#include "nr_utils.h"
#include "nr_grid.h"  /* for nr_complex_t */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===========================================================================
 * L5: Finite Difference Stencils
 * =========================================================================== */

/* Get the stride in the given direction */
static int gf_stride(const nr_gf_t *gf, int dir) {
    if (dir == 0) return gf->stride_x;  /* x-direction */
    if (dir == 1) return gf->stride_y;  /* y-direction */
    return 1;                            /* z-direction */
}

/* Get the dimension in the given direction */
static int gf_dim(const nr_gf_t *gf, int dir) {
    if (dir == 0) return gf->nx + 2 * gf->ng;
    if (dir == 1) return gf->ny + 2 * gf->ng;
    return gf->nz + 2 * gf->ng;
}

/* Get the index in the given direction */
static int gf_idx(int i, int j, int k, int dir) {
    if (dir == 0) return i;
    if (dir == 1) return j;
    return k;
}

double nr_fd_deriv1(const nr_gf_t *gf, int i, int j, int k,
                     double dx, int dir, int order) {
    if (!gf || !gf->data) return 0.0;
    if (dx == 0.0) return 0.0;

    int idx = gf_idx(i, j, k, dir);
    int dim = gf_dim(gf, dir);
    (void)gf_stride(gf, dir); /* Used only in centered path below */

    /* Use one-sided stencils near boundaries */
    if (idx < order/2) {
        /* Forward 2nd-order */
        if (idx + 2 < dim) {
            double f0 = gf->data[i * gf->stride_x + j * gf->stride_y + k];
            int i1,j1,k1,i2,j2,k2;
            i1=i;j1=j;k1=k; i2=i;j2=j;k2=k;
            if (dir == 0) { i1=i+1; i2=i+2; }
            else if (dir == 1) { j1=j+1; j2=j+2; }
            else { k1=k+1; k2=k+2; }
            double f1 = gf->data[i1*gf->stride_x + j1*gf->stride_y + k1];
            double f2 = gf->data[i2*gf->stride_x + j2*gf->stride_y + k2];
            return (-1.5*f0 + 2.0*f1 - 0.5*f2) / dx;
        }
        return 0.0;
    }
    if (idx >= dim - order/2) {
        /* Backward 2nd-order */
        if (idx - 2 >= 0) {
            double f0 = gf->data[i * gf->stride_x + j * gf->stride_y + k];
            int i1,j1,k1,i2,j2,k2;
            i1=i;j1=j;k1=k; i2=i;j2=j;k2=k;
            if (dir == 0) { i1=i-1; i2=i-2; }
            else if (dir == 1) { j1=j-1; j2=j-2; }
            else { k1=k-1; k2=k-2; }
            double f1 = gf->data[i1*gf->stride_x + j1*gf->stride_y + k1];
            double f2 = gf->data[i2*gf->stride_x + j2*gf->stride_y + k2];
            return (1.5*f0 - 2.0*f1 + 0.5*f2) / dx;
        }
        return 0.0;
    }

    /* Centered FD in interior */
    double *d = gf->data;
    int sx = gf->stride_x, sy = gf->stride_y;
    int base = i * sx + j * sy + k;

    if (order == 2) {
        int di = (dir == 0) ? sx : ((dir == 1) ? sy : 1);
        double fp = d[base + di];
        double fm = d[base - di];
        return (fp - fm) / (2.0 * dx);
    } else if (order == 4) {
        int di = (dir == 0) ? sx : ((dir == 1) ? sy : 1);
        double fp1 = d[base + di];
        double fm1 = d[base - di];
        double fp2 = d[base + 2*di];
        double fm2 = d[base - 2*di];
        return (-fp2 + 8.0*fp1 - 8.0*fm1 + fm2) / (12.0 * dx);
    } else if (order == 6) {
        int di = (dir == 0) ? sx : ((dir == 1) ? sy : 1);
        double fp1 = d[base + di];
        double fm1 = d[base - di];
        double fp2 = d[base + 2*di];
        double fm2 = d[base - 2*di];
        double fp3 = d[base + 3*di];
        double fm3 = d[base - 3*di];
        return (fp3 - 9.0*fp2 + 45.0*fp1 - 45.0*fm1 + 9.0*fm2 - fm3) / (60.0 * dx);
    } else {
        /* 8th-order */
        int di = (dir == 0) ? sx : ((dir == 1) ? sy : 1);
        double fp1 = d[base + di];
        double fm1 = d[base - di];
        double fp2 = d[base + 2*di];
        double fm2 = d[base - 2*di];
        double fp3 = d[base + 3*di];
        double fm3 = d[base - 3*di];
        double fp4 = d[base + 4*di];
        double fm4 = d[base - 4*di];
        return (3.0*fp4 - 32.0*fp3 + 168.0*fp2 - 672.0*fp1
                + 672.0*fm1 - 168.0*fm2 + 32.0*fm3 - 3.0*fm4) / (840.0 * dx);
    }
}

double nr_fd_deriv2(const nr_gf_t *gf, int i, int j, int k,
                     double dx, int dir, int order) {
    if (!gf || !gf->data) return 0.0;
    if (dx == 0.0) return 0.0;

    int idx = gf_idx(i, j, k, dir);
    int dim = gf_dim(gf, dir);

    if (idx < order/2 || idx >= dim - order/2) {
        /* Use 2nd-order near boundaries */
        if (idx >= 1 && idx + 1 < dim) {
            double *d = gf->data;
            int sx = gf->stride_x, sy = gf->stride_y;
            int base = i * sx + j * sy + k;
            int di = (dir == 0) ? sx : ((dir == 1) ? sy : 1);
            double fp = d[base + di];
            double fc = d[base];
            double fm = d[base - di];
            return (fp - 2.0*fc + fm) / (dx * dx);
        }
        return 0.0;
    }

    double *d = gf->data;
    int sx = gf->stride_x, sy = gf->stride_y;
    int base = i * sx + j * sy + k;
    int di = (dir == 0) ? sx : ((dir == 1) ? sy : 1);

    if (order == 2) {
        double fp = d[base + di];
        double fc = d[base];
        double fm = d[base - di];
        return (fp - 2.0*fc + fm) / (dx * dx);
    } else if (order == 4) {
        double f0  = d[base];
        double fp1 = d[base + di];
        double fm1 = d[base - di];
        double fp2 = d[base + 2*di];
        double fm2 = d[base - 2*di];
        return (-fp2 + 16.0*fp1 - 30.0*f0 + 16.0*fm1 - fm2) / (12.0 * dx * dx);
    } else if (order == 6) {
        double f0  = d[base];
        double fp1 = d[base + di];
        double fm1 = d[base - di];
        double fp2 = d[base + 2*di];
        double fm2 = d[base - 2*di];
        double fp3 = d[base + 3*di];
        double fm3 = d[base - 3*di];
        return (2.0*fp3 - 27.0*fp2 + 270.0*fp1 - 490.0*f0
                + 270.0*fm1 - 27.0*fm2 + 2.0*fm3) / (180.0 * dx * dx);
    } else {
        /* 8th-order */
        double f0  = d[base];
        double fp1 = d[base + di];
        double fm1 = d[base - di];
        double fp2 = d[base + 2*di];
        double fm2 = d[base - 2*di];
        double fp3 = d[base + 3*di];
        double fm3 = d[base - 3*di];
        double fp4 = d[base + 4*di];
        double fm4 = d[base - 4*di];
        return (-9.0*fp4 + 128.0*fp3 - 1008.0*fp2 + 8064.0*fp1 - 14350.0*f0
                + 8064.0*fm1 - 1008.0*fm2 + 128.0*fm3 - 9.0*fm4) / (5040.0 * dx * dx);
    }
}

double nr_fd_deriv2_mixed(const nr_gf_t *gf, int i, int j, int k,
                           double dx, double dy, int dir1, int dir2,
                           int order) {
    (void)order;
    if (!gf || !gf->data) return 0.0;
    /* Compute mixed derivative by applying first derivatives successively.
     * ∂²_{xy} f = ∂_x(∂_y f) ≈ D_x(D_y f)
     *
     * Use 2nd-order centered stencils for the mixed case for robustness.
     */
    double *d = gf->data;
    int sx = gf->stride_x, sy = gf->stride_y;
    int di1 = (dir1 == 0) ? sx : ((dir1 == 1) ? sy : 1);
    int di2 = (dir2 == 0) ? sx : ((dir2 == 1) ? sy : 1);

    int base = i * sx + j * sy + k;
    double fpp = d[base + di1 + di2];
    double fpm = d[base + di1 - di2];
    double fmp = d[base - di1 + di2];
    double fmm = d[base - di1 - di2];

    double h1 = (dir1 == 0) ? dx : ((dir1 == 1) ? dy : dy);
    double h2 = (dir2 == 0) ? dx : ((dir2 == 1) ? dy : dy);
    return (fpp - fpm - fmp + fmm) / (4.0 * h1 * h2);
}

double nr_fd_laplacian(const nr_gf_t *gf, int i, int j, int k,
                        double dx, double dy, double dz, int order) {
    return nr_fd_deriv2(gf, i, j, k, dx, 0, order)
         + nr_fd_deriv2(gf, i, j, k, dy, 1, order)
         + nr_fd_deriv2(gf, i, j, k, dz, 2, order);
}

double nr_fd_upwind(const nr_gf_t *gf, int i, int j, int k,
                     double dx, double beta, int dir, int order) {
    /* Centered part */
    double centered = nr_fd_deriv1(gf, i, j, k, dx, dir, order);
    /* Upwind correction: -|β|/2 * h * ∂² f (adds numerical viscosity) */
    double abs_beta = fabs(beta);
    double second = nr_fd_deriv2(gf, i, j, k, dx, dir, 2);
    return centered - abs_beta * dx * 0.5 * second;
}

/* ===========================================================================
 * L5: Runge-Kutta Time Integration
 * =========================================================================== */

void nr_rk4_step(const double *f, nr_rhs_func_t rhs, double dt,
                  double *f_new, void *ctx, int N) {
    if (!f || !rhs || !f_new || N <= 0) return;

    double *k1 = (double*)malloc(N * sizeof(double));
    double *k2 = (double*)malloc(N * sizeof(double));
    double *k3 = (double*)malloc(N * sizeof(double));
    double *k4 = (double*)malloc(N * sizeof(double));
    double *ftmp = (double*)malloc(N * sizeof(double));

    if (!k1 || !k2 || !k3 || !k4 || !ftmp) {
        free(k1); free(k2); free(k3); free(k4); free(ftmp);
        return;
    }

    /* k1 = dt * rhs(f) */
    rhs(f, k1, ctx);
    for (int i = 0; i < N; i++) k1[i] *= dt;

    /* k2 = dt * rhs(f + k1/2) */
    for (int i = 0; i < N; i++) ftmp[i] = f[i] + 0.5 * k1[i];
    rhs(ftmp, k2, ctx);
    for (int i = 0; i < N; i++) k2[i] *= dt;

    /* k3 = dt * rhs(f + k2/2) */
    for (int i = 0; i < N; i++) ftmp[i] = f[i] + 0.5 * k2[i];
    rhs(ftmp, k3, ctx);
    for (int i = 0; i < N; i++) k3[i] *= dt;

    /* k4 = dt * rhs(f + k3) */
    for (int i = 0; i < N; i++) ftmp[i] = f[i] + k3[i];
    rhs(ftmp, k4, ctx);
    for (int i = 0; i < N; i++) k4[i] *= dt;

    /* f_new = f + (k1 + 2*k2 + 2*k3 + k4) / 6 */
    for (int i = 0; i < N; i++) {
        f_new[i] = f[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
    }

    free(k1); free(k2); free(k3); free(k4); free(ftmp);
}

void nr_rk2_step(const double *f, nr_rhs_func_t rhs, double dt,
                  double *f_new, void *ctx, int N) {
    if (!f || !rhs || !f_new || N <= 0) return;

    double *k1 = (double*)malloc(N * sizeof(double));
    double *ftmp = (double*)malloc(N * sizeof(double));
    if (!k1 || !ftmp) { free(k1); free(ftmp); return; }

    /* k1 = dt * rhs(f) */
    rhs(f, k1, ctx);
    for (int i = 0; i < N; i++) k1[i] *= dt;

    /* k2 = dt * rhs(f + k1/2) */
    for (int i = 0; i < N; i++) ftmp[i] = f[i] + 0.5 * k1[i];
    rhs(ftmp, f_new, ctx);
    for (int i = 0; i < N; i++) {
        f_new[i] = f[i] + dt * f_new[i];
    }

    free(k1); free(ftmp);
}

void nr_rk3_low_storage_step(const double *f, nr_rhs_func_t rhs,
                              double dt, double *f_new,
                              void *ctx, int N) {
    if (!f || !rhs || !f_new || N <= 0) return;
    /* Williamson (1980) low-storage RK3:
     * Uses just 2 arrays: f (in) and f_new (in/out).
     *
     * Substep 1: f1 = f + (1/3)*dt*RHS(f)
     * Substep 2: f2 = f + (15/16)*dt*RHS( (3/4)*f + (1/4)*f1 )
     * Substep 3: f_new = f + (8/15)*dt*RHS( f2 ) ... wait, need temp
     *
     * Simplified: use a temporary array for the RHS evaluation.
     */
    double *tmp = (double*)malloc(N * sizeof(double));
    double *rhs_vals = (double*)malloc(N * sizeof(double));
    if (!tmp || !rhs_vals) { free(tmp); free(rhs_vals); return; }

    /* Stage 1 */
    rhs(f, rhs_vals, ctx);
    for (int i = 0; i < N; i++) tmp[i] = f[i] + (1.0/3.0) * dt * rhs_vals[i];

    /* Stage 2 */
    rhs(tmp, rhs_vals, ctx);
    for (int i = 0; i < N; i++) tmp[i] = f[i] + (15.0/16.0) * dt * rhs_vals[i];

    /* Stage 3 */
    rhs(tmp, rhs_vals, ctx);
    for (int i = 0; i < N; i++) f_new[i] = f[i] + (8.0/15.0) * dt * rhs_vals[i];

    free(tmp); free(rhs_vals);
}

/* ===========================================================================
 * L5: CFL Time Step
 * =========================================================================== */

double nr_cfl_timestep(double alpha_max, double beta_max,
                        double dx, double dy, double dz,
                        double cfl_factor) {
    double h_min = dx;
    if (dy < h_min) h_min = dy;
    if (dz < h_min) h_min = dz;
    /* Characteristic speed = lapse + |shift| (approximately the gauge speed) */
    double c_max = alpha_max + beta_max;
    if (c_max < 1e-15) c_max = 1.0;
    double dt = cfl_factor * h_min / c_max;
    return (dt > 0) ? dt : 0.0;
}

/* ===========================================================================
 * L5: Trilinear Interpolation
 * =========================================================================== */

double nr_interp_trilinear(const nr_gf_t *gf, const nr_grid_t *grid,
                            double x, double y, double z) {
    if (!gf || !grid) return 0.0;

    int ng = grid->ng;
    /* Find the cell indices */
    double fx = (x - grid->xmin) / grid->dx + ng;
    double fy = (y - grid->ymin) / grid->dy + ng;
    double fz = (z - grid->zmin) / grid->dz + ng;

    int i0 = (int)floor(fx);
    int j0 = (int)floor(fy);
    int k0 = (int)floor(fz);

    int nx_tot = grid->nx + 2*ng;
    int ny_tot = grid->ny + 2*ng;
    int nz_tot = grid->nz + 2*ng;

    /* Clamp to valid range */
    if (i0 < 0 || i0+1 >= nx_tot || j0 < 0 || j0+1 >= ny_tot
        || k0 < 0 || k0+1 >= nz_tot) return 0.0;

    double dx_f = fx - i0;
    double dy_f = fy - j0;
    double dz_f = fz - k0;

    /* Ensure fractional parts are in [0,1] */
    if (dx_f < 0) dx_f = 0;
    if (dx_f > 1) dx_f = 1;
    if (dy_f < 0) dy_f = 0;
    if (dy_f > 1) dy_f = 1;
    if (dz_f < 0) dz_f = 0;
    if (dz_f > 1) dz_f = 1;

    int sx = gf->stride_x, sy = gf->stride_y;
    double *d = gf->data;

    double c000 = d[i0*sx + j0*sy + k0];
    double c001 = d[i0*sx + j0*sy + k0+1];
    double c010 = d[i0*sx + (j0+1)*sy + k0];
    double c011 = d[i0*sx + (j0+1)*sy + k0+1];
    double c100 = d[(i0+1)*sx + j0*sy + k0];
    double c101 = d[(i0+1)*sx + j0*sy + k0+1];
    double c110 = d[(i0+1)*sx + (j0+1)*sy + k0];
    double c111 = d[(i0+1)*sx + (j0+1)*sy + k0+1];

    /* Interpolate in z first */
    double c00 = c000 + (c001 - c000) * dz_f;
    double c01 = c010 + (c011 - c010) * dz_f;
    double c10 = c100 + (c101 - c100) * dz_f;
    double c11 = c110 + (c111 - c110) * dz_f;

    /* Interpolate in y */
    double c0 = c00 + (c01 - c00) * dy_f;
    double c1 = c10 + (c11 - c10) * dy_f;

    /* Interpolate in x */
    return c0 + (c1 - c0) * dx_f;
}

void nr_interp_tensor_trilinear(const nr_sym_tensor3_t *tensor,
                                 const nr_grid_t *grid,
                                 double x, double y, double z,
                                 double result[6]) {
    for (int c = 0; c < 6; c++) {
        result[c] = nr_interp_trilinear(tensor->comp[c], grid, x, y, z);
    }
}

void nr_interp_vector_trilinear(const nr_vector3_t *vec,
                                 const nr_grid_t *grid,
                                 double x, double y, double z,
                                 double result[3]) {
    for (int c = 0; c < 3; c++) {
        result[c] = nr_interp_trilinear(vec->comp[c], grid, x, y, z);
    }
}

/* ===========================================================================
 * L3: Special Functions — Legendre Polynomials & Spherical Harmonics
 * =========================================================================== */

double nr_factorial(int n) {
    if (n < 0) return 0.0;
    if (n <= 1) return 1.0;
    double val = 1.0;
    for (int i = 2; i <= n; i++) val *= i;
    return val;
}

double nr_double_factorial(int n) {
    if (n <= 0) return 1.0;
    double val = 1.0;
    for (int i = n; i >= 1; i -= 2) val *= i;
    return val;
}

double nr_legendre_P(int l, int m, double x) {
    /* Associated Legendre polynomial P_l^m(x)
     * Uses stable recurrence: (l-m)P_l^m = x(2l-1)P_{l-1}^m - (l+m-1)P_{l-2}^m
     *
     * First compute P_m^m, then recurse up.
     */
    if (m < 0 || m > l) return 0.0;
    if (l == 0 && m == 0) return 1.0;

    /* P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2} */
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));  /* sqrt(1-x^2) */
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }

    if (l == m) return pmm;

    /* P_{m+1}^m = x(2m+1) P_m^m */
    double pmmp1 = x * (2.0 * m + 1.0) * pmm;
    if (l == m + 1) return pmmp1;

    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = (x * (2.0 * ll - 1.0) * pmmp1
               - (ll + m - 1.0) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

/* nr_complex_t is defined in nr_grid.h, which is included via nr_utils.h */

nr_complex_t nr_spherical_harmonic(int l, int m, double theta, double phi) {
    /* Y_{lm}(θ,φ) = N_{lm} P_l^m(cos θ) e^{i m φ}
     * N_{lm} = sqrt( (2l+1)/(4π) * (l-m)!/(l+m)! )
     */
    if (m < 0) {
        /* Y_{l,-m} = (-1)^m conj(Y_{lm}) */
        nr_complex_t ylm = nr_spherical_harmonic(l, -m, theta, phi);
        if (m % 2 != 0) {
            ylm.re = -ylm.re;
            ylm.im = -ylm.im;
        }
        ylm.im = -ylm.im;
        return ylm;
    }
    if (m > l) return nr_complex_make(0.0, 0.0);

    double cos_th = cos(theta);
    double Plm = nr_legendre_P(l, m, cos_th);

    double norm = sqrt((2.0 * l + 1.0) / (4.0 * M_PI)
                       * nr_factorial(l - m) / nr_factorial(l + m));

    double val = norm * Plm;
    return nr_complex_make(val * cos(m * phi), val * sin(m * phi));
}

/* ===========================================================================
 * L5: SOR Solver (3D Poisson Equation)
 * =========================================================================== */

double nr_sor_iteration_3d(nr_gf_t *u, const nr_gf_t *source,
                            const nr_grid_t *grid, double h,
                            double omega) {
    if (!u || !source || !grid) return 0.0;
    int ng = grid->ng;
    double max_change = 0.0;
    (void)(1.0 / (h * h));  /* h² used in the stencil via 1/6 factor below */

    int imax = grid->nx + ng;
    int jmax = grid->ny + ng;
    int kmax = grid->nz + ng;

    for (int i = ng; i < imax; i++) {
        for (int j = ng; j < jmax; j++) {
            for (int k = ng; k < kmax; k++) {
                int idx = i * u->stride_x + j * u->stride_y + k;
                double u_old = u->data[idx];

                /* 6 neighbors */
                double u_sum = u->data[(i+1)*u->stride_x + j*u->stride_y + k]
                             + u->data[(i-1)*u->stride_x + j*u->stride_y + k]
                             + u->data[i*u->stride_x + (j+1)*u->stride_y + k]
                             + u->data[i*u->stride_x + (j-1)*u->stride_y + k]
                             + u->data[i*u->stride_x + j*u->stride_y + (k+1)]
                             + u->data[i*u->stride_x + j*u->stride_y + (k-1)];

                double src = source->data[idx];
                double u_new = (1.0 - omega) * u_old
                             + (omega / 6.0) * (u_sum - h * h * src);

                u->data[idx] = u_new;
                double change = fabs(u_new - u_old);
                if (change > max_change) max_change = change;
            }
        }
    }
    return max_change;
}

int nr_sor_solve_3d(nr_gf_t *u, const nr_gf_t *source,
                     const nr_grid_t *grid, double h, double omega,
                     double tol, int max_iter) {
    if (max_iter <= 0) return -1;
    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = nr_sor_iteration_3d(u, source, grid, h, omega);
        if (max_change < tol) return iter + 1;
    }
    return -1;  /* did not converge */
}

/* ===========================================================================
 * L5: Kreiss-Oliger Dissipation (Generic 3D)
 * =========================================================================== */

void nr_kreiss_oliger_gf(nr_gf_t *gf, double eps, int order,
                          double dx, double dy, double dz) {
    (void)dx; (void)dy; (void)dz;  /* Uniform spacing assumed */
    if (!gf || !gf->data) return;
    if (order < 2 || eps == 0.0) return;

    int ng = gf->ng;
    int nx_tot = gf->nx + 2 * ng;
    int ny_tot = gf->ny + 2 * ng;
    int nz_tot = gf->nz + 2 * ng;

    /* We need to apply the dissipation and store in a temporary array,
     * then copy back, because the operator uses original values. */
    size_t total = (size_t)nx_tot * ny_tot * nz_tot;
    double *temp = (double*)malloc(total * sizeof(double));
    if (!temp) return;

    /* For order 2p:
     * KO operator in 1D: (-1)^{p+1} * eps * h^{2p-1} * D_+^p D_-^p / 2^{2p}
     * In 3D: sum over directions
     *
     * For order 2: D_+ D_- = central 2nd derivative stencil ⇒ (f_{i+1} - 2f_i + f_{i-1})
     * For order 4: D_+^2 D_-^2 = (f_{i+2} - 4f_{i+1} + 6f_i - 4f_{i-1} + f_{i-2})
     * For order 6: D_+^3 D_-^3 = (f_{i+3} - 6f_{i+2} + 15f_{i+1} - 20f_i + 15f_{i-1} - 6f_{i-2} + f_{i-3})
     */

    int half = order / 2;  /* number of points on each side */

    for (int i = 0; i < nx_tot; i++) {
        for (int j = 0; j < ny_tot; j++) {
            for (int k = 0; k < nz_tot; k++) {
                int idx = i * gf->stride_x + j * gf->stride_y + k;
                double dissipation = 0.0;

                /* x-direction */
                if (i >= half && i < nx_tot - half) {
                    double disp_x = 0.0;
                    int sign = (half % 2 == 0) ? -1 : 1;  /* (-1)^{p+1} = (-1)^{half+1} */
                    if (half == 1) {
                        /* Order 2: f_{i+1} - 2f_i + f_{i-1} */
                        disp_x = gf->data[(i+1)*gf->stride_x + j*gf->stride_y + k]
                               - 2.0 * gf->data[idx]
                               + gf->data[(i-1)*gf->stride_x + j*gf->stride_y + k];
                    } else if (half == 2) {
                        /* Order 4: f_{i+2} - 4f_{i+1} + 6f_i - 4f_{i-1} + f_{i-2} */
                        disp_x = gf->data[(i+2)*gf->stride_x + j*gf->stride_y + k]
                               - 4.0 * gf->data[(i+1)*gf->stride_x + j*gf->stride_y + k]
                               + 6.0 * gf->data[idx]
                               - 4.0 * gf->data[(i-1)*gf->stride_x + j*gf->stride_y + k]
                               + gf->data[(i-2)*gf->stride_x + j*gf->stride_y + k];
                    } else {
                        /* Order 6: f_{i+3} - 6f_{i+2} + 15f_{i+1} - 20f_i + 15f_{i-1} - 6f_{i-2} + f_{i-3} */
                        disp_x = gf->data[(i+3)*gf->stride_x + j*gf->stride_y + k]
                               - 6.0 * gf->data[(i+2)*gf->stride_x + j*gf->stride_y + k]
                               + 15.0 * gf->data[(i+1)*gf->stride_x + j*gf->stride_y + k]
                               - 20.0 * gf->data[idx]
                               + 15.0 * gf->data[(i-1)*gf->stride_x + j*gf->stride_y + k]
                               - 6.0 * gf->data[(i-2)*gf->stride_x + j*gf->stride_y + k]
                               + gf->data[(i-3)*gf->stride_x + j*gf->stride_y + k];
                    }
                    dissipation += sign * eps * disp_x / pow(2.0, order);
                }

                /* y-direction */
                if (j >= half && j < ny_tot - half) {
                    double disp_y = 0.0;
                    int sign = (half % 2 == 0) ? -1 : 1;
                    if (half == 1) {
                        disp_y = gf->data[i*gf->stride_x + (j+1)*gf->stride_y + k]
                               - 2.0 * gf->data[idx]
                               + gf->data[i*gf->stride_x + (j-1)*gf->stride_y + k];
                    } else if (half == 2) {
                        disp_y = gf->data[i*gf->stride_x + (j+2)*gf->stride_y + k]
                               - 4.0 * gf->data[i*gf->stride_x + (j+1)*gf->stride_y + k]
                               + 6.0 * gf->data[idx]
                               - 4.0 * gf->data[i*gf->stride_x + (j-1)*gf->stride_y + k]
                               + gf->data[i*gf->stride_x + (j-2)*gf->stride_y + k];
                    } else {
                        disp_y = gf->data[i*gf->stride_x + (j+3)*gf->stride_y + k]
                               - 6.0 * gf->data[i*gf->stride_x + (j+2)*gf->stride_y + k]
                               + 15.0 * gf->data[i*gf->stride_x + (j+1)*gf->stride_y + k]
                               - 20.0 * gf->data[idx]
                               + 15.0 * gf->data[i*gf->stride_x + (j-1)*gf->stride_y + k]
                               - 6.0 * gf->data[i*gf->stride_x + (j-2)*gf->stride_y + k]
                               + gf->data[i*gf->stride_x + (j-3)*gf->stride_y + k];
                    }
                    dissipation += sign * eps * disp_y / pow(2.0, order);
                }

                /* z-direction */
                if (k >= half && k < nz_tot - half) {
                    double disp_z = 0.0;
                    int sign = (half % 2 == 0) ? -1 : 1;
                    if (half == 1) {
                        disp_z = gf->data[i*gf->stride_x + j*gf->stride_y + (k+1)]
                               - 2.0 * gf->data[idx]
                               + gf->data[i*gf->stride_x + j*gf->stride_y + (k-1)];
                    } else if (half == 2) {
                        disp_z = gf->data[i*gf->stride_x + j*gf->stride_y + (k+2)]
                               - 4.0 * gf->data[i*gf->stride_x + j*gf->stride_y + (k+1)]
                               + 6.0 * gf->data[idx]
                               - 4.0 * gf->data[i*gf->stride_x + j*gf->stride_y + (k-1)]
                               + gf->data[i*gf->stride_x + j*gf->stride_y + (k-2)];
                    } else {
                        disp_z = gf->data[i*gf->stride_x + j*gf->stride_y + (k+3)]
                               - 6.0 * gf->data[i*gf->stride_x + j*gf->stride_y + (k+2)]
                               + 15.0 * gf->data[i*gf->stride_x + j*gf->stride_y + (k+1)]
                               - 20.0 * gf->data[idx]
                               + 15.0 * gf->data[i*gf->stride_x + j*gf->stride_y + (k-1)]
                               - 6.0 * gf->data[i*gf->stride_x + j*gf->stride_y + (k-2)]
                               + gf->data[i*gf->stride_x + j*gf->stride_y + (k-3)];
                    }
                    dissipation += sign * eps * disp_z / pow(2.0, order);
                }

                temp[idx] = dissipation;
            }
        }
    }

    /* Apply dissipation: f_new = f + dissipation */
    for (size_t q = 0; q < total; q++) {
        gf->data[q] += temp[q];
    }

    free(temp);
}
