/**
 * @file nr_grid.c
 * @brief 3D Cartesian grid memory management and coordinate utilities.
 *
 * Implements allocation, deallocation, access, and coordinate mapping for
 * the fundamental grid data structures used in numerical relativity.
 *
 * References:
 *   - Baumgarte & Shapiro (2010), Ch. 6
 *   - Alcubierre (2008), Ch. 4
 *
 * Knowledge: L1 (grid definitions), L5 (memory layout for FD methods)
 */

#include "nr_grid.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ===========================================================================
 * Grid Function Memory Management
 * =========================================================================== */

nr_gf_t* nr_gf_alloc(int nx, int ny, int nz, int ng) {
    nr_gf_t *gf = (nr_gf_t*)malloc(sizeof(nr_gf_t));
    if (!gf) return NULL;

    gf->nx = nx;
    gf->ny = ny;
    gf->nz = nz;
    gf->ng = ng;

    int nx_tot = nx + 2 * ng;
    int ny_tot = ny + 2 * ng;
    int nz_tot = nz + 2 * ng;

    gf->stride_y = nz_tot;
    gf->stride_x = ny_tot * nz_tot;

    size_t total = (size_t)nx_tot * ny_tot * nz_tot;
    gf->data = (double*)calloc(total, sizeof(double));
    if (!gf->data) {
        free(gf);
        return NULL;
    }

    return gf;
}

void nr_gf_free(nr_gf_t *gf) {
    if (gf) {
        free(gf->data);
        free(gf);
    }
}

void nr_gf_set_all(nr_gf_t *gf, double value) {
    if (!gf || !gf->data) return;
    int nx_tot = gf->nx + 2 * gf->ng;
    int ny_tot = gf->ny + 2 * gf->ng;
    int nz_tot = gf->nz + 2 * gf->ng;
    size_t total = (size_t)nx_tot * ny_tot * nz_tot;
    for (size_t idx = 0; idx < total; idx++) {
        gf->data[idx] = value;
    }
}

void nr_gf_copy(nr_gf_t *dst, const nr_gf_t *src) {
    if (!dst || !src || !dst->data || !src->data) return;
    int nx_tot = src->nx + 2 * src->ng;
    int ny_tot = src->ny + 2 * src->ng;
    int nz_tot = src->nz + 2 * src->ng;
    size_t total = (size_t)nx_tot * ny_tot * nz_tot;
    memcpy(dst->data, src->data, total * sizeof(double));
}

double nr_gf_l2norm(const nr_gf_t *gf) {
    if (!gf || !gf->data) return 0.0;
    double sum = 0.0;
    int ng = gf->ng;
    int count = 0;
    for (int i = ng; i < gf->nx + ng; i++) {
        for (int j = ng; j < gf->ny + ng; j++) {
            for (int k = ng; k < gf->nz + ng; k++) {
                double v = gf->data[i * gf->stride_x + j * gf->stride_y + k];
                sum += v * v;
                count++;
            }
        }
    }
    return (count > 0) ? sqrt(sum / count) : 0.0;
}

double nr_gf_maxabs(const nr_gf_t *gf) {
    if (!gf || !gf->data) return 0.0;
    double maxv = 0.0;
    int ng = gf->ng;
    for (int i = ng; i < gf->nx + ng; i++) {
        for (int j = ng; j < gf->ny + ng; j++) {
            for (int k = ng; k < gf->nz + ng; k++) {
                double v = fabs(gf->data[i * gf->stride_x + j * gf->stride_y + k]);
                if (v > maxv) maxv = v;
            }
        }
    }
    return maxv;
}

/* ===========================================================================
 * Grid Structure Management
 * =========================================================================== */

void nr_grid_init(nr_grid_t *grid, int nx, int ny, int nz,
                  double xmin, double xmax,
                  double ymin, double ymax,
                  double zmin, double zmax,
                  int ng) {
    if (!grid) return;
    grid->nx = nx; grid->ny = ny; grid->nz = nz;
    grid->xmin = xmin; grid->xmax = xmax;
    grid->ymin = ymin; grid->ymax = ymax;
    grid->zmin = zmin; grid->zmax = zmax;
    grid->ng = ng;

    grid->dx = (nx > 1) ? (xmax - xmin) / (nx - 1) : 1.0;
    grid->dy = (ny > 1) ? (ymax - ymin) / (ny - 1) : 1.0;
    grid->dz = (nz > 1) ? (zmax - zmin) / (nz - 1) : 1.0;

    int nx_tot = nx + 2 * ng;
    int ny_tot = ny + 2 * ng;
    int nz_tot = nz + 2 * ng;

    grid->x_coords = (double*)malloc(nx_tot * sizeof(double));
    grid->y_coords = (double*)malloc(ny_tot * sizeof(double));
    grid->z_coords = (double*)malloc(nz_tot * sizeof(double));

    for (int i = 0; i < nx_tot; i++) {
        grid->x_coords[i] = xmin + (i - ng) * grid->dx;
    }
    for (int j = 0; j < ny_tot; j++) {
        grid->y_coords[j] = ymin + (j - ng) * grid->dy;
    }
    for (int k = 0; k < nz_tot; k++) {
        grid->z_coords[k] = zmin + (k - ng) * grid->dz;
    }
}

void nr_grid_free(nr_grid_t *grid) {
    if (grid) {
        free(grid->x_coords); grid->x_coords = NULL;
        free(grid->y_coords); grid->y_coords = NULL;
        free(grid->z_coords); grid->z_coords = NULL;
    }
}

double nr_grid_x(const nr_grid_t *grid, int i) {
    return grid->x_coords[i];
}

double nr_grid_y(const nr_grid_t *grid, int j) {
    return grid->y_coords[j];
}

double nr_grid_z(const nr_grid_t *grid, int k) {
    return grid->z_coords[k];
}

double nr_grid_r(const nr_grid_t *grid, int i, int j, int k) {
    double x = grid->x_coords[i];
    double y = grid->y_coords[j];
    double z = grid->z_coords[k];
    return sqrt(x*x + y*y + z*z);
}

void nr_grid_spherical(const nr_grid_t *grid, int i, int j, int k,
                       double *r, double *sin_theta, double *cos_theta,
                       double *sin_phi, double *cos_phi) {
    double x = grid->x_coords[i];
    double y = grid->y_coords[j];
    double z = grid->z_coords[k];
    double rloc = sqrt(x*x + y*y + z*z);
    *r = rloc;
    if (rloc > 1e-15) {
        *cos_theta = z / rloc;
        double rho = sqrt(x*x + y*y);
        *sin_theta = rho / rloc;
        if (rho > 1e-15) {
            *cos_phi = x / rho;
            *sin_phi = y / rho;
        } else {
            *cos_phi = 1.0;
            *sin_phi = 0.0;
        }
    } else {
        *sin_theta = 0.0;
        *cos_theta = 1.0;
        *cos_phi = 1.0;
        *sin_phi = 0.0;
    }
}

/* ===========================================================================
 * Vector and Tensor Allocation
 * =========================================================================== */

nr_sym_tensor3_t* nr_sym_tensor3_alloc(int nx, int ny, int nz, int ng) {
    nr_sym_tensor3_t *t = (nr_sym_tensor3_t*)malloc(sizeof(nr_sym_tensor3_t));
    if (!t) return NULL;
    for (int c = 0; c < 6; c++) {
        t->comp[c] = nr_gf_alloc(nx, ny, nz, ng);
        if (!t->comp[c]) {
            for (int d = 0; d < c; d++) nr_gf_free(t->comp[d]);
            free(t);
            return NULL;
        }
    }
    return t;
}

void nr_sym_tensor3_free(nr_sym_tensor3_t *t) {
    if (t) {
        for (int c = 0; c < 6; c++) nr_gf_free(t->comp[c]);
        free(t);
    }
}

nr_vector3_t* nr_vector3_alloc(int nx, int ny, int nz, int ng) {
    nr_vector3_t *v = (nr_vector3_t*)malloc(sizeof(nr_vector3_t));
    if (!v) return NULL;
    for (int c = 0; c < 3; c++) {
        v->comp[c] = nr_gf_alloc(nx, ny, nz, ng);
        if (!v->comp[c]) {
            for (int d = 0; d < c; d++) nr_gf_free(v->comp[d]);
            free(v);
            return NULL;
        }
    }
    return v;
}

void nr_vector3_free(nr_vector3_t *v) {
    if (v) {
        for (int c = 0; c < 3; c++) nr_gf_free(v->comp[c]);
        free(v);
    }
}
