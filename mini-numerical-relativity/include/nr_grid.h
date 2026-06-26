/**
 * @file nr_grid.h
 * @brief 3D Cartesian grid infrastructure for numerical relativity
 *
 * Numerical relativity solves Einstein's field equations on a discretized
 * 3D spatial grid. This header defines the fundamental data structures:
 * 3D staggered grids, boundary zones, multi-component fields, and the
 * complex number type used throughout the codebase.
 *
 * References:
 *   - Baumgarte & Shapiro, "Numerical Relativity" (2010), Ch. 6
 *   - Alcubierre, "Introduction to 3+1 Numerical Relativity" (2008), Ch. 4
 *   - MIT 8.962 (General Relativity), Prof. Hughes
 *
 * Knowledge Coverage:
 *   L1: Grid point, ghost zone, staggered grid definitions
 *   L3: 3D Cartesian mesh as discrete manifold approximation
 *   L5: Memory layout, indexing schemes for FD methods
 */

#ifndef NR_GRID_H
#define NR_GRID_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L1: Complex Number Type
 *
 * Used for wave extraction quantities: Ψ_4, h_{lm}, spin-weighted harmonics.
 * Placed here (rather than in nr_wave.h) to avoid circular dependencies.
 * --------------------------------------------------------------------------- */

/**
 * Complex number type for Ψ_4, strain modes, and spherical harmonics.
 */
typedef struct {
    double re;  /**< Real part */
    double im;  /**< Imaginary part */
} nr_complex_t;

/** Create a complex number from real and imaginary parts. Complexity: O(1). */
static inline nr_complex_t nr_complex_make(double re, double im) {
    nr_complex_t z = {re, im};
    return z;
}

/** Complex addition: z = a + b. Complexity: O(1). */
static inline nr_complex_t nr_complex_add(nr_complex_t a, nr_complex_t b) {
    return nr_complex_make(a.re + b.re, a.im + b.im);
}

/** Complex multiplication: z = a * b. Complexity: O(1). */
static inline nr_complex_t nr_complex_mul(nr_complex_t a, nr_complex_t b) {
    return nr_complex_make(a.re * b.re - a.im * b.im,
                           a.re * b.im + a.im * b.re);
}

/** Complex norm squared: |z|^2. Complexity: O(1). */
static inline double nr_complex_norm2(nr_complex_t z) {
    return z.re * z.re + z.im * z.im;
}

/** Complex conjugate: z̄. Complexity: O(1). */
static inline nr_complex_t nr_complex_conj(nr_complex_t z) {
    return nr_complex_make(z.re, -z.im);
}

/* ---------------------------------------------------------------------------
 * L1: Grid Point & Ghost Zone Constants
 *
 * In numerical relativity, ghost (buffer) zones surround the physical grid
 * to permit centered finite difference stencils at boundaries. The number
 * of ghost zones equals half the stencil width for centered schemes.
 *
 * For 8th-order centered FD: 4 ghost zones needed
 * For 4th-order centered FD: 2 ghost zones needed
 * For 2nd-order centered FD: 1 ghost zone needed
 * --------------------------------------------------------------------------- */
#define NR_MAX_GHOST 6       /**< Maximum ghost zones for high-order stencils */
#define NR_DEFAULT_GHOST 3   /**< Default ghost zones (6th-order stencil) */
#define NR_MAX_DIM 512       /**< Maximum grid points per dimension */
#define NR_STENCIL_ORDER 4   /**< Default finite difference order */

/** Spatial dimensionality of the computational domain */
typedef enum {
    NR_DIM_1D = 1,    /**< 1D problems (spherical symmetry) */
    NR_DIM_2D = 2,    /**< 2D problems (axisymmetry) */
    NR_DIM_3D = 3     /**< 3D problems (full GR) */
} nr_dimension_t;

/** Boundary condition type at each grid face */
typedef enum {
    NR_BC_FLAT      = 0,  /**< Asymptotic flatness: f → 0 at boundary */
    NR_BC_RADIATIVE = 1,  /**< Sommerfeld outgoing radiation condition */
    NR_BC_REFLECTION = 2, /**< Reflection symmetry (e.g., z→-z in octant) */
    NR_BC_PERIODIC  = 3,  /**< Periodic boundary (test problems) */
    NR_BC_STATIC     = 4, /**< Static boundary: ∂_t f = 0 at boundary */
    NR_BC_EXCISION   = 5  /**< Excision: no boundary condition applied inside */
} nr_bc_type_t;

/** Parity of a grid function under reflection */
typedef enum {
    NR_PARITY_EVEN = 0, /**< f(-x) = f(x) — scalars, diagonal metric components */
    NR_PARITY_ODD  = 1  /**< f(-x) = -f(x) — shift vector components */
} nr_parity_t;

/**
 * L1: 3D Cartesian Grid
 *
 * Represents the spatial hypersurface discretized as N_x × N_y × N_z points
 * with uniform spacing dx. Ghost zones surround the physical domain.
 *
 * Total grid extent: [-ghost*dx, (N+ghost)*dx] approx.
 * Physical domain: [0, (N-1)*dx] in index space [ghost, N+ghost-1]
 */
typedef struct {
    int    nx, ny, nz;       /**< Number of interior grid points per dimension */
    int    ng;               /**< Number of ghost zones on each side */
    double xmin, xmax;       /**< Physical domain bounds in x */
    double ymin, ymax;       /**< Physical domain bounds in y */
    double zmin, zmax;       /**< Physical domain bounds in z */
    double dx, dy, dz;       /**< Grid spacing in each direction */
    double *x_coords;        /**< Pre-computed x-coordinates (size: nx+2*ng) */
    double *y_coords;        /**< Pre-computed y-coordinates (size: ny+2*ng) */
    double *z_coords;        /**< Pre-computed z-coordinates (size: nz+2*ng) */
} nr_grid_t;

/**
 * L3: Grid Function — A scalar field on the 3D grid
 *
 * Grid functions are the fundamental dynamical variables in NR:
 * lapse α, conformal factor φ, trace of extrinsic curvature K, etc.
 * Memory layout is row-major: f[i][j][k] = data[i * stride_y + j * stride_z + k]
 */
typedef struct {
    double *data;            /**< Flat array of size (nx+2*ng)*(ny+2*ng)*(nz+2*ng) */
    int     nx, ny, nz;      /**< Grid dimensions including ghosts */
    int     ng;              /**< Number of ghost zones */
    int     stride_y;        /**< = nz + 2*ng */
    int     stride_x;        /**< = (ny + 2*ng) * (nz + 2*ng) */
} nr_gf_t;

/**
 * L3: Symmetric 3-tensor on the grid
 *
 * For the 3-metric γ_{ij} and extrinsic curvature K_{ij}, which are
 * symmetric rank-2 tensors with 6 independent components.
 * Stored as: components[0..5] = {xx, xy, xz, yy, yz, zz}
 */
typedef struct {
    nr_gf_t *comp[6];  /**< Six independent symmetric tensor components */
} nr_sym_tensor3_t;

/**
 * L3: Vector on the grid (3 independent components)
 *
 * For the shift vector β^i and conformal connection Γ̃^i.
 */
typedef struct {
    nr_gf_t *comp[3];  /**< Three vector components {x, y, z} */
} nr_vector3_t;

/* ===========================================================================
 * L1: Grid Function Memory Management
 * =========================================================================== */

/**
 * Allocate a 3D grid function with ghost zones.
 * Complexity: O(N^3) memory. All data initialized to zero.
 *
 * @param nx, ny, nz  Number of interior points per dimension
 * @param ng          Number of ghost zones per side
 * @return            Initialized grid function, or NULL on allocation failure
 */
nr_gf_t* nr_gf_alloc(int nx, int ny, int nz, int ng);

/** Free a grid function and its data array */
void nr_gf_free(nr_gf_t *gf);

/**
 * Access grid function at 3D index (i, j, k), where indices range
 * from 0 to (dim+2*ng-1), inclusive.
 * Complexity: O(1) — direct array access.
 */
static inline double* nr_gf_at(nr_gf_t *gf, int i, int j, int k) {
    return &gf->data[i * gf->stride_x + j * gf->stride_y + k];
}

/**
 * Set all values of a grid function to a constant.
 * Complexity: O(N^3).
 */
void nr_gf_set_all(nr_gf_t *gf, double value);

/**
 * Copy grid function src → dst (must have identical dimensions).
 * Complexity: O(N^3).
 */
void nr_gf_copy(nr_gf_t *dst, const nr_gf_t *src);

/**
 * L2-norm of a grid function over the interior (excluding ghost zones).
 * ||f||_2 = sqrt( Σ_i f_i^2 / N_total ).
 * Complexity: O(N^3).
 */
double nr_gf_l2norm(const nr_gf_t *gf);

/**
 * Maximum absolute value over the interior.
 * ||f||_∞ = max_i |f_i|.
 * Complexity: O(N^3).
 */
double nr_gf_maxabs(const nr_gf_t *gf);

/* ===========================================================================
 * L1: Grid Structure Management
 * =========================================================================== */

/**
 * Initialize a grid structure with uniform spacing.
 * Computes coordinate arrays for all three dimensions.
 * Complexity: O(N) in total grid points.
 */
void nr_grid_init(nr_grid_t *grid, int nx, int ny, int nz,
                  double xmin, double xmax,
                  double ymin, double ymax,
                  double zmin, double zmax,
                  int ng);

/** Free grid coordinate arrays */
void nr_grid_free(nr_grid_t *grid);

/**
 * Get physical x-coordinate at index i (including ghost zone offset).
 * Complexity: O(1).
 */
double nr_grid_x(const nr_grid_t *grid, int i);

/**
 * Get physical y-coordinate at index j.
 * Complexity: O(1).
 */
double nr_grid_y(const nr_grid_t *grid, int j);

/**
 * Get physical z-coordinate at index k.
 * Complexity: O(1).
 */
double nr_grid_z(const nr_grid_t *grid, int k);

/**
 * Compute radial coordinate r = sqrt(x^2 + y^2 + z^2) at grid point (i,j,k).
 * Complexity: O(1).
 */
double nr_grid_r(const nr_grid_t *grid, int i, int j, int k);

/**
 * Compute sin(theta), cos(theta) at grid point for spherical projections.
 * theta = arccos(z/r), returns sin_theta, cos_theta via pointers.
 * Complexity: O(1).
 */
void nr_grid_spherical(const nr_grid_t *grid, int i, int j, int k,
                       double *r, double *sin_theta, double *cos_theta,
                       double *sin_phi, double *cos_phi);

/* ===========================================================================
 * L3: Vector and Tensor Allocation
 * =========================================================================== */

/**
 * Allocate all 6 components of a symmetric 3-tensor.
 * Each component is a grid function with identical ghost structure.
 * Complexity: O(N^3) memory.
 */
nr_sym_tensor3_t* nr_sym_tensor3_alloc(int nx, int ny, int nz, int ng);

/** Free a symmetric 3-tensor */
void nr_sym_tensor3_free(nr_sym_tensor3_t *t);

/**
 * Allocate all 3 components of a spatial vector.
 * Complexity: O(N^3) memory.
 */
nr_vector3_t* nr_vector3_alloc(int nx, int ny, int nz, int ng);

/** Free a spatial vector */
void nr_vector3_free(nr_vector3_t *v);

#ifdef __cplusplus
}
#endif

#endif /* NR_GRID_H */
