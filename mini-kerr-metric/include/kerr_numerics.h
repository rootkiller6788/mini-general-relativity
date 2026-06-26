/**
 * @file kerr_numerics.h
 * @brief Numerical methods for Kerr geometry: root finding, integration,
 *        spectral methods, and ray tracing
 *
 * Computing observables in Kerr geometry requires robust numerical methods:
 * - Root finding for horizons, turning points, impact parameters
 * - ODE integration for geodesics and radiation transport
 * - Spectral methods for Teukolsky equation
 * - Ray tracing for black hole imaging
 *
 * MIT 8.962, Caltech Ph 205, ETH 402-0891
 */

#ifndef KERR_NUMERICS_H
#define KERR_NUMERICS_H

#include "kerr_metric.h"
#include "kerr_horizons.h"
#include "kerr_geodesics.h"

/* ==========================================================================
 * L1 -- Numerical Integration Structures
 * ========================================================================== */

/** ODE integrator configuration */
typedef struct {
    double rel_tol;       /* Relative tolerance */
    double abs_tol;       /* Absolute tolerance */
    double initial_step;  /* Initial step size */
    double max_step;      /* Maximum step size */
    double min_step;      /* Minimum step size */
    int max_iterations;   /* Maximum integration steps */
    int stiff_check;      /* 1 to enable stiffness detection */
} KerrODESettings;

/** Ray tracing configuration for BH imaging */
typedef struct {
    double camera_distance;  /* Distance from BH to camera */
    double iota;             /* Inclination angle */
    int pixels_x;            /* Image resolution */
    int pixels_y;
    double fov;              /* Field of view in GM/c^2 units */
    int max_steps;           /* Max steps per ray */
    double step_size;        /* Integration step size */
} KerrRayTraceConfig;

/** Ray tracing result for a single ray */
typedef struct {
    double alpha;          /* Impact parameter alpha (image plane) */
    double beta;           /* Impact parameter beta */
    int captured;          /* 1 if captured by BH */
    double r_min;          /* Minimum r reached */
    double deflection;     /* Total deflection angle */
    double redshift;       /* Gravitational + Doppler redshift */
    double path_length;    /* Optical path length */
} KerrRayResult;

/** Numerical derivative configuration */
typedef struct {
    double h;              /* Step size for finite difference */
    int method;            /* 0=forward, 1=central, 2=Richardson */
    int order;             /* Order of accuracy */
} KerrDiffSettings;

/* ==========================================================================
 * L2 -- Root Finding
 * ========================================================================== */

/**
 * Find the radial coordinate r where f(r)=0 using Brent's method.
 * Used for: horizon radii, turning points, photon orbit impact params.
 *
 * @param f    Function to find root of
 * @param ctx  Context pointer passed to f
 * @param a    Left bracket
 * @param b    Right bracket
 * @param tol  Tolerance
 * @return Root value, or NAN if not found
 */
double kerr_brent_root(double (*f)(double, void*), void *ctx,
                       double a, double b, double tol);

/**
 * Find all real positive roots of a cubic equation a3 x^3 + a2 x^2 + a1 x + a0 = 0.
 * Used for solving radial turning point equations.
 *
 * Returns number of real positive roots. Roots stored in roots[0..2].
 */
int kerr_cubic_real_roots(double a3, double a2, double a1, double a0,
                          double roots[3]);

/**
 * Find all real positive roots of a quartic a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0.
 * Used for KS->BL coordinate conversion.
 *
 * Returns number of real positive roots.
 */
int kerr_quartic_real_roots(double a4, double a3, double a2,
                            double a1, double a0, double roots[4]);

/**
 * Newton-Raphson root finding with analytical derivative.
 */
double kerr_newton_raphson(double (*f)(double, void*),
                           double (*df)(double, void*),
                           void *ctx, double x0, double tol,
                           int max_iter);

/* ==========================================================================
 * L3 -- ODE Integration
 * ========================================================================== */

/**
 * Initialize ODE settings with sensible defaults for Kerr geodesics.
 */
KerrODESettings kerr_ode_settings_default(void);

/**
 * Adaptive Runge-Kutta-Fehlberg 4(5) step for systems of ODEs.
 *
 * dy/dt = f(t, y), where y is an n-dimensional state vector.
 * Returns the suggested next step size.
 */
double kerr_rkf45_step(void (*f)(double, const double*, double*, void*),
                       void *ctx, double t, double *y, int n,
                       double *dt, const KerrODESettings *settings);

/**
 * Integrate ODE system with adaptive step size control from t0 to t_end.
 *
 * Returns number of steps taken.
 */
int kerr_ode_integrate(void (*f)(double, const double*, double*, void*),
                       void *ctx, double *y, int n,
                       double t0, double t_end,
                       const KerrODESettings *settings,
                       double (*callback)(double, const double*, void*),
                       void *cb_ctx);

/* ==========================================================================
 * L4 -- Spectral Methods
 * ========================================================================== */

/**
 * Compute Chebyshev polynomial T_n(x) of degree n at point x in [-1, 1].
 *
 * Using recurrence: T_0=1, T_1=x, T_n=2x T_{n-1} - T_{n-2}
 */
double kerr_chebyshev(int n, double x);

/**
 * Chebyshev differentiation matrix D_N for N+1 Chebyshev points.
 * Given function values at Chebyshev points, compute derivative values.
 *
 * @param N    Number of Chebyshev points (order N approximation)
 * @param D    [out] (N+1)x(N+1) differentiation matrix (row-major)
 */
int kerr_chebyshev_diff_matrix(int N, double *D);

/**
 * Chebyshev collocation points: x_j = cos(pi j / N), j = 0..N.
 */
int kerr_chebyshev_points(int N, double *x);

/**
 * Barycentric interpolation using Chebyshev points.
 *
 * Given function values at Chebyshev points, evaluate at arbitrary x.
 */
double kerr_chebyshev_interpolate(const double *x, const double *f,
                                  int N, double x_eval);

/**
 * Solve a radial Teukolsky equation eigenproblem using
 * a Chebyshev spectral method (continued fraction method).
 *
 * The Teukolsky equation: d^2 R / dr*^2 + V(r*) R = 0
 * Describes perturbations of Kerr geometry.
 *
 * @param s     Spin weight (-2 for gravitational, -1 for EM, 0 for scalar)
 * @param l     Angular harmonic index
 * @param m     Azimuthal harmonic index
 * @param omega Complex frequency
 * @return Radial Teukolsky function value at r=r_plus
 */
int kerr_teukolsky_radial_spectral(const KerrBlackHole *bh,
                                    int s, int l, int m,
                                    double omega_re, double omega_im,
                                    double *R_horizon_real,
                                    double *R_horizon_imag);

/* ==========================================================================
 * L5 -- Ray Tracing
 * ========================================================================== */

/**
 * Initialize a ray tracing configuration for BH imaging.
 */
KerrRayTraceConfig kerr_raytrace_config_default(const KerrBlackHole *bh);

/**
 * Trace a single null geodesic from the observer''s image plane backward
 * in time to determine if it hits the BH or escapes to infinity.
 *
 * @param bh     Black hole parameters
 * @param config Ray tracing configuration
 * @param alpha  Horizontal coordinate on image plane
 * @param beta   Vertical coordinate on image plane
 * @param result [out] Ray tracing result
 * @return 0 on success
 */
int kerr_raytrace_single(const KerrBlackHole *bh,
                         const KerrRayTraceConfig *config,
                         double alpha, double beta,
                         KerrRayResult *result);

/**
 * Compute the redshift of a photon from an accretion disk element
 * at location (r, theta, phi) to a distant observer.
 *
 * Combines gravitational redshift and relativistic Doppler shift
 * from Keplerian disk motion.
 */
double kerr_photon_redshift(const KerrBlackHole *bh,
                            const KerrBLPoint *emitter,
                            double iota, int prograde);

/**
 * Compute the light bending angle for a photon passing near a Kerr BH.
 *
 * Uses the full null geodesic equation with the Carter constant.
 */
double kerr_light_bending(const KerrBlackHole *bh,
                          double impact_param, int prograde);

/* ==========================================================================
 * L6 -- Numerical Derivatives and Integration
 * ========================================================================== */

/**
 * Compute numerical partial derivative of metric component using
 * Richardson extrapolation with central differences.
 *
 * Accuracy: O(h^4) with Richardson extrapolation.
 */
double kerr_numerical_metric_derivative(const KerrBlackHole *bh,
                                        const KerrBLPoint *pt,
                                        int mu, int nu, int coord_dir,
                                        const KerrDiffSettings *settings);

/**
 * Compute the Levi-Civita connection numerically and compare with
 * analytical Christoffel symbols. Returns max absolute error.
 */
double kerr_christoffel_verify(const KerrBlackHole *bh,
                               const KerrBLPoint *pt, double h);

/**
 * 4th-order Runge-Kutta integration step (fixed step size).
 * Generic implementation for Kerr geodesics.
 */
int kerr_rk4_fixed_step(void (*f)(double, const double*, double*, void*),
                        void *ctx, double t, double *y, int n, double h);

/**
 * Adams-Bashforth-Moulton predictor-corrector method (order 4).
 * Requires past 4 steps.
 */
int kerr_adams_bashforth_moulton(void (*f)(double, const double*, double*, void*),
                                  void *ctx, double t, double *y, int n,
                                  double h, double **history);

/* ==========================================================================
 * L7 -- Monte Carlo and Optimization
 * ========================================================================== */

/**
 * Compute geodesic deviation via numerical integration of the
 * Jacobi equation along a reference geodesic.
 *
 * d^2 xi^mu / dlambda^2 + 2 Gamma^mu_nurho (dx^nu/dlambda) (d xi^rho/dlambda)
 * + (partial_nu Gamma^mu_rhosigma) (dx^rho/dlambda)(dx^sigma/dlambda) xi^nu = 0
 */
int kerr_geodesic_deviation(const KerrBlackHole *bh,
                            const KerrGeodesicState *ref_geodesic,
                            double initial_dev[4],
                            double dlambda, int n_steps,
                            double deviation[4]);

/**
 * Compute tidal forces on an extended body falling into Kerr BH.
 * The tidal tensor C_ij = R_i0j0 (components of Riemann).
 */
int kerr_tidal_tensor(const KerrBlackHole *bh, const KerrBLPoint *pt,
                      double tidal[3][3]);

/**
 * Find the relativistic image of a point source behind a Kerr BH
 * using the photon shell method (multiple images due to strong lensing).
 *
 * Returns number of images found.
 */
int kerr_strong_lensing_images(const KerrBlackHole *bh,
                               double source_alpha, double source_beta,
                               double iota, int max_images,
                               double *img_alpha, double *img_beta);

#endif /* KERR_NUMERICS_H */