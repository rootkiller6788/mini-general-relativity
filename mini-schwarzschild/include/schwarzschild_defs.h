/**
 * schwarzschild_defs.h — Core definitions for Schwarzschild spacetime
 *
 * Reference: Wald (1984) Ch.6, Carroll (2004) Ch.5, MTW (1973) Ch.23
 * MIT 8.962 — General Relativity
 *
 * Covers:
 *   L1 — Definitions: metric tensor, coordinates, constants
 *   L2 — Core Concepts: event horizon, coordinate singularity, asymptotically flat
 */

#ifndef SCHWARZSCHILD_DEFS_H
#define SCHWARZSCHILD_DEFS_H

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * Physical constants (SI units)
 * ========================================================================== */

#define SCHW_G_N  (6.67430e-11)
#define SCHW_C    (299792458.0)
#define SCHW_M_SUN (1.98847e30)
#define SCHW_AU   (1.495978707e11)
#define SCHW_PC   (3.085677581e16)
#define SCHW_M_PL (2.176434e-8)
#define SCHW_L_PL (1.616255e-35)
#define SCHW_GEO_MASS_FACTOR (SCHW_G_N / (SCHW_C * SCHW_C))
#define SCHW_GEO_ENERGY_FACTOR (SCHW_G_N / (SCHW_C * SCHW_C * SCHW_C * SCHW_C))

/* ==========================================================================
 * Coordinate system enum
 * ========================================================================== */

typedef enum {
    SCHW_COORD_STANDARD,
    SCHW_COORD_ISOTROPIC,
    SCHW_COORD_EDDINGTON_FINKELSTEIN_IN,
    SCHW_COORD_EDDINGTON_FINKELSTEIN_OUT,
    SCHW_COORD_KRUSKAL_SZEKERES,
    SCHW_COORD_PAINLEVE_GULLSTRAND,
    SCHW_COORD_LEMAITRE,
    SCHW_COORD_TORTOISE
} SchwarzschildCoordType;

/* ==========================================================================
 * 4D coordinate representation
 * ========================================================================== */

typedef struct {
    double t;
    double r;
    double theta;
    double phi;
} SchwarzschildPoint;

typedef struct {
    double v0;
    double v1;
    double v2;
    double v3;
} SchwarzschildVector;

typedef struct {
    double g00, g01, g02, g03;
    double g11, g12, g13;
    double g22, g23;
    double g33;
} SchwarzschildSymmetric4x4;

/* ==========================================================================
 * Schwarzschild spacetime parameters
 * ========================================================================== */

typedef struct {
    double mass;
    double rs;
    double geometric_mass;
    double horizon_area;
    double surface_gravity;
    double hawking_temp;
    double bekenstein_entropy;
} SchwarzschildBlackHole;

/* ==========================================================================
 * Metric component computation
 * ========================================================================== */

void schwarzschild_metric_covariant(const SchwarzschildPoint *p, double rs,
                                     SchwarzschildSymmetric4x4 *metric);
void schwarzschild_metric_contravariant(const SchwarzschildPoint *p, double rs,
                                         SchwarzschildSymmetric4x4 *inverse_metric);
double schwarzschild_g_tt(double r, double rs);
double schwarzschild_g_rr(double r, double rs);
double schwarzschild_g_theta_theta(double r);
double schwarzschild_g_phi_phi(double r, double theta);
double schwarzschild_metric_determinant_root(double r, double theta);

/* ==========================================================================
 * Schwarzschild radius and horizon properties
 * ========================================================================== */

double schwarzschild_radius(double mass);
double geometric_mass(double mass);
double schwarzschild_surface_gravity(double mass);
double schwarzschild_horizon_area(double mass);

/* ==========================================================================
 * Coordinate transformations
 * ========================================================================== */

double schwarzschild_to_tortoise(double r, double rs);
double tortoise_to_schwarzschild(double r_star, double rs);
double schwarzschild_to_isotropic(double r_schw, double rs);
double isotropic_to_schwarzschild(double r_iso, double rs);
void schwarzschild_to_kruskal(double t, double r, double rs, double *T, double *X);
double schwarzschild_to_eddington_finkelstein_v(double t, double r, double rs);
double schwarzschild_to_eddington_finkelstein_u(double t, double r, double rs);
double schwarzschild_to_painleve_gullstrand_time(double t, double r, double rs);
double schwarzschild_proper_distance(double r1, double r2, double rs);
double schwarzschild_volume_element(const SchwarzschildPoint *p, double rs);

/* ==========================================================================
 * Tetrad (vierbein) for local inertial frame
 * ========================================================================== */

void schwarzschild_tetrad(const SchwarzschildPoint *p, double rs, double tetrad[4][4]);
void schwarzschild_project_to_tetrad(const SchwarzschildVector *v_coord,
                                      const double tetrad[4][4],
                                      SchwarzschildVector *v_tetrad);
void schwarzschild_black_hole_init(double mass, SchwarzschildBlackHole *bh);
int schwarzschild_is_inside_horizon(double r, double rs);
double schwarzschild_lapse(double r, double rs);
void schwarzschild_shift(const SchwarzschildPoint *p, double rs, double beta[3]);

#endif /* SCHWARZSCHILD_DEFS_H */
