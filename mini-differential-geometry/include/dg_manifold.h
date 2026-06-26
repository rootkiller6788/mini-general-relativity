/**
 * @file dg_manifold.h
 * @brief Manifold, Chart, and Atlas definitions for differential geometry
 *
 * Knowledge Coverage: L1 Definitions (Manifold, Chart, Atlas, Tangent Space)
 * Reference: Wald §2, Carroll §2, do Carmo §0-§1
 * Courses: MIT 8.962, Stanford PHYSICS 230, Cambridge Part III GR
 *
 * A smooth (C^\infty) manifold M of dimension n is a topological space
 * locally homeomorphic to R^n, equipped with a maximal smooth atlas.
 */

#ifndef DG_MANIFOLD_H
#define DG_MANIFOLD_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Maximums for static allocation in GR context */
#define DG_MAX_DIM      4
#define DG_MAX_CHARTS  16
#define DG_MAX_NAME    64
#define DG_MAX_COORD_NAME 16

/**
 * @brief A chart (coordinate system) on a manifold.
 *
 * (U, phi) where U is an open set in M and phi: U -> R^n is a homeomorphism.
 */
typedef struct {
    char   coord_names[DG_MAX_DIM][DG_MAX_COORD_NAME];
    double domain_min[DG_MAX_DIM];
    double domain_max[DG_MAX_DIM];
    int    dim;
    bool   is_periodic[DG_MAX_DIM];
    double period[DG_MAX_DIM];
} Chart;

/**
 * @brief A smooth manifold defined by an atlas of compatible charts.
 */
typedef struct {
    int     dim;
    char    name[DG_MAX_NAME];
    int     num_charts;
    Chart   charts[DG_MAX_CHARTS];
    bool    is_compact;
    bool    is_oriented;
    bool    has_boundary;
} Manifold;

/**
 * @brief Tangent vector X = X^mu d_mu at a point p in M.
 *
 * T_pM is an n-dimensional real vector space. In coordinates,
 * tangent vectors transform as X'^mu = (dx'^mu / dx^nu) X^nu.
 */
typedef struct {
    double components[DG_MAX_DIM];
    int    dim;
    int    chart_index;
} TangentVector;

/**
 * @brief Cotangent vector (1-form) omega = omega_mu dx^mu at a point.
 *
 * T_p^*M is the dual space. 1-forms transform as
 * omega'_mu = (dx^nu / dx'^mu) omega_nu.
 */
typedef struct {
    double components[DG_MAX_DIM];
    int    dim;
    int    chart_index;
} CotangentVector;

/* ---- API ---- */

Manifold manifold_create(int dim, const char *name);

int manifold_add_chart(Manifold *m, const Chart *c);

Chart chart_create(int dim,
                   const char coord_names[][DG_MAX_COORD_NAME],
                   const double domain_min[],
                   const double domain_max[]);

bool chart_contains_point(const Chart *c, const double x[]);

/**
 * @brief Transition map Jacobian between two overlapping charts.
 *
 * Given charts (U_alpha, phi_alpha) and (U_beta, phi_beta) and a point
 * p in U_alpha ∩ U_beta, compute the Jacobian
 * J^mu_nu = d(x_beta^mu) / d(x_alpha^nu)
 * of the transition map phi_beta ∘ phi_alpha^{-1}.
 *
 * @param m          Manifold
 * @param chart_src  Source chart index
 * @param chart_dst  Destination chart index
 * @param x_src      Point coordinates in source chart
 * @param jacobian   Output: dim x dim Jacobian matrix (row-major)
 * @return           0 on success, -1 on failure
 *
 * Uses registered transition function callback.
 * Course Ref: Wald §2.1 — Transition Functions
 */
int manifold_transition_jacobian(const Manifold *m,
                                  int chart_src, int chart_dst,
                                  const double x_src[],
                                  double jacobian[]);

/** Callback type for transition map x_dst = f(x_src) */
typedef void (*transition_func_t)(const double *x_src, double *x_dst, int dim);

int manifold_set_transition_map(Manifold *m,
                                 int chart_src, int chart_dst,
                                 transition_func_t func);

bool manifold_atlas_is_complete(const Manifold *m);

TangentVector tangent_vector_create(int dim, const double components[],
                                     int chart_index);

CotangentVector cotangent_vector_create(int dim, const double components[],
                                         int chart_index);

/**
 * @brief Action of a 1-form on a tangent vector: omega(X) = omega_mu X^mu.
 * @param omega Cotangent vector
 * @param X     Tangent vector
 * @return      The scalar omega_mu X^mu
 *
 * Both must reference the same chart and dimension.
 * Course Ref: Wald §2.3 — Dual Vector Spaces
 */
double cotangent_act_on_tangent(const CotangentVector *omega,
                                 const TangentVector *X);

void manifold_destroy(Manifold *m);

#ifdef __cplusplus
}
#endif

#endif /* DG_MANIFOLD_H */
