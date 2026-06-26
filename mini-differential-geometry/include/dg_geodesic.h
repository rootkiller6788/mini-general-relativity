/**
 * @file dg_geodesic.h
 * @brief Geodesics, parallel transport, and geodesic deviation
 *
 * Knowledge: L1 Geodesic eq, L2 Parallel transport/deviation,
 *            L3 Affine parameter, L4 From action principle,
 *            L5 Numerical RK4 integration, L6 Schwarzschild/FRW geodesics
 * Reference: Wald Sec.3.3, Carroll Sec.3, Schutz Sec.6-7
 * Courses: MIT 8.962 Sec.3, Caltech Ph 205, Oxford CMT
 */

#ifndef DG_GEODESIC_H
#define DG_GEODESIC_H

#include "dg_curvature.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DG_MAX_GEO_STEPS 100000

typedef struct {
    double pos[DG_MAX_DIM];
    double vel[DG_MAX_DIM];
    double tau;
    int    dim;
    int    chart_index;
} GeodesicPoint;

typedef struct {
    double step_size;
    int    max_steps;
    double tolerance;
    int    record_every;
} GeodesicParams;

void geodesic_rhs(double tau, const double *y, double *dydt, void *params);
GeodesicPoint* geodesic_integrate(const GeodesicPoint *initial,
                                   const Connection *conn,
                                   const GeodesicParams *params,
                                   int *n_steps);
int geodesic_step_rk4(GeodesicPoint *point, const Connection *conn, double dtau);

double geodesic_velocity_norm(const GeodesicPoint *point, const Metric *metric);
double proper_time_total(const GeodesicPoint *points, int n, const Metric *metric);

int parallel_transport(const double *V_initial, const GeodesicPoint *points,
                        int n, const Connection *conn, double *V_result);

void geodesic_deviation_rhs(const double *xi, const double *dxi_dtau,
                             const GeodesicPoint *point, const RiemannTensor *R,
                             double *ddxi);

void geodesic_schwarzschild_init(double r0, double b, double M, GeodesicPoint *point);
void geodesic_frw_null_init(double a0, GeodesicPoint *point);

void geodesic_point_print(const GeodesicPoint *point);
void geodesic_points_free(GeodesicPoint *points);

#ifdef __cplusplus
}
#endif
#endif