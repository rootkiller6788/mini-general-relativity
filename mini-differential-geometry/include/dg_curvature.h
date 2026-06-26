/**
 * @file dg_curvature.h
 * @brief Riemann, Ricci, Einstein, and Weyl curvature tensors
 *
 * Knowledge: L1 Riemann/Ricci/Einstein, L2 Geodesic deviation,
 *            L3 Sectional curvature/Weyl, L4 Bianchi identities
 * Reference: Wald Sec.3.2-3.4, Carroll Sec.3, Schutz Sec.6
 * Courses: MIT 8.962 Sec.3, Cambridge Part III, ETH 402-0891
 */

#ifndef DG_CURVATURE_H
#define DG_CURVATURE_H

#include "dg_connection.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int     dim;
    double *R;
    int     chart_index;
} RiemannTensor;

typedef struct {
    int     dim;
    double *R;
    int     chart_index;
} RicciTensor;

RiemannTensor* riemann_alloc(int dim, int chart_index);
RiemannTensor* riemann_compute(const Connection *conn, const double *dGamma);
double riemann_get(const RiemannTensor *R, int rho, int sigma, int mu, int nu);
double* riemann_lower_first_index(const RiemannTensor *R, const Metric *g);

RicciTensor* ricci_alloc(int dim, int chart_index);
RicciTensor* ricci_compute(const RiemannTensor *R);
double ricci_get(const RicciTensor *Ric, int mu, int nu);

double ricci_scalar_compute(const RicciTensor *Ric, const Metric *g);
double* einstein_tensor_compute(const RicciTensor *Ric, const Metric *g, double R);
double einstein_verify_divergence_free(const double *G, const Connection *conn,
                                        const double *dG, const Metric *g);

double* weyl_tensor_compute(const RiemannTensor *R, const RicciTensor *Ric,
                             const Metric *g, double R_sc);
double sectional_curvature(const double *R_cov, const Metric *g,
                            const double *X, const double *Y, int dim);
double bianchi_identity_first_verify(const RiemannTensor *R);
double bianchi_identity_second_verify(const RiemannTensor *R,
                                       const Connection *conn, const double *dR);
double kretschmann_scalar(const double *R_cov, const Metric *g, int dim);

void riemann_free(RiemannTensor *R);
void ricci_free(RicciTensor *Ric);
void riemann_print(const RiemannTensor *R);
void ricci_print(const RicciTensor *Ric);

#ifdef __cplusplus
}
#endif
#endif