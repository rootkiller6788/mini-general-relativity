/**
 * @file dg_connection.h
 * @brief Affine connection and covariant derivative
 *
 * Knowledge: L1 Christoffel symbols, L2 Levi-Civita connection,
 *            L3 Connection 1-forms, L4 Covariant derivative of tensors
 * Reference: Wald Sec.3.1-3.2, Carroll Sec.3, Schutz Sec.6
 * Courses: MIT 8.962 Sec.3, Oxford CMT, Caltech Ph 205
 */

#ifndef DG_CONNECTION_H
#define DG_CONNECTION_H

#include "dg_metric.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int     dim;
    double *Gamma;
    int     chart_index;
} Connection;

Connection* connection_alloc(int dim, int chart_index);
Connection* connection_from_metric(const Metric *metric, const double *dmetric);
Connection* connection_from_metric_fd(const Metric *metric, double h);
void        connection_free(Connection *conn);

double connection_get(const Connection *conn, int rho, int mu, int nu);
void   connection_set(Connection *conn, int rho, int mu, int nu, double value);

double covariant_derivative_scalar(const double *dphi, int dim, int mu);
double covariant_derivative_vector(const Connection *conn,
                                    const double *V, const double *dV,
                                    int mu, int nu);
double covariant_derivative_covector(const Connection *conn,
                                      const double *omega, const double *domega,
                                      int mu, int nu);
double covariant_derivative_tensor02(const Connection *conn,
                                      const double *T, const double *dT,
                                      int rho, int mu, int nu);
void covariant_derivative_vector_all(const Connection *conn,
                                      const double *V, const double *dV,
                                      double *nabla_V);

double connection_verify_metric_compatibility(const Connection *conn,
                                               const Metric *metric,
                                               const double *dg);
double connection_compute_torsion(const Connection *conn);
void connection_print(const Connection *conn);

#ifdef __cplusplus
}
#endif
#endif