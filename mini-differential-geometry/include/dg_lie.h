/**
 * @file dg_lie.h
 * @brief Lie derivative and Killing vector fields
 *
 * Knowledge: L1 Lie derivative, L2 Killing vectors/isometries,
 *            L3 Lie bracket/Frobenius, L4 Symmetries and conserved quantities
 * Reference: Wald Appendix C, Carroll Sec.3, Schutz Sec.6
 * Courses: MIT 8.962, Princeton PHY 535, Cambridge Part III
 */

#ifndef DG_LIE_H
#define DG_LIE_H

#include "dg_connection.h"
#include "dg_geodesic.h"

#ifdef __cplusplus
extern "C" {
#endif

double lie_derivative_scalar(const double *X, const double *df, int dim);
void lie_derivative_vector(const double *X, const double *Y,
                            const double *dX, const double *dY,
                            int dim, double *LXY);
void lie_derivative_covector(const double *X, const double *omega,
                              const double *dX, const double *domega,
                              int dim, double *Lw);
void lie_derivative_tensor02(const double *X, const double *T,
                              const double *dX, const double *dT,
                              int dim, double *L_T);

double killing_equation_check(const double *X, const double *dX,
                               const Metric *g, const double *dg);
double* killing_schwarzschild(int dim, int *n_killing);
double killing_conserved_quantity(const double *xi,
                                   const GeodesicPoint *point,
                                   const Metric *metric);
int max_killing_vectors(int dim);

#ifdef __cplusplus
}
#endif
#endif