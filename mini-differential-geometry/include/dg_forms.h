/**
 * @file dg_forms.h
 * @brief Differential forms, exterior calculus, Hodge duality
 *
 * Knowledge: L1 p-forms/wedge/d, L2 Hodge star/codifferential,
 *            L3 de Rham cohomology, L4 Stokes theorem/Cartan eqs
 * Reference: Wald Appendix B, Nakahara Sec.5-6
 * Courses: MIT 8.962, Oxford CMT, Cambridge Part III
 */

#ifndef DG_FORMS_H
#define DG_FORMS_H

#include "dg_metric.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int     degree;
    int     dim;
    double *components;
    int     chart_index;
} DifferentialForm;

DifferentialForm* form_alloc(int degree, int dim, int chart_index);
DifferentialForm* form_create(int degree, int dim,
                               const double *components, int chart_index);
DifferentialForm* form_clone(const DifferentialForm *omega);
void form_free(DifferentialForm *omega);

double form_get(const DifferentialForm *omega, const int *indices);
void   form_set(DifferentialForm *omega, const int *indices, double value);

DifferentialForm* form_wedge(const DifferentialForm *alpha,
                              const DifferentialForm *beta);
DifferentialForm* form_exterior_derivative(const DifferentialForm *omega,
                                            const double *domega);
DifferentialForm* form_exterior_derivative_fd(const DifferentialForm *omega,
                                               double h);
DifferentialForm* form_interior_product(const double *X,
                                         const DifferentialForm *omega);
DifferentialForm* form_hodge_star(const DifferentialForm *omega,
                                   const Metric *metric);
DifferentialForm* form_laplace_de_rham(const DifferentialForm *omega,
                                        const double *domega,
                                        const Metric *metric);

double form_verify_stokes(const DifferentialForm *omega, const double *domega,
                           const Metric *metric);
int  form_is_closed(const DifferentialForm *omega, const double *domega,
                     double tolerance);
int  form_is_exact(const DifferentialForm *omega, double tolerance);
void form_print(const DifferentialForm *omega);

#ifdef __cplusplus
}
#endif
#endif