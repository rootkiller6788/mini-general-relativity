/**
 * @file dg_metric.h
 * @brief Metric tensor g_{mu nu} - the fundamental field of GR
 *
 * Knowledge: L1 Metric def/signature, L2 Proper time/length/volume element,
 *            L3 Pseudo-Riemannian structure, L4 Signature verification
 * Reference: Wald Sec.3, Carroll Sec.2-3, Schutz Sec.2
 * Courses: MIT 8.962 Sec.3, Princeton PHY 535, ETH 402-0891
 */

#ifndef DG_METRIC_H
#define DG_METRIC_H

#include "dg_tensor.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SIG_LORENTZIAN_MPPP = 0,
    SIG_LORENTZIAN_PMMM = 1,
    SIG_RIEMANNIAN      = 2
} MetricSignature;

typedef struct {
    int             dim;
    MetricSignature signature;
    double         *g;
    double         *g_inv;
    double          det;
    int             chart_index;
} Metric;

Metric* metric_alloc(int dim, MetricSignature signature);
Metric* metric_create(int dim, MetricSignature signature,
                       const double *components, int chart_index);
Metric* metric_clone(const Metric *m);
void    metric_free(Metric *m);

double metric_get(const Metric *m, int mu, int nu);
void   metric_set(Metric *m, int mu, int nu, double value);
double metric_inv_get(const Metric *m, int mu, int nu);
int    metric_recompute_inv_det(Metric *m);

double metric_line_element(const Metric *m, const double dx[]);
double metric_proper_time_inc(const Metric *m, const double dx[]);
double metric_proper_length_inc(const Metric *m, const double dx[]);
int    metric_classify_dx(const Metric *m, const double dx[]);

double metric_volume_element(const Metric *m);
double metric_determinant(const Metric *m);
int    metric_verify_signature(const Metric *m);

void   metric_print(const Metric *m);

#ifdef __cplusplus
}
#endif
#endif