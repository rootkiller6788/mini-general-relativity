/**
 * @file dg_tensor.h
 * @brief Tensor algebra on manifolds - type (r,s) tensors and operations
 *
 * Knowledge: L1 Tensor def, L2 Index ops, L3 Product/contraction/symmetrization
 * Reference: Wald Sec.2.3-2.4, Carroll Sec.2, Schutz Sec.3
 * Courses: MIT 8.962, Stanford PHYSICS 230, Cambridge Part III
 */

#ifndef DG_TENSOR_H
#define DG_TENSOR_H

#include "dg_manifold.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DG_MAX_RANK 4

typedef struct {
    int     rank_r;
    int     rank_s;
    int     dim;
    int     chart_index;
    double *components;
    int     owns_memory;
} Tensor;

Tensor* tensor_alloc(int r, int s, int dim);
Tensor* tensor_create(int r, int s, int dim, const double *components, int chart_index);
Tensor* tensor_clone(const Tensor *src);
void    tensor_free(Tensor *t);

int    tensor_flat_index(const Tensor *t, const int *idx);
double tensor_get(const Tensor *t, const int *idx);
void   tensor_set(Tensor *t, const int *idx, double value);

Tensor* tensor_add(const Tensor *a, const Tensor *b);
Tensor* tensor_subtract(const Tensor *a, const Tensor *b);
Tensor* tensor_scale(const Tensor *t, double scalar);

Tensor* tensor_product(const Tensor *a, const Tensor *b);
Tensor* tensor_contract(const Tensor *t, int contra_pos, int cov_pos);
Tensor* tensor_raise_index(const Tensor *t, const Tensor *g_inv, int cov_pos);
Tensor* tensor_lower_index(const Tensor *t, const Tensor *g, int contra_pos);

Tensor* tensor_transform(const Tensor *t, const double *jacobian,
                          const double *jacobian_inv, int new_chart);

Tensor* tensor_symmetrize_02(const Tensor *t);
Tensor* tensor_antisymmetrize_02(const Tensor *t);
int     tensor_is_symmetric_02(const Tensor *t);
int     tensor_is_antisymmetric_02(const Tensor *t);

double tensor_norm(const Tensor *t);
double tensor_trace_11(const Tensor *t);
void   tensor_print(const Tensor *t);

#ifdef __cplusplus
}
#endif
#endif /* DG_TENSOR_H */