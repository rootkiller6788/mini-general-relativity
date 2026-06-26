#include "dg_tensor.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static int tensor_ncomp(int r, int s, int dim) {
    int n = 1;
    for (int i = 0; i < r + s; i++) n *= dim;
    return n;
}

Tensor* tensor_alloc(int r, int s, int dim) {
    if (r < 0 || s < 0 || dim < 1 || dim > DG_MAX_DIM) return NULL;
    if (r + s > DG_MAX_RANK) return NULL;
    Tensor *t = (Tensor*)malloc(sizeof(Tensor));
    if (!t) return NULL;
    int n = tensor_ncomp(r, s, dim);
    t->components = (double*)calloc(n, sizeof(double));
    if (!t->components) { free(t); return NULL; }
    t->rank_r = r; t->rank_s = s; t->dim = dim;
    t->chart_index = 0; t->owns_memory = 1;
    return t;
}

Tensor* tensor_create(int r, int s, int dim, const double *comp, int chart) {
    Tensor *t = tensor_alloc(r, s, dim);
    if (!t) return NULL;
    int n = tensor_ncomp(r, s, dim);
    if (comp) memcpy(t->components, comp, n * sizeof(double));
    t->chart_index = chart;
    return t;
}

Tensor* tensor_clone(const Tensor *src) {
    if (!src) return NULL;
    return tensor_create(src->rank_r, src->rank_s, src->dim,
                          src->components, src->chart_index);
}

void tensor_free(Tensor *t) {
    if (!t) return;
    if (t->owns_memory && t->components) free(t->components);
    free(t);
}

int tensor_flat_index(const Tensor *t, const int *idx) {
    if (!t || !idx) return -1;
    int flat = 0, rank = t->rank_r + t->rank_s;
    for (int k = 0; k < rank; k++) {
        if (idx[k] < 0 || idx[k] >= t->dim) return -1;
        flat = flat * t->dim + idx[k];
    }
    return flat;
}

double tensor_get(const Tensor *t, const int *idx) {
    int fi = tensor_flat_index(t, idx);
    return (fi >= 0) ? t->components[fi] : 0.0;
}

void tensor_set(Tensor *t, const int *idx, double val) {
    int fi = tensor_flat_index(t, idx);
    if (fi >= 0) t->components[fi] = val;
}

Tensor* tensor_add(const Tensor *a, const Tensor *b) {
    if (!a || !b) return NULL;
    if (a->rank_r != b->rank_r || a->rank_s != b->rank_s || a->dim != b->dim)
        return NULL;
    Tensor *c = tensor_alloc(a->rank_r, a->rank_s, a->dim);
    if (!c) return NULL;
    int n = tensor_ncomp(a->rank_r, a->rank_s, a->dim);
    for (int i = 0; i < n; i++)
        c->components[i] = a->components[i] + b->components[i];
    return c;
}

Tensor* tensor_subtract(const Tensor *a, const Tensor *b) {
    if (!a || !b) return NULL;
    if (a->rank_r != b->rank_r || a->rank_s != b->rank_s || a->dim != b->dim)
        return NULL;
    Tensor *c = tensor_alloc(a->rank_r, a->rank_s, a->dim);
    if (!c) return NULL;
    int n = tensor_ncomp(a->rank_r, a->rank_s, a->dim);
    for (int i = 0; i < n; i++)
        c->components[i] = a->components[i] - b->components[i];
    return c;
}

Tensor* tensor_scale(const Tensor *t, double scalar) {
    if (!t) return NULL;
    Tensor *c = tensor_clone(t);
    if (!c) return NULL;
    int n = tensor_ncomp(t->rank_r, t->rank_s, t->dim);
    for (int i = 0; i < n; i++)
        c->components[i] *= scalar;
    return c;
}

Tensor* tensor_product(const Tensor *a, const Tensor *b) {
    if (!a || !b || a->dim != b->dim) return NULL;
    int rC = a->rank_r + b->rank_r, sC = a->rank_s + b->rank_s;
    if (rC + sC > DG_MAX_RANK) return NULL;
    Tensor *c = tensor_alloc(rC, sC, a->dim);
    if (!c) return NULL;
    int nA = tensor_ncomp(a->rank_r, a->rank_s, a->dim);
    int nB = tensor_ncomp(b->rank_r, b->rank_s, b->dim);
    for (int ia = 0; ia < nA; ia++) {
        double va = a->components[ia];
        if (va == 0.0) continue;
        int base = ia * nB;
        for (int ib = 0; ib < nB; ib++)
            c->components[base + ib] = va * b->components[ib];
    }
    return c;
}

Tensor* tensor_contract(const Tensor *t, int cpos, int vpos) {
    if (!t || t->rank_r < 1 || t->rank_s < 1) return NULL;
    if (cpos < 0 || cpos >= t->rank_r || vpos < 0 || vpos >= t->rank_s) return NULL;
    int rC = t->rank_r - 1, sC = t->rank_s - 1, dim = t->dim;
    Tensor *c = tensor_alloc(rC, sC, dim);
    if (!c) return NULL;
    int nC = tensor_ncomp(rC, sC, dim);
    int rankC = rC + sC;
    for (int fc = 0; fc < nC; fc++) {
        int ires[DG_MAX_RANK], isrc[DG_MAX_RANK];
        int rem = fc;
        for (int k = rankC - 1; k >= 0; k--) { ires[k] = rem % dim; rem /= dim; }
        int si = 0;
        for (int p = 0; p < t->rank_r; p++) {
            if (p == cpos) isrc[si] = -1;
            else isrc[si] = ires[(p < cpos) ? p : p - 1];
            si++;
        }
        for (int p = 0; p < t->rank_s; p++) {
            if (p == vpos) isrc[si] = -1;
            else isrc[si] = ires[rC + ((p < vpos) ? p : p - 1)];
            si++;
        }
        double sum = 0.0;
        for (int k = 0; k < dim; k++) {
            int tmp[DG_MAX_RANK];
            for (int q = 0; q < t->rank_r + t->rank_s; q++)
                tmp[q] = (isrc[q] == -1) ? k : isrc[q];
            int fi = tensor_flat_index(t, tmp);
            if (fi >= 0) sum += t->components[fi];
        }
        c->components[fc] = sum;
    }
    return c;
}

Tensor* tensor_raise_index(const Tensor *t, const Tensor *ginv, int cov_pos) {
    if (!t || !ginv || t->rank_s < 1 || ginv->rank_r != 2 || ginv->rank_s != 0) return NULL;
    if (cov_pos < 0 || cov_pos >= t->rank_s) return NULL;
    int rC = t->rank_r + 1, sC = t->rank_s - 1, dim = t->dim;
    Tensor *c = tensor_alloc(rC, sC, dim);
    if (!c) return NULL;
    int nC = tensor_ncomp(rC, sC, dim);
    for (int fc = 0; fc < nC; fc++) {
        int ires[DG_MAX_RANK], rem = fc;
        for (int k = rC + sC - 1; k >= 0; k--) { ires[k] = rem % dim; rem /= dim; }
        int new_mu = ires[t->rank_r];
        double sum = 0.0;
        for (int sig = 0; sig < dim; sig++) {
            int isrc[DG_MAX_RANK], si = 0;
            for (int p = 0; p < t->rank_r; p++)
                isrc[si++] = ires[(p < t->rank_r) ? p : 0];
            for (int p = 0; p < t->rank_s; p++) {
                if (p == cov_pos) isrc[si++] = sig;
                else {
                    int rp = (p < cov_pos) ? p : p - 1;
                    isrc[si++] = ires[rC + rp];
                }
            }
            int fi = tensor_flat_index(t, isrc);
            double tv = (fi >= 0) ? t->components[fi] : 0.0;
            int gi[2] = {new_mu, sig};
            sum += tensor_get(ginv, gi) * tv;
        }
        c->components[fc] = sum;
    }
    return c;
}

Tensor* tensor_lower_index(const Tensor *t, const Tensor *g, int contra_pos) {
    if (!t || !g || t->rank_r < 1 || g->rank_r != 0 || g->rank_s != 2) return NULL;
    if (contra_pos < 0 || contra_pos >= t->rank_r) return NULL;
    int rC = t->rank_r - 1, sC = t->rank_s + 1, dim = t->dim;
    Tensor *c = tensor_alloc(rC, sC, dim);
    if (!c) return NULL;
    int nC = tensor_ncomp(rC, sC, dim);
    for (int fc = 0; fc < nC; fc++) {
        int ires[DG_MAX_RANK], rem = fc;
        for (int k = rC + sC - 1; k >= 0; k--) { ires[k] = rem % dim; rem /= dim; }
        int new_mu = ires[rC];
        double sum = 0.0;
        for (int sig = 0; sig < dim; sig++) {
            int isrc[DG_MAX_RANK], si = 0;
            for (int p = 0; p < t->rank_r; p++) {
                if (p == contra_pos) isrc[si++] = sig;
                else isrc[si++] = ires[(p < contra_pos) ? p : p - 1];
            }
            for (int p = 0; p < t->rank_s; p++)
                isrc[si++] = ires[rC + 1 + p];
            int fi = tensor_flat_index(t, isrc);
            double tv = (fi >= 0) ? t->components[fi] : 0.0;
            int gi[2] = {sig, new_mu};
            sum += tensor_get(g, gi) * tv;
        }
        c->components[fc] = sum;
    }
    return c;
}

Tensor* tensor_transform(const Tensor *t, const double *J, const double *Jinv,
                          int new_chart) {
    if (!t || !J || !Jinv) return NULL;
    int r = t->rank_r, s = t->rank_s, dim = t->dim;
    Tensor *res = tensor_alloc(r, s, dim);
    if (!res) return NULL;
    res->chart_index = new_chart;
    int ncomp = tensor_ncomp(r, s, dim);
    int ncomp_old = ncomp;
    for (int fc = 0; fc < ncomp; fc++) {
        int inew[DG_MAX_RANK], iold[DG_MAX_RANK];
        int rem = fc;
        for (int k = r + s - 1; k >= 0; k--) { inew[k] = rem % dim; rem /= dim; }
        double sum = 0.0;
        for (int fo = 0; fo < ncomp_old; fo++) {
            int rem_o = fo;
            for (int k = r + s - 1; k >= 0; k--) { iold[k] = rem_o % dim; rem_o /= dim; }
            double coeff = 1.0;
            for (int a = 0; a < r; a++)
                coeff *= J[inew[a] * dim + iold[a]];
            for (int b = 0; b < s; b++)
                coeff *= Jinv[iold[r + b] * dim + inew[r + b]];
            int fio = tensor_flat_index(t, iold);
            if (fio >= 0 && coeff != 0.0)
                sum += coeff * t->components[fio];
        }
        res->components[fc] = sum;
    }
    return res;
}

Tensor* tensor_symmetrize_02(const Tensor *t) {
    if (!t || t->rank_r != 0 || t->rank_s != 2) return NULL;
    int dim = t->dim;
    Tensor *s = tensor_alloc(0, 2, dim);
    if (!s) return NULL;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++) {
            int i1[2] = {mu, nu}, i2[2] = {nu, mu};
            tensor_set(s, i1, 0.5 * (tensor_get(t, i1) + tensor_get(t, i2)));
        }
    return s;
}

Tensor* tensor_antisymmetrize_02(const Tensor *t) {
    if (!t || t->rank_r != 0 || t->rank_s != 2) return NULL;
    int dim = t->dim;
    Tensor *a = tensor_alloc(0, 2, dim);
    if (!a) return NULL;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++) {
            int i1[2] = {mu, nu}, i2[2] = {nu, mu};
            tensor_set(a, i1, 0.5 * (tensor_get(t, i1) - tensor_get(t, i2)));
        }
    return a;
}

int tensor_is_symmetric_02(const Tensor *t) {
    if (!t || t->rank_r != 0 || t->rank_s != 2) return 0;
    int dim = t->dim;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = mu + 1; nu < dim; nu++) {
            int i1[2] = {mu, nu}, i2[2] = {nu, mu};
            if (fabs(tensor_get(t, i1) - tensor_get(t, i2)) > 1e-12) return 0;
        }
    return 1;
}

int tensor_is_antisymmetric_02(const Tensor *t) {
    if (!t || t->rank_r != 0 || t->rank_s != 2) return 0;
    int dim = t->dim;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = mu; nu < dim; nu++) {
            int i1[2] = {mu, nu}, i2[2] = {nu, mu};
            if (mu == nu && fabs(tensor_get(t, i1)) > 1e-12) return 0;
            if (fabs(tensor_get(t, i1) + tensor_get(t, i2)) > 1e-12) return 0;
        }
    return 1;
}

double tensor_norm(const Tensor *t) {
    if (!t) return 0.0;
    double sum = 0.0;
    int n = tensor_ncomp(t->rank_r, t->rank_s, t->dim);
    for (int i = 0; i < n; i++)
        sum += t->components[i] * t->components[i];
    return sqrt(sum);
}

double tensor_trace_11(const Tensor *t) {
    if (!t || t->rank_r != 1 || t->rank_s != 1) return 0.0;
    double tr = 0.0;
    for (int mu = 0; mu < t->dim; mu++) {
        int idx[2] = {mu, mu};
        tr += tensor_get(t, idx);
    }
    return tr;
}

void tensor_print(const Tensor *t) {
    if (!t) { printf("Tensor: (null)\n"); return; }
    int n = tensor_ncomp(t->rank_r, t->rank_s, t->dim);
    printf("Tensor type (%d,%d) dim=%d total=%d:\n", t->rank_r, t->rank_s, t->dim, n);
    int printed = 0;
    for (int i = 0; i < n && printed < 80; i++) {
        if (fabs(t->components[i]) > 1e-14) {
            int idx[DG_MAX_RANK], rem = i, rank = t->rank_r + t->rank_s;
            for (int k = rank - 1; k >= 0; k--) { idx[k] = rem % t->dim; rem /= t->dim; }
            printf("  T[");
            for (int k = 0; k < rank; k++) printf("%d%s", idx[k], k < rank-1 ? "," : "");
            printf("] = %g\n", t->components[i]);
            printed++;
        }
    }
    if (printed >= 80) printf("  ... (%d more non-zero)\n", n - printed);
}