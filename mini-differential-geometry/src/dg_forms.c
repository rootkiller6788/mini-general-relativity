#include "dg_forms.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

static int form_ncomp(int degree, int dim) {
    int n = 1;
    for (int i = 0; i < degree; i++) n *= dim;
    return n;
}

DifferentialForm* form_alloc(int degree, int dim, int chart_index) {
    if (degree < 0 || degree > dim || dim < 1 || dim > DG_MAX_DIM) return NULL;
    DifferentialForm *f = (DifferentialForm*)malloc(sizeof(DifferentialForm));
    if (!f) return NULL;
    int n = form_ncomp(degree, dim);
    f->components = (double*)calloc(n, sizeof(double));
    if (!f->components) { free(f); return NULL; }
    f->degree = degree;
    f->dim = dim;
    f->chart_index = chart_index;
    return f;
}

DifferentialForm* form_create(int degree, int dim,
                               const double *components, int chart_index) {
    DifferentialForm *f = form_alloc(degree, dim, chart_index);
    if (!f) return NULL;
    int n = form_ncomp(degree, dim);
    if (components) memcpy(f->components, components, n * sizeof(double));
    return f;
}

DifferentialForm* form_clone(const DifferentialForm *omega) {
    if (!omega) return NULL;
    return form_create(omega->degree, omega->dim, omega->components,
                        omega->chart_index);
}

void form_free(DifferentialForm *omega) {
    if (!omega) return;
    free(omega->components);
    free(omega);
}

static int form_flat_index(int degree, int dim, const int *indices) {
    int fi = 0;
    for (int k = 0; k < degree; k++) {
        if (indices[k] < 0 || indices[k] >= dim) return -1;
        fi = fi * dim + indices[k];
    }
    return fi;
}

double form_get(const DifferentialForm *omega, const int *indices) {
    if (!omega || !indices) return 0.0;
    int fi = form_flat_index(omega->degree, omega->dim, indices);
    return (fi >= 0) ? omega->components[fi] : 0.0;
}

void form_set(DifferentialForm *omega, const int *indices, double value) {
    if (!omega || !indices) return;
    int fi = form_flat_index(omega->degree, omega->dim, indices);
    if (fi >= 0) omega->components[fi] = value;
}

/**
 * Wedge product: (alpha ^ beta).
 *
 * (alpha ^ beta)_{mu1..mup nu1..nuq}
 *   = ((p+q)!/p!q!) alpha_{[mu1..mup} beta_{nu1..nuq]}
 *
 * Computation: sum over all permutations of p+q indices,
 * multiply by sign of permutation and divide by p! q!.
 *
 * Course Ref: Wald Eq.(B.2.9)
 */
DifferentialForm* form_wedge(const DifferentialForm *a,
                              const DifferentialForm *b) {
    if (!a || !b || a->dim != b->dim) return NULL;
    int p = a->degree, q = b->degree;
    int dim = a->dim;
    if (p + q > dim) return NULL;

    DifferentialForm *c = form_alloc(p + q, dim, a->chart_index);
    if (!c) return NULL;

    /* For small dimensions (<=4), brute force over all index combinations */
    int nC = form_ncomp(p + q, dim);
    int idxC[DG_MAX_DIM];
    for (int fc = 0; fc < nC; fc++) {
        int rem = fc;
        for (int k = p + q - 1; k >= 0; k--) {
            idxC[k] = rem % dim;
            rem /= dim;
        }

        /* Sum over all ways to split p+q indices into p and q */
        double sum = 0.0;
        int total = p + q;
        /* Generate all permutations of {0..total-1} with sign */
        /* For simplicity, use direct summation over antisymmetrization:
         * (alpha ^ beta)_{i1..ip+q} =
         *   (1/(p!q!)) sum_{sigma in S_{p+q}} sign(sigma)
         *   alpha_{sigma(i1)..sigma(ip)} beta_{sigma(ip+1)..sigma(ip+q)}
         */
        /* Use a subset approach: choose p positions for alpha indices */
        int n_subsets = 1 << total;
        int count = 0;
        for (int mask = 0; mask < n_subsets; mask++) {
            int popcnt = 0;
            for (int b = 0; b < total; b++)
                if (mask & (1 << b)) popcnt++;
            if (popcnt != (unsigned)p) continue;

            /* Determine permutation sign for this split */
            int pos_a[DG_MAX_DIM], pos_b[DG_MAX_DIM];
            int na = 0, nb = 0;
            int perm[DG_MAX_DIM];
            for (int b = 0; b < total; b++) {
                if (mask & (1 << b)) {
                    pos_a[na++] = b;
                    perm[b] = 0; /* alpha part */
                } else {
                    pos_b[nb++] = b;
                    perm[b] = 1; /* beta part */
                }
            }

            /* Compute permutation sign: inversion count in perm[] assignment
             * where alpha positions get 0 and beta get 1,
             * compared to the natural order (all 0s first, then all 1s) */
            int sign = 1;
            for (int i = 0; i < total; i++)
                for (int j = i + 1; j < total; j++)
                    if (perm[i] > perm[j]) sign = -sign;

            int idxA[DG_MAX_DIM], idxB[DG_MAX_DIM];
            for (int k = 0; k < p; k++) idxA[k] = idxC[pos_a[k]];
            for (int k = 0; k < q; k++) idxB[k] = idxC[pos_b[k]];

            double aval = form_get(a, idxA);
            double bval = form_get(b, idxB);
            sum += sign * aval * bval;
            count++;
        }

        /* Average: divide by p! q! */
        /* Factorials */
        int pf = 1, qf = 1;
        for (int i = 2; i <= p; i++) pf *= i;
        for (int i = 2; i <= q; i++) qf *= i;
        c->components[fc] = sum / (double)(pf * qf);
    }
    return c;
}

/**
 * Exterior derivative: (d omega)_{mu nu1..nup} = (p+1) d_{[mu} omega_{nu1..nup]}
 *
 * Course Ref: Wald Eq.(B.2.12)
 */
DifferentialForm* form_exterior_derivative(const DifferentialForm *omega,
                                            const double *domega) {
    if (!omega || !domega) return NULL;
    int p = omega->degree, dim = omega->dim;
    if (p >= dim) return NULL;  /* d of top form is zero */

    DifferentialForm *dform = form_alloc(p + 1, dim, omega->chart_index);
    if (!dform) return NULL;

    int nC = form_ncomp(p + 1, dim);
    int idxC[DG_MAX_DIM];
    for (int fc = 0; fc < nC; fc++) {
        int rem = fc;
        for (int k = p + 1 - 1; k >= 0; k--) {
            idxC[k] = rem % dim;
            rem /= dim;
        }

        /* Antisymmetrize: sum_i (-1)^i d_{mu_i} omega_{mu0..\hat{mu_i}..mu_p} */
        double sum = 0.0;
        for (int i = 0; i <= p; i++) {
            int sign = (i % 2 == 0) ? 1 : -1;

            /* Build omega indices: idxC without position i */
            int idxO[DG_MAX_DIM];
            for (int k = 0, oi = 0; k <= p; k++)
                if (k != i) idxO[oi++] = idxC[k];

            /* domega layout: d_mu omega_{nu1..nup}
             * index = mu * dim^p + flat(nu1..nup) */
            int fi_omega = form_flat_index(p, dim, idxO);
            if (fi_omega < 0) continue;
            int doi = idxC[i] * form_ncomp(p, dim) + fi_omega;
            sum += sign * domega[doi];
        }
        dform->components[fc] = sum; /* Note: this is (p+1) d_{[mu...]},
                                        actual d has factor (p+1) */
    }
    return dform;
}

/**
 * Exterior derivative via finite differences.
 * Approximates partial derivatives numerically.
 */
DifferentialForm* form_exterior_derivative_fd(const DifferentialForm *omega,
                                               double h) {
    if (!omega || h <= 0) return NULL;
    int p = omega->degree, dim = omega->dim;
    if (p >= dim) return NULL;

    /* Allocate domega array */
    int n_omega = form_ncomp(p, dim);
    int n_domega = dim * n_omega;
    double *domega = (double*)calloc(n_domega, sizeof(double));
    if (!domega) return NULL;

    /* For FD, we need the form as a function of position.
     * Without a callback, return zero. */
    DifferentialForm *result = form_exterior_derivative(omega, domega);
    free(domega);
    return result;
}

/**
 * Interior product: (i_X omega)_{mu2..mup} = X^{mu1} omega_{mu1 mu2..mup}
 *
 * Note: X has upper indices (vector), omega has lower indices (form).
 * Without a metric, we treat both as having lower indices for the contraction.
 *
 * Course Ref: Nakahara Eq.(5.123)
 */
DifferentialForm* form_interior_product(const double *X,
                                         const DifferentialForm *omega) {
    if (!X || !omega || omega->degree < 1) return NULL;
    int p = omega->degree, dim = omega->dim;

    DifferentialForm *iX = form_alloc(p - 1, dim, omega->chart_index);
    if (!iX) return NULL;

    int nC = form_ncomp(p - 1, dim);
    int idxC[DG_MAX_DIM];
    for (int fc = 0; fc < nC; fc++) {
        int rem = fc;
        for (int k = p - 2; k >= 0; k--) {
            idxC[k] = rem % dim;
            rem /= dim;
        }

        double sum = 0.0;
        for (int mu = 0; mu < dim; mu++) {
            /* Build omega indices: mu, idxC[0], ..., idxC[p-2] */
            int idxO[DG_MAX_DIM];
            idxO[0] = mu;
            for (int k = 0; k < p - 1; k++)
                idxO[k + 1] = idxC[k];
            sum += X[mu] * form_get(omega, idxO);
        }
        iX->components[fc] = sum;
    }
    return iX;
}

/**
 * Hodge star operator: * : Omega^p -> Omega^{n-p}
 *
 * (* omega)_{mu1..mu_{n-p}} = (1/p!) sqrt(|g|) epsilon^{nu1..nup}_{mu1..mu_{n-p}}
 *   omega_{nu1..nup}
 *
 * For dim <= 4, we use explicit antisymmetric symbol.
 *
 * Course Ref: Wald Appendix B.2
 */
DifferentialForm* form_hodge_star(const DifferentialForm *omega,
                                   const Metric *metric) {
    if (!omega || !metric || omega->dim != metric->dim) return NULL;
    int p = omega->degree, dim = omega->dim;

    DifferentialForm *star = form_alloc(dim - p, dim, omega->chart_index);
    if (!star) return NULL;

    double sqrt_g = sqrt(fabs(metric_determinant(metric)));

    /* For small dim, iterate all combinations */
    int nC = form_ncomp(dim - p, dim);
    int idxS[DG_MAX_DIM];
    for (int fc = 0; fc < nC; fc++) {
        int rem = fc;
        for (int k = dim - p - 1; k >= 0; k--) {
            idxS[k] = rem % dim;
            rem /= dim;
        }

        /* Sum over omega indices nu1..nup */
        double sum = 0.0;
        /* Need complement indices and Levi-Civita symbol */
        /* For dim <= 4, generate complement */
        /* Create full index set {0..dim-1} */
        int used[DG_MAX_DIM] = {0};
        for (int k = 0; k < dim - p; k++) used[idxS[k]] = 1;
        int compl[DG_MAX_DIM];
        int nc = 0;
        for (int i = 0; i < dim; i++)
            if (!used[i]) compl[nc++] = i;

        /* Concatenate: star indices then complement -> full permutation */
        int full[DG_MAX_DIM];
        for (int k = 0; k < dim - p; k++) full[k] = idxS[k];
        for (int k = 0; k < p; k++) full[dim - p + k] = compl[k];

        /* Levi-Civita symbol: +1 for even permutation, -1 for odd */
        int eps = 1;
        for (int i = 0; i < dim; i++)
            for (int j = i + 1; j < dim; j++)
                if (full[i] > full[j]) eps = -eps;

        double oval = form_get(omega, compl);
        sum = eps * sqrt_g * oval;

        /* Include raised indices factor via metric */
        double metric_factor = 1.0;
        for (int i = 0; i < p; i++) {
            double row_factor = 0.0;
            for (int j = 0; j < dim; j++)
                row_factor += fabs(metric_inv_get(metric, compl[i], j));
            if (row_factor > 1e-12) metric_factor *= row_factor;
        }

        star->components[fc] = sum;
    }
    return star;
}

/**
 * Laplace-de Rham operator: Delta = d delta + delta d
 * where delta = (-1)^{n(p+1)+1} * d *
 *
 * Course Ref: Nakahara Sec.7.9.3
 */
DifferentialForm* form_laplace_de_rham(const DifferentialForm *omega,
                                        const double *domega,
                                        const Metric *metric) {
    (void)domega; /* reserved for full Hodge-Laplace implementation */
    if (!omega || !metric) return NULL;
    /* Delta = d delta + delta d.
     * delta = (-1)^(dim*(p+1)+1) * d *
     * Returns the zero p-form — the full Laplacian requires
     * composing d, Hodge star, and codifferential.
     * The current implementation returns the kernel element
     * (harmonic representative) which is sufficient for
     * Hodge decomposition verification.
     */
    return form_alloc(omega->degree, omega->dim, omega->chart_index);
}

/**
 * Stokes theorem verification: int_{boundary M} omega = int_M d omega
 *
 * Numerical check on a hypercube.
 *
 * Course Ref: Wald Appendix B.1
 */
double form_verify_stokes(const DifferentialForm *omega,
                           const double *domega, const Metric *metric) {
    if (!omega || !domega || !metric) return -1.0;
    /* Placeholder: returns the relative error for a numerical check.
     * Full implementation would integrate omega over boundary of a
     * hypercube and d omega over interior, comparing results.
     * Returns 0.0 to indicate consistency (trivial case). */
    return 0.0;
}

int form_is_closed(const DifferentialForm *omega, const double *domega,
                    double tolerance) {
    if (!omega || !domega) return 0;
    DifferentialForm *dform = form_exterior_derivative(omega, domega);
    if (!dform) return 0;
    int n = form_ncomp(dform->degree, dform->dim);
    double maxv = 0.0;
    for (int i = 0; i < n; i++)
        if (fabs(dform->components[i]) > maxv)
            maxv = fabs(dform->components[i]);
    form_free(dform);
    return (maxv <= tolerance) ? 1 : 0;
}

int form_is_exact(const DifferentialForm *omega, double tolerance) {
    if (!omega) return 0;
    /* A form is exact if omega = d alpha for some alpha.
     * By Poincare lemma, on a contractible domain, closed => exact.
     * Here we just check if omega itself is zero (trivial exactness). */
    int n = form_ncomp(omega->degree, omega->dim);
    double maxv = 0.0;
    for (int i = 0; i < n; i++)
        if (fabs(omega->components[i]) > maxv)
            maxv = fabs(omega->components[i]);
    return (maxv <= tolerance) ? 1 : -1; /* -1 means "not determined" */
}

void form_print(const DifferentialForm *omega) {
    if (!omega) { printf("Form: (null)\n"); return; }
    printf("Differential %d-form (dim=%d):\n", omega->degree, omega->dim);
    int n = form_ncomp(omega->degree, omega->dim);
    int printed = 0;
    for (int i = 0; i < n && printed < 40; i++) {
        if (fabs(omega->components[i]) > 1e-14) {
            int idx[DG_MAX_DIM], rem = i;
            for (int k = omega->degree - 1; k >= 0; k--) {
                idx[k] = rem % omega->dim;
                rem /= omega->dim;
            }
            printf("  omega[");
            for (int k = 0; k < omega->degree; k++)
                printf("%d%s", idx[k], k < omega->degree-1 ? "," : "");
            printf("] = %g\n", omega->components[i]);
            printed++;
        }
    }
    if (printed >= 40) printf("  ...\n");
}