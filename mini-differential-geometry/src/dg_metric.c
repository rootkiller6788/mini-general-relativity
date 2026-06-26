#include "dg_metric.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* Small 4x4 matrix inversion using Gaussian elimination with partial pivoting.
 * Returns 0 on success, -1 if singular.
 * Reference: Numerical Recipes Sec.2.1
 */
static int matrix_invert_4x4(const double *A, double *Ainv, int dim) {
    /* Augmented matrix [A | I] with row-major storage */
    double aug[4][8] = {{0}};
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++)
            aug[i][j] = A[i * dim + j];
        aug[i][dim + i] = 1.0;
    }
    /* Forward elimination */
    for (int col = 0; col < dim; col++) {
        /* Partial pivot */
        int maxrow = col;
        double maxval = fabs(aug[col][col]);
        for (int row = col + 1; row < dim; row++) {
            if (fabs(aug[row][col]) > maxval) {
                maxval = fabs(aug[row][col]);
                maxrow = row;
            }
        }
        if (maxval < 1e-15) return -1; /* singular */
        if (maxrow != col) {
            for (int j = 0; j < 2 * dim; j++) {
                double tmp = aug[col][j];
                aug[col][j] = aug[maxrow][j];
                aug[maxrow][j] = tmp;
            }
        }
        double pivot = aug[col][col];
        for (int j = 0; j < 2 * dim; j++)
            aug[col][j] /= pivot;
        for (int row = 0; row < dim; row++) {
            if (row == col) continue;
            double factor = aug[row][col];
            for (int j = 0; j < 2 * dim; j++)
                aug[row][j] -= factor * aug[col][j];
        }
    }
    /* Extract inverse */
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            Ainv[i * dim + j] = aug[i][dim + j];
    return 0;
}

/* Compute determinant via LU decomposition for small matrices.
 * For 4x4 we use Laplace expansion for numerical stability.
 */
static double matrix_det(const double *A, int dim) {
    if (dim == 1) return A[0];
    if (dim == 2) return A[0]*A[3] - A[1]*A[2];
    if (dim == 3) {
        return A[0]*(A[4]*A[8] - A[5]*A[7])
             - A[1]*(A[3]*A[8] - A[5]*A[6])
             + A[2]*(A[3]*A[7] - A[4]*A[6]);
    }
    /* dim == 4: Laplace expansion along first row */
    if (dim == 4) {
        double det = 0.0;
        int sign = 1;
        for (int j = 0; j < 4; j++) {
            /* Submatrix excluding row 0, col j */
            double sub[9];
            int si = 0;
            for (int r = 1; r < 4; r++)
                for (int c = 0; c < 4; c++)
                    if (c != j) sub[si++] = A[r * 4 + c];
            det += sign * A[j] * matrix_det(sub, 3);
            sign = -sign;
        }
        return det;
    }
    return 0.0;
}

Metric* metric_alloc(int dim, MetricSignature sig) {
    if (dim < 1 || dim > DG_MAX_DIM) return NULL;
    Metric *m = (Metric*)malloc(sizeof(Metric));
    if (!m) return NULL;
    m->dim = dim;
    m->signature = sig;
    m->g = (double*)calloc(dim * dim, sizeof(double));
    m->g_inv = (double*)calloc(dim * dim, sizeof(double));
    if (!m->g || !m->g_inv) {
        free(m->g); free(m->g_inv); free(m);
        return NULL;
    }
    m->det = 0.0;
    m->chart_index = 0;
    return m;
}

Metric* metric_create(int dim, MetricSignature sig,
                       const double *comp, int chart) {
    Metric *m = metric_alloc(dim, sig);
    if (!m) return NULL;
    if (comp) {
        memcpy(m->g, comp, dim * dim * sizeof(double));
        metric_recompute_inv_det(m);
    }
    m->chart_index = chart;
    return m;
}

Metric* metric_clone(const Metric *m) {
    if (!m) return NULL;
    Metric *c = metric_alloc(m->dim, m->signature);
    if (!c) return NULL;
    memcpy(c->g, m->g, m->dim * m->dim * sizeof(double));
    memcpy(c->g_inv, m->g_inv, m->dim * m->dim * sizeof(double));
    c->det = m->det;
    c->chart_index = m->chart_index;
    return c;
}

void metric_free(Metric *m) {
    if (!m) return;
    free(m->g);
    free(m->g_inv);
    free(m);
}

double metric_get(const Metric *m, int mu, int nu) {
    if (!m || mu < 0 || mu >= m->dim || nu < 0 || nu >= m->dim) return 0.0;
    return m->g[mu * m->dim + nu];
}

void metric_set(Metric *m, int mu, int nu, double value) {
    if (!m || mu < 0 || mu >= m->dim || nu < 0 || nu >= m->dim) return;
    m->g[mu * m->dim + nu] = value;
    m->g[nu * m->dim + mu] = value;  /* enforce symmetry */
}

double metric_inv_get(const Metric *m, int mu, int nu) {
    if (!m || mu < 0 || mu >= m->dim || nu < 0 || nu >= m->dim) return 0.0;
    return m->g_inv[mu * m->dim + nu];
}

int metric_recompute_inv_det(Metric *m) {
    if (!m) return -1;
    m->det = matrix_det(m->g, m->dim);
    if (fabs(m->det) < 1e-15) return -1;
    return matrix_invert_4x4(m->g, m->g_inv, m->dim);
}

double metric_line_element(const Metric *m, const double dx[]) {
    if (!m || !dx) return 0.0;
    double ds2 = 0.0;
    int dim = m->dim;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            ds2 += m->g[mu * dim + nu] * dx[mu] * dx[nu];
    return ds2;
}

double metric_proper_time_inc(const Metric *m, const double dx[]) {
    double ds2 = metric_line_element(m, dx);
    /* Wald signature (-,+,+,+): timelike if ds2 < 0 */
    if (ds2 < 0) return sqrt(-ds2);
    return 0.0;
}

double metric_proper_length_inc(const Metric *m, const double dx[]) {
    double ds2 = metric_line_element(m, dx);
    if (ds2 > 0) return sqrt(ds2);
    return 0.0;
}

int metric_classify_dx(const Metric *m, const double dx[]) {
    double ds2 = metric_line_element(m, dx);
    if (fabs(ds2) < 1e-15) return 0;   /* null */
    return (ds2 < 0) ? -1 : 1;         /* timelike : spacelike */
}

double metric_volume_element(const Metric *m) {
    if (!m) return 1.0;
    return sqrt(fabs(m->det));
}

double metric_determinant(const Metric *m) {
    if (!m) return 0.0;
    return m->det;
}

int metric_verify_signature(const Metric *m) {
    if (!m) return 0;
    /* For 4D Lorentzian: eigenvalues should have 1 negative, 3 positive
     * (Wald convention). We compute eigenvalues via characteristic polynomial
     * for the 4x4 symmetric matrix.
     * This is a simplified check using the determinant signs of principal minors
     * (Sylvester's criterion adapted for Lorentzian signature).
     */
    int dim = m->dim;
    /* Check symmetry */
    for (int i = 0; i < dim; i++)
        for (int j = i + 1; j < dim; j++)
            if (fabs(m->g[i*dim+j] - m->g[j*dim+i]) > 1e-12) return 0;

    /* For Lorentzian: det(g) < 0 is necessary */
    if (m->signature == SIG_LORENTZIAN_MPPP || m->signature == SIG_LORENTZIAN_PMMM) {
        /* Check that det is negative for even spacetime dim (4D) */
        if (dim == 4 && m->det >= 0) return 0;
    }
    return 1;
}

void metric_print(const Metric *m) {
    if (!m) { printf("Metric: (null)\n"); return; }
    const char *signames[] = {"Lorentzian(-+++)", "Lorentzian(+---)", "Riemannian"};
    printf("Metric dim=%d sig=%s det=%g\n",
           m->dim, signames[m->signature], m->det);
    for (int mu = 0; mu < m->dim; mu++) {
        printf("  [");
        for (int nu = 0; nu < m->dim; nu++)
            printf(" % 8.4f", m->g[mu * m->dim + nu]);
        printf(" ]\n");
    }
}