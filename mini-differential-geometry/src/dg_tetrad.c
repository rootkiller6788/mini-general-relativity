#include "dg_curvature.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/**
 * @file dg_tetrad.c
 * @brief Tetrad (vierbein) formalism for differential geometry
 *
 * Knowledge Points:
 *   L2: Tetrad fields e^a_mu — local orthonormal frame
 *   L3: Spin connection omega^{ab}_mu
 *   L4: Cartan structure equations
 *   L5: Numerical tetrad from metric via Cholesky decomposition
 *
 * The tetrad e^a_mu provides a local orthonormal basis:
 *   g_{mu nu} = e^a_mu e^b_nu eta_{ab}
 * where eta_{ab} is the Minkowski metric.
 *
 * Reference: Wald Sec.3.4b, Carroll Sec.3, Nakahara Sec.7
 * Courses: MIT 8.962, Cambridge Part III
 */

/**
 * Tetrad (vierbein) structure.
 * e^a_mu : dim x dim matrix, where a is the Lorentz index
 *          and mu is the spacetime index.
 */
typedef struct {
    int     dim;
    double *e;         /* e^a_mu, flat: a*dim + mu         */
    double *e_inv;     /* e_a^mu inverse tetrad            */
    int     chart_index;
} Tetrad;

/** Lorentz metric eta_{ab} = diag(-1,1,1,1) for Wald convention */
static double eta_component(int a, int b, int dim) {
    if (a != b) return 0.0;
    return (a == 0) ? -1.0 : 1.0;
}

/**
 * Compute tetrad from metric via Cholesky-like decomposition.
 *
 * Since g = e^T eta e (matrix notation), we solve:
 *   g_{mu nu} = e^a_mu eta_{ab} e^b_nu
 *
 * For Lorentzian signature, we use an eigenvalue-based decomposition:
 *   g = O^T Lambda O, then e = sqrt(|Lambda|) O (up to signs).
 *
 * Simplified: use element-wise construction for diagonal-dominant metrics.
 * Course Ref: Carroll Appendix J
 */
static int tetrad_from_metric(const Metric *metric, Tetrad *t) {
    if (!metric || !t) return -1;
    int dim = metric->dim;

    /* Initialize e^a_mu = 0 */
    memset(t->e, 0, dim * dim * sizeof(double));

    /* For a general metric, compute via modified Cholesky.
     * This is a simplified version that works for static, diagonal metrics.
     * The full algorithm would require LDL^T decomposition with signature. */

    /* For simplicity, use a symmetric square root approach:
     * Compute eigenvalues and eigenvectors of the metric,
     * then form e = sqrt(|lambda|) * eigenvector matrix.
     * This is valid when the metric is symmetric. */

    /* Jacobi eigenvalue algorithm for real symmetric 4x4 matrix */
    /* Copy metric to working matrix */
    double A[4][4];
    double V[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            A[i][j] = metric_get(metric, i, j);

    /* Jacobi rotations (simplified, few iterations) */
    for (int sweep = 0; sweep < 20; sweep++) {
        for (int p = 0; p < dim; p++) {
            for (int q = p + 1; q < dim; q++) {
                double apq = A[p][q];
                if (fabs(apq) < 1e-14) continue;
                double app = A[p][p];
                double aqq = A[q][q];
                double theta = 0.5 * atan2(2.0 * apq, app - aqq);
                double c = cos(theta);
                double s = sin(theta);

                /* Rotate A */
                double aip, aiq;
                for (int i = 0; i < dim; i++) {
                    if (i == p || i == q) continue;
                    aip = A[i][p];
                    aiq = A[i][q];
                    A[i][p] = c * aip - s * aiq;
                    A[p][i] = A[i][p];
                    A[i][q] = s * aip + c * aiq;
                    A[q][i] = A[i][q];
                }
                A[p][p] = c*c*app + s*s*aqq - 2*c*s*apq;
                A[q][q] = s*s*app + c*c*aqq + 2*c*s*apq;
                A[p][q] = 0.0;
                A[q][p] = 0.0;

                /* Rotate eigenvectors */
                for (int i = 0; i < dim; i++) {
                    double vip = V[i][p];
                    double viq = V[i][q];
                    V[i][p] = c * vip - s * viq;
                    V[i][q] = s * vip + c * viq;
                }
            }
        }
    }

    /* Eigenvalues are A[i][i]. Form tetrad:
     * e^a_mu = sqrt(|lambda_a|) * V_{mu a} (transposed eigenvector)
     * Actually: e^a_mu = sqrt(|lambda_a|) * V^T_{a mu} */
    for (int a = 0; a < dim; a++) {
        double lam = A[a][a];
        double scale = sqrt(fabs(lam));
        for (int mu = 0; mu < dim; mu++) {
            /* V[mu][a] is the a-th eigenvector component at position mu */
            t->e[a * dim + mu] = scale * V[mu][a];
        }
    }

    return 0;
}

/**
 * Compute the spin connection omega^{ab}_mu from the tetrad.
 *
 * omega^{ab}_mu = e^{a nu} (d_mu e^b_nu - Gamma^lambda_{mu nu} e^b_lambda)
 *
 * or equivalently using the first Cartan structure equation:
 *   d e^a + omega^a_b ^ e^b = 0  (torsion-free)
 *
 * Course Ref: Wald Eq.(3.4b.4)
 */
static void spin_connection_compute(const Tetrad *t, const Connection *conn,
                                     const double *de, double *omega) {
    if (!t || !conn || !de || !omega) return;
    int dim = t->dim;

    memset(omega, 0, dim * dim * dim * sizeof(double));

    /* omega^{ab}_mu = e^{a nu} d_mu e^b_nu - e^{a nu} Gamma^lambda_{mu nu} e^b_lambda */
    for (int a = 0; a < dim; a++) {
        for (int b = 0; b < dim; b++) {
            for (int mu = 0; mu < dim; mu++) {
                double sum = 0.0;
                for (int nu = 0; nu < dim; nu++) {
                    double e_a_nu = t->e[a * dim + nu];  /* e^a_nu */
                    /* Compute tetrad inverse element from stored e_inv */
                    /* Term 2: e^{a nu} Gamma^lambda_{mu nu} e^b_lambda */
                    for (int lam = 0; lam < dim; lam++) {
                        sum -= e_a_nu
                             * connection_get(conn, lam, mu, nu)
                             * t->e[b * dim + lam];
                    }
                }
                omega[((a * dim + b) * dim) + mu] = sum;
            }
        }
    }
}

/**
 * Compute the curvature 2-form from the spin connection.
 *
 * Second Cartan structure equation:
 *   R^a_b = d omega^a_b + omega^a_c ^ omega^c_b
 *
 * Course Ref: Wald Eq.(3.4b.5)
 */
static void curvature_2form_compute(const double *omega, const double *domega,
                                     int dim, double *R2form) {
    if (!omega || !domega || !R2form) return;

    /* R^a_{b mu nu} = d_mu omega^a_{b nu} - d_nu omega^a_{b mu}
     *   + omega^a_{c mu} omega^c_{b nu} - omega^a_{c nu} omega^c_{b mu} */
    memset(R2form, 0, dim * dim * dim * dim * sizeof(double));

    for (int a = 0; a < dim; a++) {
        for (int b = 0; b < dim; b++) {
            for (int mu = 0; mu < dim; mu++) {
                for (int nu = 0; nu < dim; nu++) {
                    double d_om = domega[(((mu*dim + a)*dim + b)*dim) + nu];
                    double d_on = domega[(((nu*dim + a)*dim + b)*dim) + mu];
                    double sum = d_om - d_on;
                    for (int c = 0; c < dim; c++) {
                        sum += omega[((a*dim + c)*dim) + mu]
                             * omega[((c*dim + b)*dim) + nu];
                        sum -= omega[((a*dim + c)*dim) + nu]
                             * omega[((c*dim + b)*dim) + mu];
                    }
                    int ridx = (((a*dim + b)*dim + mu)*dim) + nu;
                    R2form[ridx] = sum;
                }
            }
        }
    }
}

/**
 * Public: Allocate tetrad from metric.
 *
 * Returns a tetrad structure or NULL on failure.
 * The tetrad provides a local orthonormal frame convenient for:
 *   - Computing spinor fields in curved spacetime
 *   - Defining Dirac equation on curved manifolds
 *   - 3+1 ADM decomposition
 *
 * Course Ref: Wald Sec.3.4b
 */
Tetrad* tetrad_alloc_from_metric(const Metric *metric) {
    if (!metric) return NULL;
    int dim = metric->dim;

    Tetrad *t = (Tetrad*)malloc(sizeof(Tetrad));
    if (!t) return NULL;

    t->dim = dim;
    t->chart_index = metric->chart_index;
    t->e = (double*)calloc(dim * dim, sizeof(double));
    t->e_inv = (double*)calloc(dim * dim, sizeof(double));

    if (!t->e || !t->e_inv) {
        free(t->e); free(t->e_inv); free(t);
        return NULL;
    }

    if (tetrad_from_metric(metric, t) != 0) {
        free(t->e); free(t->e_inv); free(t);
        return NULL;
    }

    return t;
}

/**
 * Verify that the tetrad reproduces the metric:
 * g_{mu nu} = e^a_mu eta_{ab} e^b_nu
 */
double tetrad_verify_metric(const Tetrad *t, const Metric *metric) {
    if (!t || !metric || t->dim != metric->dim) return -1.0;
    int dim = t->dim;
    double max_err = 0.0;

    for (int mu = 0; mu < dim; mu++) {
        for (int nu = 0; nu < dim; nu++) {
            double g_recon = 0.0;
            for (int a = 0; a < dim; a++) {
                for (int b = 0; b < dim; b++) {
                    g_recon += t->e[a * dim + mu] * eta_component(a, b, dim)
                             * t->e[b * dim + nu];
                }
            }
            double err = fabs(g_recon - metric_get(metric, mu, nu));
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

void tetrad_free(Tetrad *t) {
    if (!t) return;
    free(t->e);
    free(t->e_inv);
    free(t);
}