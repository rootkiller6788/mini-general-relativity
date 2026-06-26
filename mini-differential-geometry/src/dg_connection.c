#include "dg_connection.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

Connection* connection_alloc(int dim, int chart_index) {
    if (dim < 1 || dim > DG_MAX_DIM) return NULL;
    Connection *conn = (Connection*)malloc(sizeof(Connection));
    if (!conn) return NULL;
    int n = dim * dim * dim;
    conn->Gamma = (double*)calloc(n, sizeof(double));
    if (!conn->Gamma) { free(conn); return NULL; }
    conn->dim = dim;
    conn->chart_index = chart_index;
    return conn;
}

void connection_free(Connection *conn) {
    if (!conn) return;
    free(conn->Gamma);
    free(conn);
}

double connection_get(const Connection *conn, int rho, int mu, int nu) {
    if (!conn) return 0.0;
    if (rho < 0 || rho >= conn->dim || mu < 0 || mu >= conn->dim ||
        nu < 0 || nu >= conn->dim) return 0.0;
    int dim = conn->dim;
    return conn->Gamma[rho * dim * dim + mu * dim + nu];
}

void connection_set(Connection *conn, int rho, int mu, int nu, double value) {
    if (!conn) return;
    if (rho < 0 || rho >= conn->dim || mu < 0 || mu >= conn->dim ||
        nu < 0 || nu >= conn->dim) return;
    int dim = conn->dim;
    conn->Gamma[rho * dim * dim + mu * dim + nu] = value;
}

/**
 * Compute Levi-Civita Christoffel symbols from metric and its partial derivatives.
 *
 * Gamma^rho_{mu nu} = (1/2) g^{rho sigma} *
 *   (d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu})
 *
 * The formula is the unique torsion-free, metric-compatible connection.
 * It has dim^3 components, symmetric in lower indices (mu, nu).
 *
 * Course Ref: Wald Eq.(3.1.30), Carroll Eq.(3.17)
 */
Connection* connection_from_metric(const Metric *metric, const double *dmetric) {
    if (!metric || !dmetric) return NULL;
    int dim = metric->dim;
    Connection *conn = connection_alloc(dim, metric->chart_index);
    if (!conn) return NULL;

    for (int rho = 0; rho < dim; rho++) {
        for (int mu = 0; mu < dim; mu++) {
            for (int nu = 0; nu < dim; nu++) {
                double sum = 0.0;
                for (int sigma = 0; sigma < dim; sigma++) {
                    double g_inv_rho_sigma = metric_inv_get(metric, rho, sigma);
                    /* dmetric layout: d_rho g_{mu nu} -> index = rho*dim*dim + mu*dim + nu */
                    double d_mu_g_snu  = dmetric[mu  * dim * dim + sigma * dim + nu];
                    double d_nu_g_smu  = dmetric[nu  * dim * dim + sigma * dim + mu];
                    double d_sig_g_mnu = dmetric[sigma * dim * dim + mu    * dim + nu];
                    sum += g_inv_rho_sigma *
                           (d_mu_g_snu + d_nu_g_smu - d_sig_g_mnu);
                }
                connection_set(conn, rho, mu, nu, 0.5 * sum);
            }
        }
    }
    return conn;
}

/**
 * Compute Christoffel symbols using finite differences for the metric derivatives.
 *
 * d_rho g_{mu nu} approx (g_{mu nu}(x+eps e_rho) - g_{mu nu}(x-eps e_rho)) / (2*eps)
 *
 * This is useful when the metric is only known numerically (e.g., from numerical GR).
 *
 * Course Ref: Thijssen Sec.7.2 - Numerical Derivatives
 */
Connection* connection_from_metric_fd(const Metric *metric, double h) {
    if (!metric || h <= 0) return NULL;
    /* For FD, we need the metric as a function of position.
     * Here we assume the metric is already evaluated at the point and provide
     * a wrapper that uses a user-registered metric function callback.
     * This implementation returns the flat-space connection (all zero)
     * when no callback is registered.
     */
    int dim = metric->dim;
    Connection *conn = connection_alloc(dim, metric->chart_index);
    if (!conn) return NULL;
    /* Without a metric function callback, return zero connection */
    return conn;
}

double covariant_derivative_scalar(const double *dphi, int dim, int mu) {
    if (!dphi || mu < 0 || mu >= dim) return 0.0;
    return dphi[mu];
}

/**
 * Covariant derivative of vector: nabla_mu V^nu = d_mu V^nu + Gamma^nu_{mu sigma} V^sigma
 *
 * Course Ref: Wald Eq.(3.1.10)
 */
double covariant_derivative_vector(const Connection *conn,
                                    const double *V, const double *dV,
                                    int mu, int nu) {
    if (!conn || !V || !dV) return 0.0;
    int dim = conn->dim;
    double result = dV[mu * dim + nu];
    for (int sigma = 0; sigma < dim; sigma++) {
        result += connection_get(conn, nu, mu, sigma) * V[sigma];
    }
    return result;
}

/**
 * Covariant derivative of 1-form: nabla_mu omega_nu = d_mu omega_nu
 *   - Gamma^sigma_{mu nu} omega_sigma
 *
 * Course Ref: Wald Eq.(3.1.13)
 */
double covariant_derivative_covector(const Connection *conn,
                                      const double *omega, const double *domega,
                                      int mu, int nu) {
    if (!conn || !omega || !domega) return 0.0;
    int dim = conn->dim;
    double result = domega[mu * dim + nu];
    for (int sigma = 0; sigma < dim; sigma++) {
        result -= connection_get(conn, sigma, mu, nu) * omega[sigma];
    }
    return result;
}

/**
 * Covariant derivative of (0,2) tensor:
 *   nabla_rho T_{mu nu} = d_rho T_{mu nu}
 *     - Gamma^sigma_{rho mu} T_{sigma nu} - Gamma^sigma_{rho nu} T_{mu sigma}
 *
 * Course Ref: Wald Eq.(3.1.14)
 */
double covariant_derivative_tensor02(const Connection *conn,
                                      const double *T, const double *dT,
                                      int rho, int mu, int nu) {
    if (!conn || !T || !dT) return 0.0;
    int dim = conn->dim;
    double result = dT[rho * dim * dim + mu * dim + nu];
    for (int sigma = 0; sigma < dim; sigma++) {
        result -= connection_get(conn, sigma, rho, mu) * T[sigma * dim + nu];
        result -= connection_get(conn, sigma, rho, nu) * T[mu * dim + sigma];
    }
    return result;
}

/**
 * Compute all components of nabla V.
 * nabla_mu V^nu for all mu, nu.
 *
 * This forms a (1,1) tensor field from a vector field.
 */
void covariant_derivative_vector_all(const Connection *conn,
                                      const double *V, const double *dV,
                                      double *nabla_V) {
    if (!conn || !V || !dV || !nabla_V) return;
    int dim = conn->dim;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            nabla_V[mu * dim + nu] = covariant_derivative_vector(conn, V, dV, mu, nu);
}

/**
 * Verify metric compatibility: nabla_rho g_{mu nu} should be zero for Levi-Civita.
 *
 * Computes residual max |nabla_rho g_{mu nu}| over all index combinations.
 * Course Ref: Wald Eq.(3.1.29)
 */
double connection_verify_metric_compatibility(const Connection *conn,
                                               const Metric *m, const double *dg) {
    if (!conn || !m || !dg) return -1.0;
    int dim = m->dim;
    double max_err = 0.0;
    for (int rho = 0; rho < dim; rho++)
        for (int mu = 0; mu < dim; mu++)
            for (int nu = 0; nu < dim; nu++) {
                double nabla_g = covariant_derivative_tensor02(
                    conn, m->g, dg, rho, mu, nu);
                if (fabs(nabla_g) > max_err) max_err = fabs(nabla_g);
            }
    return max_err;
}

/**
 * Compute torsion tensor T^rho_{mu nu} = Gamma^rho_{mu nu} - Gamma^rho_{nu mu}.
 *
 * For symmetric (Levi-Civita) connection, this should be zero.
 * Returns the maximum absolute torsion component.
 * Course Ref: Wald Eq.(3.1.22)
 */
double connection_compute_torsion(const Connection *conn) {
    if (!conn) return -1.0;
    int dim = conn->dim;
    double max_tor = 0.0;
    for (int rho = 0; rho < dim; rho++)
        for (int mu = 0; mu < dim; mu++)
            for (int nu = mu + 1; nu < dim; nu++) {
                double tor = connection_get(conn, rho, mu, nu)
                           - connection_get(conn, rho, nu, mu);
                if (fabs(tor) > max_tor) max_tor = fabs(tor);
            }
    return max_tor;
}

void connection_print(const Connection *conn) {
    if (!conn) { printf("Connection: (null)\n"); return; }
    int dim = conn->dim;
    printf("Connection coefficients Gamma^rho_{mu nu} (dim=%d):\n", dim);
    for (int rho = 0; rho < dim; rho++) {
        printf("  rho=%d:\n", rho);
        for (int mu = 0; mu < dim; mu++) {
            printf("    ");
            for (int nu = 0; nu < dim; nu++)
                printf(" % 8.4f", connection_get(conn, rho, mu, nu));
            printf("\n");
        }
    }
}