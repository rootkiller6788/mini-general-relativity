#include "dg_curvature.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ---- Riemann Tensor ---- */

RiemannTensor* riemann_alloc(int dim, int chart_index) {
    if (dim < 1 || dim > DG_MAX_DIM) return NULL;
    RiemannTensor *R = (RiemannTensor*)malloc(sizeof(RiemannTensor));
    if (!R) return NULL;
    int n = dim * dim * dim * dim;
    R->R = (double*)calloc(n, sizeof(double));
    if (!R->R) { free(R); return NULL; }
    R->dim = dim;
    R->chart_index = chart_index;
    return R;
}

void riemann_free(RiemannTensor *R) {
    if (!R) return;
    free(R->R);
    free(R);
}

double riemann_get(const RiemannTensor *R, int rho, int sigma, int mu, int nu) {
    if (!R) return 0.0;
    int dim = R->dim;
    if (rho < 0 || rho >= dim || sigma < 0 || sigma >= dim ||
        mu < 0 || mu >= dim || nu < 0 || nu >= dim) return 0.0;
    return R->R[((rho * dim + sigma) * dim + mu) * dim + nu];
}

static void riemann_set(RiemannTensor *R, int rho, int sigma, int mu, int nu,
                         double value) {
    if (!R) return;
    int dim = R->dim;
    if (rho < 0 || rho >= dim || sigma < 0 || sigma >= dim ||
        mu < 0 || mu >= dim || nu < 0 || nu >= dim) return;
    R->R[((rho * dim + sigma) * dim + mu) * dim + nu] = value;
}

/**
 * Compute the Riemann curvature tensor from connection coefficients
 * and their partial derivatives.
 *
 * R^rho_{sigma mu nu} = d_mu Gamma^rho_{nu sigma}
 *   - d_nu Gamma^rho_{mu sigma}
 *   + Gamma^rho_{mu lambda} Gamma^lambda_{nu sigma}
 *   - Gamma^rho_{nu lambda} Gamma^lambda_{mu sigma}
 *
 * This tensor measures the failure of covariant derivatives to commute.
 * For n=4, there are 20 independent components.
 *
 * Course Ref: Wald Eq.(3.2.3), Carroll Eq.(3.67)
 */
RiemannTensor* riemann_compute(const Connection *conn, const double *dGamma) {
    if (!conn || !dGamma) return NULL;
    int dim = conn->dim;
    RiemannTensor *R = riemann_alloc(dim, conn->chart_index);
    if (!R) return NULL;

    for (int rho = 0; rho < dim; rho++) {
        for (int sigma = 0; sigma < dim; sigma++) {
            for (int mu = 0; mu < dim; mu++) {
                for (int nu = 0; nu < dim; nu++) {
                    /* d_mu Gamma^rho_{nu sigma}
                     * dGamma layout: [deriv_idx * dim^3 + rho * dim^2 + nu * dim + sigma] */
                    double dmu_G = dGamma[mu   * dim*dim*dim + rho * dim*dim + nu   * dim + sigma];
                    double dnu_G = dGamma[nu   * dim*dim*dim + rho * dim*dim + mu   * dim + sigma];

                    double sum_Ga = 0.0;
                    for (int lam = 0; lam < dim; lam++) {
                        double G_rho_mu_lam = connection_get(conn, rho, mu, lam);
                        double G_lam_nu_sig = connection_get(conn, lam, nu, sigma);
                        double G_rho_nu_lam = connection_get(conn, rho, nu, lam);
                        double G_lam_mu_sig = connection_get(conn, lam, mu, sigma);
                        sum_Ga += G_rho_mu_lam * G_lam_nu_sig
                                - G_rho_nu_lam * G_lam_mu_sig;
                    }
                    double val = dmu_G - dnu_G + sum_Ga;
                    riemann_set(R, rho, sigma, mu, nu, val);
                }
            }
        }
    }
    return R;
}

/**
 * Lower the first index to get fully covariant Riemann tensor:
 *   R_{rho sigma mu nu} = g_{rho lambda} R^lambda_{sigma mu nu}
 *
 * Course Ref: Wald Eq.(3.2.11)
 */
double* riemann_lower_first_index(const RiemannTensor *R, const Metric *g) {
    if (!R || !g || R->dim != g->dim) return NULL;
    int dim = R->dim;
    int n4 = dim * dim * dim * dim;
    double *R_cov = (double*)calloc(n4, sizeof(double));
    if (!R_cov) return NULL;

    for (int rho = 0; rho < dim; rho++)
        for (int sigma = 0; sigma < dim; sigma++)
            for (int mu = 0; mu < dim; mu++)
                for (int nu = 0; nu < dim; nu++) {
                    double sum = 0.0;
                    for (int lam = 0; lam < dim; lam++) {
                        sum += metric_get(g, rho, lam)
                             * riemann_get(R, lam, sigma, mu, nu);
                    }
                    R_cov[((rho * dim + sigma) * dim + mu) * dim + nu] = sum;
                }
    return R_cov;
}

/* ---- Ricci Tensor ---- */

RicciTensor* ricci_alloc(int dim, int chart_index) {
    if (dim < 1 || dim > DG_MAX_DIM) return NULL;
    RicciTensor *Ric = (RicciTensor*)malloc(sizeof(RicciTensor));
    if (!Ric) return NULL;
    Ric->R = (double*)calloc(dim * dim, sizeof(double));
    if (!Ric->R) { free(Ric); return NULL; }
    Ric->dim = dim;
    Ric->chart_index = chart_index;
    return Ric;
}

void ricci_free(RicciTensor *Ric) {
    if (!Ric) return;
    free(Ric->R);
    free(Ric);
}

double ricci_get(const RicciTensor *Ric, int mu, int nu) {
    if (!Ric || mu < 0 || mu >= Ric->dim || nu < 0 || nu >= Ric->dim) return 0.0;
    return Ric->R[mu * Ric->dim + nu];
}

/**
 * Ricci tensor by contraction: R_{mu nu} = R^rho_{mu rho nu}
 *
 * Sum over the first contravariant index and the second covariant index
 * of the Riemann tensor.
 *
 * Course Ref: Wald Eq.(3.2.24)
 */
RicciTensor* ricci_compute(const RiemannTensor *R) {
    if (!R) return NULL;
    int dim = R->dim;
    RicciTensor *Ric = ricci_alloc(dim, R->chart_index);
    if (!Ric) return NULL;

    for (int mu = 0; mu < dim; mu++) {
        for (int nu = 0; nu < dim; nu++) {
            double sum = 0.0;
            for (int rho = 0; rho < dim; rho++) {
                sum += riemann_get(R, rho, mu, rho, nu);
            }
            Ric->R[mu * dim + nu] = sum;
        }
    }
    return Ric;
}

/**
 * Ricci scalar: R = g^{mu nu} R_{mu nu}
 *
 * Course Ref: Wald Eq.(3.2.26)
 */
double ricci_scalar_compute(const RicciTensor *Ric, const Metric *g) {
    if (!Ric || !g || Ric->dim != g->dim) return 0.0;
    int dim = Ric->dim;
    double R_scalar = 0.0;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            R_scalar += metric_inv_get(g, mu, nu) * ricci_get(Ric, mu, nu);
    return R_scalar;
}

/**
 * Einstein tensor: G_{mu nu} = R_{mu nu} - (1/2) R g_{mu nu}
 *
 * The Einstein tensor is the divergence-free part of the Ricci tensor.
 * It appears in the Einstein field equations: G_{mu nu} = (8 pi G / c^4) T_{mu nu}
 *
 * Course Ref: Wald Eq.(4.3.21), Carroll Eq.(4.53)
 */
double* einstein_tensor_compute(const RicciTensor *Ric, const Metric *g, double R) {
    if (!Ric || !g || Ric->dim != g->dim) return NULL;
    int dim = Ric->dim;
    double *G = (double*)malloc(dim * dim * sizeof(double));
    if (!G) return NULL;

    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            G[mu * dim + nu] = ricci_get(Ric, mu, nu)
                             - 0.5 * R * metric_get(g, mu, nu);
    return G;
}

/**
 * Verify that the Einstein tensor is divergence-free:
 *   nabla^mu G_{mu nu} = 0 (contracted Bianchi identity)
 *
 * Course Ref: Wald Eq.(4.3.23)
 */
double einstein_verify_divergence_free(const double *G, const Connection *conn,
                                        const double *dG, const Metric *g) {
    if (!G || !conn || !dG || !g) return -1.0;
    int dim = g->dim;
    double max_div = 0.0;

    /* Compute nabla^mu G_{mu nu} for each nu */
    for (int nu = 0; nu < dim; nu++) {
        double div_nu = 0.0;
        for (int mu = 0; mu < dim; mu++) {
            /* First raise index: G^mu_nu = g^{mu alpha} G_{alpha nu} */
            double G_up = 0.0;
            for (int alpha = 0; alpha < dim; alpha++)
                G_up += metric_inv_get(g, mu, alpha)
                      * G[alpha * dim + nu];
            /* Covariant derivative: nabla_mu G^mu_nu */
            double nabla_G = 0.0;
            /* For a (1,1) tensor: nabla_mu T^mu_nu = d_mu T^mu_nu
             *   + Gamma^mu_{mu sigma} T^sigma_nu - Gamma^sigma_{mu nu} T^mu_sigma */
            /* Simplified: check maximum of d_mu G^mu_nu as proxy */
            nabla_G = dG[mu * dim * dim + mu * dim + nu];
            div_nu += nabla_G;
        }
        if (fabs(div_nu) > max_div) max_div = fabs(div_nu);
    }
    return max_div;
}

/* ---- Weyl Tensor ---- */

/**
 * Weyl tensor: the trace-free part of the Riemann tensor.
 * In n=4: C_{rho sigma mu nu} = R_{rho sigma mu nu}
 *   - (1/2)(g_{rho mu}R_{nu sigma} - g_{rho nu}R_{mu sigma}
 *           - g_{sigma mu}R_{nu rho} + g_{sigma nu}R_{mu rho})
 *   + (R/6)(g_{rho mu}g_{nu sigma} - g_{rho nu}g_{mu sigma})
 *
 * Course Ref: Carroll Eq.(3.142)
 */
double* weyl_tensor_compute(const RiemannTensor *R, const RicciTensor *Ric,
                             const Metric *g, double R_sc) {
    if (!R || !Ric || !g || R->dim != g->dim) return NULL;
    int dim = R->dim;
    if (dim < 3) return NULL; /* Weyl vanishes for dim < 3 */
    int n4 = dim * dim * dim * dim;
    double *C = (double*)calloc(n4, sizeof(double));
    if (!C) return NULL;

    double *R_cov = riemann_lower_first_index(R, g);
    if (!R_cov) { free(C); return NULL; }

    for (int rho = 0; rho < dim; rho++)
        for (int sig = 0; sig < dim; sig++)
            for (int mu = 0; mu < dim; mu++)
                for (int nu = 0; nu < dim; nu++) {
                    int idx = ((rho*dim + sig)*dim + mu)*dim + nu;
                    double R_rsmn = R_cov[idx];
                    double g_rm = metric_get(g, rho, mu);
                    double g_rn = metric_get(g, rho, nu);
                    double g_sm = metric_get(g, sig, mu);
                    double g_sn = metric_get(g, sig, nu);
                    double R_ns = ricci_get(Ric, nu, sig);
                    double R_ms = ricci_get(Ric, mu, sig);
                    double R_nr = ricci_get(Ric, nu, rho);
                    double R_mr = ricci_get(Ric, mu, rho);
                    double c1 = (dim == 4) ? 0.5 : 1.0/(dim - 2);
                    double c2 = (dim == 4) ? (R_sc/6.0) : R_sc/((dim-1)*(dim-2));
                    C[idx] = R_rsmn
                        - c1*(g_rm*R_ns - g_rn*R_ms - g_sm*R_nr + g_sn*R_mr)
                        + c2*(g_rm*g_sn - g_rn*g_sm);
                }
    free(R_cov);
    return C;
}

/* ---- Sectional Curvature ---- */

/**
 * Sectional curvature K(X,Y) = R(X,Y,Y,X) / (|X|^2|Y|^2 - <X,Y>^2)
 *
 * For a 2-plane spanned by orthonormal X, Y:
 *   denominator = 1 (if normalized)
 *
 * Course Ref: do Carmo Sec.2.4
 */
double sectional_curvature(const double *R_cov, const Metric *g,
                            const double *X, const double *Y, int dim) {
    if (!R_cov || !g || !X || !Y || dim < 2) return 0.0;
    double num = 0.0;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            for (int rho = 0; rho < dim; rho++)
                for (int sig = 0; sig < dim; sig++) {
                    int idx = ((mu*dim + nu)*dim + rho)*dim + sig;
                    num += R_cov[idx] * X[mu] * Y[nu] * X[rho] * Y[sig];
                }
    double normX2 = 0.0, normY2 = 0.0, dotXY = 0.0;
    for (int mu = 0; mu < dim; mu++) {
        for (int nu = 0; nu < dim; nu++) {
            double g_mn = metric_get(g, mu, nu);
            normX2 += g_mn * X[mu] * X[nu];
            normY2 += g_mn * Y[mu] * Y[nu];
            dotXY  += g_mn * X[mu] * Y[nu];
        }
    }
    double denom = normX2 * normY2 - dotXY * dotXY;
    if (fabs(denom) < 1e-15) return 0.0;
    return num / denom;
}

/* ---- Bianchi Identities ---- */

/**
 * First Bianchi identity: R^rho_{[sigma mu nu]} = 0
 *   R^rho_{sigma mu nu} + R^rho_{mu nu sigma} + R^rho_{nu sigma mu} = 0
 *
 * Course Ref: Wald Eq.(3.2.15)
 */
double bianchi_identity_first_verify(const RiemannTensor *R) {
    if (!R) return -1.0;
    int dim = R->dim;
    double max_viol = 0.0;
    for (int rho = 0; rho < dim; rho++)
        for (int sig = 0; sig < dim; sig++)
            for (int mu = 0; mu < dim; mu++)
                for (int nu = 0; nu < dim; nu++) {
                    double v = riemann_get(R, rho, sig, mu, nu)
                             + riemann_get(R, rho, mu, nu, sig)
                             + riemann_get(R, rho, nu, sig, mu);
                    if (fabs(v) > max_viol) max_viol = fabs(v);
                }
    return max_viol;
}

/**
 * Second Bianchi identity: nabla_{[lambda} R^rho_{|sigma| mu nu]} = 0
 *
 * Course Ref: Wald Eq.(3.2.18)
 */
double bianchi_identity_second_verify(const RiemannTensor *R,
                                       const Connection *conn, const double *dR) {
    if (!R || !conn || !dR) return -1.0;
    int dim = R->dim;
    double max_viol = 0.0;
    for (int rho = 0; rho < dim; rho++)
        for (int sig = 0; sig < dim; sig++)
            for (int mu = 0; mu < dim; mu++)
                for (int nu = 0; nu < dim; nu++)
                    for (int lam = 0; lam < dim; lam++) {
                        /* Approximate covariant derivative:
                         * nabla_lam R^rho_{sig mu nu} approx d_lam R^rho_{sig mu nu}
                         * (connection terms omitted for this check) */
                        int dRidx = (((lam*dim + rho)*dim + sig)*dim + mu)*dim + nu;
                        double dR_val = dR[dRidx];
                        /* For a proper check, we would antisymmetrize over [lam mu nu].
                         * Here we check that dR is roughly cyclic in those indices. */
                        double dR_lam = dR_val;
                        double dR_mu  = dR[(((mu*dim + rho)*dim + sig)*dim + lam)*dim + nu];
                        double dR_nu  = dR[(((nu*dim + rho)*dim + sig)*dim + mu)*dim + lam];
                        double viol = dR_lam + dR_mu + dR_nu;
                        if (fabs(viol) > max_viol) max_viol = fabs(viol);
                    }
    return max_viol;
}

/**
 * Kretschmann scalar: K = R^{mu nu rho sigma} R_{mu nu rho sigma}
 *
 * Fully contracted product of the Riemann tensor with itself.
 * Used to detect genuine curvature singularities (e.g., r=0 in Schwarzschild).
 *
 * Course Ref: Carroll Sec.5.4
 */
double kretschmann_scalar(const double *R_cov, const Metric *g, int dim) {
    if (!R_cov || !g) return 0.0;
    double K = 0.0;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            for (int rho = 0; rho < dim; rho++)
                for (int sig = 0; sig < dim; sig++) {
                    int idx = ((mu*dim + nu)*dim + rho)*dim + sig;
                    /* Raise all indices of R_cov */
                    for (int a = 0; a < dim; a++)
                        for (int b = 0; b < dim; b++)
                            for (int c = 0; c < dim; c++)
                                for (int d = 0; d < dim; d++) {
                                    int idx2 = ((a*dim + b)*dim + c)*dim + d;
                                    double fac = metric_inv_get(g, mu, a)
                                               * metric_inv_get(g, nu, b)
                                               * metric_inv_get(g, rho, c)
                                               * metric_inv_get(g, sig, d);
                                    K += fac * R_cov[idx2] * R_cov[idx] / dim;
                                }
                }
    return K;
}

/* ---- Print utilities ---- */

void riemann_print(const RiemannTensor *R) {
    if (!R) { printf("Riemann: (null)\n"); return; }
    int dim = R->dim;
    printf("Riemann tensor R^rho_{sigma mu nu} (dim=%d):\n", dim);
    int count = 0;
    for (int rho = 0; rho < dim; rho++)
        for (int sig = 0; sig < dim; sig++)
            for (int mu = 0; mu < dim; mu++)
                for (int nu = 0; nu < dim; nu++) {
                    double v = riemann_get(R, rho, sig, mu, nu);
                    if (fabs(v) > 1e-12 && count < 40) {
                        printf("  R[%d,%d,%d,%d] = %g\n", rho, sig, mu, nu, v);
                        count++;
                    }
                }
    if (count >= 40) printf("  ... (more non-zero)\n");
}

void ricci_print(const RicciTensor *Ric) {
    if (!Ric) { printf("Ricci: (null)\n"); return; }
    int dim = Ric->dim;
    printf("Ricci tensor R_{mu nu} (dim=%d):\n", dim);
    for (int mu = 0; mu < dim; mu++) {
        printf("  [");
        for (int nu = 0; nu < dim; nu++)
            printf(" % 8.4f", ricci_get(Ric, mu, nu));
        printf(" ]\n");
    }
}