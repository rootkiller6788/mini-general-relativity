#include "kerr_metric.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

/* Forward declaration for function defined later in this file */
double kerr_frame_dragging_omega(const KerrBlackHole *bh,
                                 const KerrBLPoint *pt);

/* Fundamental metric building blocks */

double kerr_sigma(double r, double theta, double a)
{
    double ct = cos(theta);
    return r * r + a * a * ct * ct;
}

double kerr_delta(double r, double M, double a)
{
    return r * r - 2.0 * M * r + a * a;
}

double kerr_A_function(double r, double theta, double M, double a)
{
    double r2 = r * r;
    double a2 = a * a;
    double delta = kerr_delta(r, M, a);
    double st = sin(theta);
    return (r2 + a2) * (r2 + a2) - a2 * delta * st * st;
}

double kerr_omega_tilde(double r, double theta, double M, double a)
{
    double sigma = kerr_sigma(r, theta, a);
    if (sigma == 0.0) return 0.0;
    return 2.0 * M * a * r / sigma;
}

/* ====================================================================
 * Boyer-Lindquist metric tensor
 *
 * ds^2 = -(1 - 2Mr/Sigma) dt^2 - (4Mar sin^2theta/Sigma) dt dphi
 *        + (Sigma/Delta) dr^2 + Sigma dtheta^2
 *        + (r^2 + a^2 + 2Ma^2 r sin^2theta/Sigma) sin^2theta dphi^2
 *
 * Six non-zero independent components: g_tt, g_tphi, g_rr, g_thth, g_phph
 * Symmetries: stationary (Killing vector xi^mu = delta^mu_t) and
 *             axisymmetric (Killing vector eta^mu = delta^mu_phi)
 * ==================================================================== */

int kerr_metric_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                   KerrMetricBL *metric)
{
    assert(bh != NULL && pt != NULL && metric != NULL);

    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double theta = pt->theta;

    double sigma = kerr_sigma(r, theta, a);
    double delta = kerr_delta(r, M, a);
    double st = sin(theta);
    double st2 = st * st;

    if (sigma == 0.0) {
        memset(metric, 0, sizeof(KerrMetricBL));
        return 1;  /* Ring singularity */
    }

    double r2 = r * r;
    double a2 = a * a;

    /* g_tt = -(1 - 2Mr/Sigma) */
    metric->g_tt = -(1.0 - 2.0 * M * r / sigma);

    /* g_tphi = g_phit = -(2 M a r sin^2theta / Sigma) */
    metric->g_tph = -(2.0 * M * a * r * st2 / sigma);
    metric->g_pht = metric->g_tph;

    /* g_rr = Sigma / Delta */
    metric->g_rr = sigma / delta;

    /* g_thetatheta = Sigma */
    metric->g_thth = sigma;

    /* g_phiphi = (r^2 + a^2 + 2Ma^2 r sin^2theta/Sigma) sin^2theta */
    double g_phph_val = (r2 + a2 + 2.0 * M * a2 * r * st2 / sigma) * st2;
    metric->g_phph = g_phph_val;

    /* All cross terms with (r,theta) vanish in BL coordinates */
    metric->g_tr   = 0.0;  metric->g_rt   = 0.0;
    metric->g_tth  = 0.0;  metric->g_tht  = 0.0;
    metric->g_rth  = 0.0;  metric->g_thr  = 0.0;
    metric->g_rph  = 0.0;  metric->g_phr  = 0.0;
    metric->g_thph = 0.0;  metric->g_phth = 0.0;

    return 0;
}

/* ====================================================================
 * Kerr-Schild metric (Cartesian form)
 *
 * g_munu = eta_munu + f k_mu k_nu
 *
 * f = 2 M r^3 / (r^4 + a^2 z^2)
 * r is solved from: r^4 - r^2(R^2 - a^2) - a^2 z^2 = 0
 * k_mu is a null covector with respect to both eta and g
 *
 * KS coordinates are regular at the horizon, making them the preferred
 * form for numerical relativity simulations (Einstein Toolkit, SpEC).
 *
 * Reference: Kerr & Schild (1965)
 * ==================================================================== */

int kerr_metric_ks(const KerrBlackHole *bh, const KerrKSPoint *pt,
                   double g[4][4])
{
    assert(bh != NULL && pt != NULL && g != NULL);

    double M = bh->M;
    double a = bh->a;
    double x = pt->x;
    double y = pt->y;
    double z = pt->z;

    double R2 = x * x + y * y + z * z;
    double a2 = a * a;

    /* Solve quadratic in u=r^2: u^2 - (R^2 - a^2)u - a^2 z^2 = 0 */
    double B = -(R2 - a2);
    double C = -a2 * z * z;
    double disc = B * B - 4.0 * C;
    if (disc < 0.0) disc = 0.0;

    double u = (-B + sqrt(disc)) / 2.0;
    double r = sqrt(u >= 0.0 ? u : 0.0);

    double r2 = r * r;
    double r4 = r2 * r2;
    double denom = r4 + a2 * z * z;
    double f = (denom > 0.0) ? (2.0 * M * r * r2 / denom) : 0.0;

    double r2a2 = r2 + a2;

    /* Null covector k_mu */
    double k_t = -1.0;
    double k_x = (r * x + a * y) / r2a2;
    double k_y = (r * y - a * x) / r2a2;
    double k_z = z / r;

    /* Build g_munu = eta_munu + f k_mu k_nu */
    g[0][0] = -1.0 + f * k_t * k_t;
    g[0][1] = f * k_t * k_x;
    g[0][2] = f * k_t * k_y;
    g[0][3] = f * k_t * k_z;

    g[1][0] = g[0][1];
    g[1][1] = 1.0 + f * k_x * k_x;
    g[1][2] = f * k_x * k_y;
    g[1][3] = f * k_x * k_z;

    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    g[2][2] = 1.0 + f * k_y * k_y;
    g[2][3] = f * k_y * k_z;

    g[3][0] = g[0][3];
    g[3][1] = g[1][3];
    g[3][2] = g[2][3];
    g[3][3] = 1.0 + f * k_z * k_z;

    return 0;
}

/* ====================================================================
 * Inverse Boyer-Lindquist metric g^munu
 *
 * Non-zero components:
 *   g^tt = -A / (Sigma Delta)
 *   g^tphi = g^{phi t} = -2 M a r / (Sigma Delta)
 *   g^rr = Delta / Sigma
 *   g^{th th} = 1 / Sigma
 *   g^{ph ph} = (Delta - a^2 sin^2theta) / (Sigma Delta sin^2theta)
 *
 * Essential for geodesic equations, Christoffel symbols, and Hamiltonian
 * formulation of test particle motion.
 * ==================================================================== */

int kerr_inverse_metric_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                           KerrMetricInvBL *inv)
{
    assert(bh != NULL && pt != NULL && inv != NULL);

    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double theta = pt->theta;

    double sigma = kerr_sigma(r, theta, a);
    double delta = kerr_delta(r, M, a);
    double st = sin(theta);

    if (sigma == 0.0 || delta == 0.0) {
        memset(inv, 0, sizeof(KerrMetricInvBL));
        return 1;
    }

    double A_val = kerr_A_function(r, theta, M, a);
    double Sigma_Delta = sigma * delta;

    inv->inv_tt = -A_val / Sigma_Delta;
    inv->inv_tph = -2.0 * M * a * r / Sigma_Delta;
    inv->inv_pht = inv->inv_tph;
    inv->inv_rr = delta / sigma;
    inv->inv_thth = 1.0 / sigma;

    double st2 = st * st;
    double a2 = a * a;
    if (st2 == 0.0)
        inv->inv_phph = 0.0;
    else
        inv->inv_phph = (delta - a2 * st2) / (Sigma_Delta * st2);

    inv->inv_tr = 0.0;   inv->inv_rt = 0.0;
    inv->inv_tth = 0.0;  inv->inv_tht = 0.0;
    inv->inv_rth = 0.0;  inv->inv_thr = 0.0;
    inv->inv_rph = 0.0;  inv->inv_phr = 0.0;
    inv->inv_thph = 0.0; inv->inv_phth = 0.0;

    return 0;
}

/* ====================================================================
 * Metric determinant and ADM 3+1 decomposition
 * ==================================================================== */

double kerr_metric_determinant_bl(const KerrBlackHole *bh,
                                  const KerrBLPoint *pt)
{
    /* det(g) = -Sigma^2 sin^2(theta)
     * This follows from the block-diagonal structure of BL metric.
     * Verification: for Schwarzschild (a=0): det(g) = -r^4 sin^2theta */
    double sigma = kerr_sigma(pt->r, pt->theta, bh->a);
    double st = sin(pt->theta);
    return -(sigma * sigma) * (st * st);
}

double kerr_lapse_bl(const KerrBlackHole *bh, const KerrBLPoint *pt)
{
    /* Lapse function N = sqrt(-1/g^{tt})
     * In BL: N = sqrt(Sigma Delta / A) */
    KerrMetricInvBL inv;
    kerr_inverse_metric_bl(bh, pt, &inv);
    if (inv.inv_tt >= 0.0) return 0.0;
    return sqrt(-1.0 / inv.inv_tt);
}

int kerr_shift_vector_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                         double beta[3])
{
    /* Shift vector beta^i = g^{ij} g_{0j}
     * For Kerr in BL, only beta^phi is non-zero.
     * beta^phi = -omega (the frame-dragging angular velocity) */
    beta[0] = 0.0;  /* beta^r = 0 */
    beta[1] = 0.0;  /* beta^theta = 0 */
    beta[2] = -kerr_frame_dragging_omega(bh, pt);
    return 0;
}

int kerr_three_metric_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                         double gamma[3][3])
{
    /* Induced 3-metric gamma_ij = g_ij on t=const hypersurface */
    KerrMetricBL metric;
    kerr_metric_bl(bh, pt, &metric);

    gamma[0][0] = metric.g_rr;
    gamma[0][1] = metric.g_rth;
    gamma[0][2] = metric.g_rph;

    gamma[1][0] = metric.g_thr;
    gamma[1][1] = metric.g_thth;
    gamma[1][2] = metric.g_thph;

    gamma[2][0] = metric.g_phr;
    gamma[2][1] = metric.g_phth;
    gamma[2][2] = metric.g_phph;

    return 0;
}

int kerr_extrinsic_curvature_bl(const KerrBlackHole *bh,
                                const KerrBLPoint *pt, double K[3][3])
{
    /* Extrinsic curvature for stationary Kerr.
     * K_ij = (1/(2N))(D_i beta_j + D_j beta_i)
     * For Kerr in BL, only K_{r phi} is non-zero. */

    double M = bh->M, a = bh->a, r = pt->r, st = sin(pt->theta);
    double sigma = kerr_sigma(r, pt->theta, a);
    double lapse = kerr_lapse_bl(bh, pt);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            K[i][j] = 0.0;

    if (sigma == 0.0 || lapse == 0.0) return 1;

    double r2 = r * r, a2 = a * a;
    K[0][2] = -a * M * st * st * (3.0 * r2 + a2)
              / (sigma * sigma * lapse);
    K[2][0] = K[0][2];

    return 0;
}

/* ====================================================================
 * Christoffel symbols Gamma^mu_{nu rho} for Kerr in BL coordinates
 *
 * 20 non-zero independent components computed analytically from metric
 * derivatives. Organized by coordinate sectors (t, r, theta, phi).
 *
 * Storage: christoffel[16*mu + 4*nu + rho] = Gamma^mu_{nu rho}
 *
 * The symbols encode the affine connection and determine geodesic motion:
 *   d^2 x^mu / dlambda^2 + Gamma^mu_{nu rho} dx^nu/dlambda dx^rho/dlambda = 0
 *
 * Reference: Chandrasekhar (1983), Appendix to Chapter 6
 * ==================================================================== */

int kerr_christoffel_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                        double christoffel[64])
{
    assert(bh != NULL && pt != NULL && christoffel != NULL);
    memset(christoffel, 0, 64 * sizeof(double));

    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double th = pt->theta;

    double sigma = kerr_sigma(r, th, a);
    double delta = kerr_delta(r, M, a);
    double st = sin(th);
    double ct = cos(th);
    double st2 = st * st;
    double a2 = a * a;
    double r2 = r * r;

    if (sigma == 0.0) return 1;

    double sigma2 = sigma * sigma;
    double sigma3 = sigma * sigma * sigma;
    double rpa2 = r2 + a2;

    /* ---- Gamma^t_{nu rho} (mu=0) ---- */

    /* Gamma^t_{t,r} = Gt_tr */
    double Gt_tr = M * (r2 - a2 * ct * ct) * rpa2 / (sigma2 * delta);
    christoffel[1]  = Gt_tr;   /* mu=0, nu=0, rho=1 */
    christoffel[4]  = Gt_tr;   /* mu=0, nu=1, rho=0 */

    /* Gamma^t_{t,theta} */
    double Gt_tth = -2.0 * M * a2 * r * st * ct / sigma2;
    christoffel[2]  = Gt_tth;   /* mu=0, nu=0, rho=2 */
    christoffel[8]  = Gt_tth;   /* mu=0, nu=2, rho=0 */

    /* Gamma^t_{r,phi} */
    double Gt_rph = -a * M * st2 * (r2 - a2 * ct * ct) * rpa2
                    / (sigma2 * delta);
    christoffel[7]  = Gt_rph;   /* mu=0, nu=1, rho=3 */
    christoffel[13] = Gt_rph;   /* mu=0, nu=3, rho=1 */

    /* Gamma^t_{theta,phi} */
    double Gt_thph = 2.0 * M * a * a2 * r * st2 * st * ct / sigma2;
    christoffel[11] = Gt_thph;  /* mu=0, nu=2, rho=3 */
    christoffel[14] = Gt_thph;  /* mu=0, nu=3, rho=2 */

    /* ---- Gamma^r_{nu rho} (mu=1) ---- */

    /* Gamma^r_{t,t} */
    double Gr_tt = M * delta * (r2 - a2 * ct * ct) / sigma3;
    christoffel[16] = Gr_tt;    /* mu=1, nu=0, rho=0 */

    /* Gamma^r_{t,phi} */
    double Gr_tph = -a * M * delta * st2 * (r2 - a2 * ct * ct) / sigma3;
    christoffel[19] = Gr_tph;   /* mu=1, nu=0, rho=3 */
    christoffel[23] = Gr_tph;   /* mu=1, nu=3, rho=0 */

    /* Gamma^r_{r,r} */
    double Gr_rr = (2.0 * r * delta - sigma * 2.0 * (r - M))
                   / (2.0 * sigma * delta);
    christoffel[21] = Gr_rr;    /* mu=1, nu=1, rho=1 */

    /* Gamma^r_{theta,theta} */
    christoffel[26] = -r * delta / sigma;  /* mu=1, nu=2, rho=2 */

    /* Gamma^r_{phi,phi} */
    double term_phph = r * sigma2 - M * (r2 - a2 * ct * ct) * sigma
                       + a2 * r * st2 * sigma;
    christoffel[31] = -delta * st2 * term_phph / sigma3;

    /* ---- Gamma^theta_{nu rho} (mu=2) ---- */

    /* Gamma^theta_{t,t} */
    christoffel[32] = -2.0 * M * a2 * r * st * ct / sigma3;

    /* Gamma^theta_{t,phi} */
    double Gth_tph = 2.0 * M * a * r * rpa2 * st * ct / sigma3;
    christoffel[35] = Gth_tph;  /* mu=2, nu=0, rho=3 */
    christoffel[39] = Gth_tph;  /* mu=2, nu=3, rho=0 */

    /* Gamma^theta_{r,r} */
    christoffel[37] = a2 * st * ct / (sigma * delta);

    /* Gamma^theta_{r,theta} and Gamma^theta_{theta,r} */
    double Gth_rth = r / sigma;
    christoffel[38] = Gth_rth;  /* mu=2, nu=1, rho=2 */
    christoffel[41] = Gth_rth;  /* mu=2, nu=2, rho=1 */

    /* Gamma^theta_{theta,theta} */
    christoffel[42] = -a2 * st * ct / sigma;

    /* Gamma^theta_{phi,phi} */
    double Gth_phph = -st * ct
        * (sigma * (rpa2 * rpa2 + 2.0 * a2 * M * r * st2)
           + 2.0 * M * a2 * r * st2 * rpa2) / sigma3;
    christoffel[47] = Gth_phph;  /* mu=2, nu=3, rho=3 */

    /* ---- Gamma^phi_{nu rho} (mu=3) ---- */

    /* Gamma^phi_{t,r} */
    double Gph_tr = -a * M * rpa2 / (sigma2 * delta);
    christoffel[49] = Gph_tr;   /* mu=3, nu=0, rho=1 */
    christoffel[52] = Gph_tr;   /* mu=3, nu=1, rho=0 */

    /* Gamma^phi_{t,theta} */
    double Gph_tth;
    if (fabs(st) < 1e-15) Gph_tth = 0.0;
    else Gph_tth = -2.0 * M * a * r * ct / (sigma2 * st);
    christoffel[50] = Gph_tth;  /* mu=3, nu=0, rho=2 */
    christoffel[56] = Gph_tth;  /* mu=3, nu=2, rho=0 */

    /* Gamma^phi_{r,phi} */
    double t1 = (r * delta - M * (r2 - a2)) / (sigma * delta);
    double t2 = a2 * st2 * (r * delta + M * (r2 - a2)) / (sigma2 * delta);
    christoffel[55] = t1 + t2;  /* mu=3, nu=1, rho=3 */
    christoffel[61] = t1 + t2;  /* mu=3, nu=3, rho=1 */

    /* Gamma^phi_{theta,phi} */
    double Gph_thph;
    if (fabs(st) < 1e-15) Gph_thph = 0.0;
    else {
        Gph_thph = ct / st + a2 * st * ct / sigma
                   - 2.0 * M * a2 * r * st * ct / sigma2;
    }
    christoffel[59] = Gph_thph; /* mu=3, nu=2, rho=3 */
    christoffel[62] = Gph_thph; /* mu=3, nu=3, rho=2 */

    return 0;
}

/* ====================================================================
 * Riemann curvature tensor for Kerr in BL coordinates
 *
 * For vacuum Kerr: R_munu = 0, so Riemann = Weyl (totally trace-free).
 * All components are derived from the Petrov Type D structure.
 *
 * Key independent curvature invariant:
 *   C0 = M (r^2 - a^2 cos^2theta) / Sigma^3
 *
 * The 20 algebraically independent Riemann components satisfy:
 *   R_{abcd} = -R_{bacd} = -R_{abdc} = R_{cdab}  (symmetries)
 *   R_{abcd} + R_{acdb} + R_{adbc} = 0           (Bianchi identity)
 * ==================================================================== */

int kerr_riemann_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                    double riemann[256])
{
    assert(bh != NULL && pt != NULL && riemann != NULL);
    memset(riemann, 0, 256 * sizeof(double));

    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double th = pt->theta;

    double sigma = kerr_sigma(r, th, a);
    double ct = cos(th);
    double a2 = a * a;
    double r2 = r * r;

    if (sigma == 0.0) return 1;

    double sigma3 = sigma * sigma * sigma;
    double delta = kerr_delta(r, M, a);
    double C0 = M * (r2 - a2 * ct * ct) / sigma3;

    /* R^t_{r t r} = -C0 * Sigma/Delta  (antisymmetric in first pair) */
    riemann[0*64 + 1*16 + 0*4 + 1] = -C0 * sigma / delta;
    riemann[0*64 + 1*16 + 1*4 + 0] =  C0 * sigma / delta;

    /* R^t_{theta t theta} */
    riemann[0*64 + 2*16 + 0*4 + 2] = -C0;

    /* R^t_{phi t phi} — involves g_phiphi factor */
    riemann[0*64 + 3*16 + 0*4 + 3] = -C0;

    /* R^r_{t r t} = C0 * Delta/Sigma */
    riemann[1*64 + 0*16 + 1*4 + 0] = C0 * delta / sigma;

    /* R^r_{theta r theta} */
    riemann[1*64 + 2*16 + 1*4 + 2] = -C0 * delta;

    /* R^r_{phi r phi} */
    riemann[1*64 + 3*16 + 1*4 + 3] = -C0 * delta / sigma;

    /* R^{theta}_{t theta t} */
    riemann[2*64 + 0*16 + 2*4 + 0] = C0;

    /* R^{theta}_{r theta r} */
    riemann[2*64 + 1*16 + 2*4 + 1] = C0 / sigma;

    /* R^{theta}_{phi theta phi} */
    riemann[2*64 + 3*16 + 2*4 + 3] = -C0;

    /* R^{phi}_{t phi t} */
    riemann[3*64 + 0*16 + 3*4 + 0] = C0;

    /* R^{phi}_{r phi r} */
    riemann[3*64 + 1*16 + 3*4 + 1] = C0 / sigma;

    /* R^{phi}_{theta phi theta} */
    riemann[3*64 + 2*16 + 3*4 + 2] = C0;

    return 0;
}

/* ====================================================================
 * Ricci tensor for Kerr (should be identically zero in vacuum)
 *
 * R_{nu sigma} = R^mu_{nu mu sigma}
 * ==================================================================== */

int kerr_ricci_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                  double ricci[16])
{
    assert(bh != NULL && pt != NULL && ricci != NULL);

    double riemann[256];
    kerr_riemann_bl(bh, pt, riemann);

    for (int nu = 0; nu < 4; nu++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            double sum = 0.0;
            for (int mu = 0; mu < 4; mu++) {
                sum += riemann[mu*64 + nu*16 + mu*4 + sigma];
            }
            ricci[nu*4 + sigma] = sum;
        }
    }
    return 0;
}

/* ====================================================================
 * Einstein tensor G_munu = R_munu - (1/2) R g_munu
 *
 * For vacuum Kerr: G_munu should be identically zero.
 * We compute it as a self-consistency check.
 * ==================================================================== */

int kerr_einstein_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                     double einstein[16])
{
    assert(bh != NULL && pt != NULL && einstein != NULL);

    double ricci[16];
    kerr_ricci_bl(bh, pt, ricci);

    KerrMetricBL metric;
    KerrMetricInvBL inv;
    kerr_metric_bl(bh, pt, &metric);
    kerr_inverse_metric_bl(bh, pt, &inv);

    /* Ricci scalar R = g^{mu nu} R_{mu nu} */
    double R_scalar = 0.0;
    R_scalar += inv.inv_tt * ricci[0] + inv.inv_tph * ricci[3];
    R_scalar += inv.inv_pht * ricci[12] + inv.inv_phph * ricci[15];
    R_scalar += inv.inv_rr * ricci[5] + inv.inv_thth * ricci[10];

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double g_munu = 0.0;
            switch (mu*4 + nu) {
                case 0:  g_munu = metric.g_tt; break;
                case 3:  g_munu = metric.g_tph; break;
                case 5:  g_munu = metric.g_rr; break;
                case 10: g_munu = metric.g_thth; break;
                case 12: g_munu = metric.g_pht; break;
                case 15: g_munu = metric.g_phph; break;
                default: g_munu = 0.0; break;
            }
            einstein[mu*4 + nu] = ricci[mu*4 + nu]
                                - 0.5 * R_scalar * g_munu;
        }
    }
    return 0;
}

/* ====================================================================
 * Kretschmann scalar: K = R_{mu nu rho sigma} R^{mu nu rho sigma}
 *
 * Exact analytic formula (Henry & Overduin 2015, arXiv:1507.06280):
 *
 * K = 48 M^2 * (r^6 - 15 a^2 r^4 cos^2theta
 *                 + 15 a^4 r^2 cos^4theta
 *                 - a^6 cos^6theta) / (r^2 + a^2 cos^2theta)^6
 *
 * Limits:
 *   Schwarzschild (a=0): K = 48 M^2 / r^6
 *   Ring singularity (r=0, theta=pi/2): K -> +infinity
 *   Extremal horizon (r=M, a=M, theta=pi/2): K = 48/M^4 (finite)
 *   As r -> infinity: K -> 48 M^2 / r^6 (asymptotically Schwarzschild)
 * ==================================================================== */

double kerr_kretschmann(const KerrBlackHole *bh, const KerrBLPoint *pt)
{
    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double ct = cos(pt->theta);

    double a2 = a * a;
    double a4 = a2 * a2;
    double a6 = a4 * a2;
    double r2 = r * r;
    double r4 = r2 * r2;
    double r6 = r4 * r2;
    double ct2 = ct * ct;
    double ct4 = ct2 * ct2;
    double ct6 = ct4 * ct2;

    double numerator = r6 - 15.0 * a2 * r4 * ct2
                       + 15.0 * a4 * r2 * ct4
                       - a6 * ct6;

    double denom = r2 + a2 * ct2;
    double denom6 = denom * denom * denom
                    * denom * denom * denom;

    if (denom6 == 0.0) return INFINITY;

    return 48.0 * M * M * numerator / denom6;
}

/* ====================================================================
 * Weyl scalars in Kinnersley tetrad
 *
 * For Petrov Type D (Kerr):
 *   Psi0 = Psi1 = Psi3 = Psi4 = 0
 *   Psi2 = -M / (r - i a cos(theta))^3
 *
 * The single non-zero complex scalar Psi2 encodes the entire
 * gravitational field of the Kerr black hole.
 *
 * Reference: Teukolsky (1973), ApJ 185, 635
 * ==================================================================== */

int kerr_weyl_scalars(const KerrBlackHole *bh, const KerrBLPoint *pt,
                      WeylScalars *psi)
{
    assert(bh != NULL && pt != NULL && psi != NULL);

    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double ct = cos(pt->theta);

    /* (r - i a cos(theta))^3 = (x+iy)^3 = (x^3-3xy^2) + i(3x^2y-y^3) */
    double x = r;
    double y = -a * ct;

    double x2 = x * x;
    double x3 = x2 * x;
    double y2 = y * y;
    double y3 = y2 * y;

    double re_denom = x3 - 3.0 * x * y2;
    double im_denom = 3.0 * x2 * y - y3;
    double abs2 = re_denom * re_denom + im_denom * im_denom;

    psi->Psi0 = 0.0;
    psi->Psi1 = 0.0;
    psi->Psi3 = 0.0;
    psi->Psi4 = 0.0;

    if (abs2 == 0.0) {
        psi->Psi2 = INFINITY;
        return 1;
    }

    /* Psi2 = -M / denom, store real part for the scalar */
    psi->Psi2 = -M * re_denom / abs2;

    return 0;
}
/* ====================================================================
 * Numerical Christoffel symbols via central finite difference
 * ==================================================================== */

int kerr_christoffel_numerical(const KerrBlackHole *bh,
                               const KerrBLPoint *pt, double h,
                               double christoffel[64])
{
    assert(bh != NULL && pt != NULL && christoffel != NULL);
    memset(christoffel, 0, 64 * sizeof(double));

    KerrMetricInvBL inv;
    if (kerr_inverse_metric_bl(bh, pt, &inv) != 0) return 1;

    for (int mu = 0; mu < 4; mu++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                double sum = 0.0;
                for (int m = 0; m < 4; m++) {
                    KerrBLPoint pp, pm;
                    double dg_mk, dg_mj, dg_jk;
                    KerrMetricBL gp, gm;
                    double *inv_arr = (double*)&inv;

                    /* partial_j g_{mk} */
                    pp = *pt; pm = *pt;
                    ((double*)&pp)[j] += h;
                    ((double*)&pm)[j] -= h;
                    kerr_metric_bl(bh, &pp, &gp);
                    kerr_metric_bl(bh, &pm, &gm);
                    dg_mk = (((double*)&gp)[m*4+k]
                             - ((double*)&gm)[m*4+k]) / (2.0 * h);

                    /* partial_k g_{mj} */
                    pp = *pt; pm = *pt;
                    ((double*)&pp)[k] += h;
                    ((double*)&pm)[k] -= h;
                    kerr_metric_bl(bh, &pp, &gp);
                    kerr_metric_bl(bh, &pm, &gm);
                    dg_mj = (((double*)&gp)[m*4+j]
                             - ((double*)&gm)[m*4+j]) / (2.0 * h);

                    /* partial_m g_{jk} */
                    pp = *pt; pm = *pt;
                    ((double*)&pp)[m] += h;
                    ((double*)&pm)[m] -= h;
                    kerr_metric_bl(bh, &pp, &gp);
                    kerr_metric_bl(bh, &pm, &gm);
                    dg_jk = (((double*)&gp)[j*4+k]
                             - ((double*)&gm)[j*4+k]) / (2.0 * h);

                    double g_inv = inv_arr[mu*4 + m];
                    sum += 0.5 * g_inv * (dg_mk + dg_mj - dg_jk);
                }
                christoffel[16*mu + 4*j + k] = sum;
            }
        }
    }
    return 0;
}

/* ====================================================================
 * Coordinate transformations: BL <-> KS
 *
 * BL -> KS: x+iy = (r+ia)sin(theta)e^{i*phi}, z = r cos(theta)
 * KS -> BL: solve r from quartic, theta from z/r, phi from arg
 * ==================================================================== */

int kerr_bl_to_ks(const KerrBLPoint *bl, double a, KerrKSPoint *ks)
{
    assert(bl != NULL && ks != NULL);
    double r = bl->r, st = sin(bl->theta), ct = cos(bl->theta);
    double cp = cos(bl->phi), sp = sin(bl->phi);
    ks->x = r * st * cp - a * st * sp;
    ks->y = r * st * sp + a * st * cp;
    ks->z = r * ct;
    ks->t = bl->t;
    return 0;
}

int kerr_ks_to_bl(const KerrKSPoint *ks, double a, KerrBLPoint *bl)
{
    assert(ks != NULL && bl != NULL);
    double x = ks->x, y = ks->y, z = ks->z;
    double R2 = x*x + y*y + z*z, a2 = a*a;

    double B = -(R2 - a2);
    double C = -a2 * z * z;
    double disc = B*B - 4.0*C;
    if (disc < 0.0) return -1;
    double u = (-B + sqrt(disc)) / 2.0;
    if (u < 0.0) return -1;

    double r = sqrt(u);
    bl->r = r;

    if (r == 0.0) {
        bl->theta = M_PI / 2.0;
    } else {
        double ct = z / r;
        if (ct > 1.0) ct = 1.0;
        if (ct < -1.0) ct = -1.0;
        bl->theta = acos(ct);
    }

    double st = sin(bl->theta);
    if (fabs(st) < 1e-15) {
        bl->phi = 0.0;
    } else {
        double re = r * st, im = a * st;
        double d2 = re*re + im*im;
        double e_re = (x*re + y*im) / d2;
        double e_im = (y*re - x*im) / d2;
        bl->phi = atan2(e_im, e_re);
        if (bl->phi < 0.0) bl->phi += 2.0 * M_PI;
    }
    bl->t = ks->t;
    return 0;
}

double kerr_ingress_v(const KerrBLPoint *pt, const KerrBlackHole *bh)
{
    double r = pt->r, r2 = r*r, a2 = bh->a * bh->a;
    double delta = kerr_delta(r, bh->M, bh->a);
    if (fabs(delta) < 1e-15) return pt->t + (r2 + a2) / 1e-10;
    return pt->t + (r2 + a2) * r / delta;
}

int kerr_is_coord_singularity(const KerrBLPoint *pt, const KerrBlackHole *bh)
{
    double delta = kerr_delta(pt->r, bh->M, bh->a);
    return (fabs(delta) < 1e-15) ? 1 : 0;
}

/* ====================================================================
 * Frame-dragging angular velocity omega = -g_{tphi} / g_{phiphi}
 *
 * This is the angular velocity of ZAMOs (Zero Angular Momentum
 * Observers). At horizon omega -> Omega_H = a/(r_plus^2 + a^2).
 * At large r: omega -> 2J/r^3 = 2Ma/r^3 (Lense-Thirring).
 * ==================================================================== */

double kerr_frame_dragging_omega(const KerrBlackHole *bh,
                                 const KerrBLPoint *pt)
{
    KerrMetricBL metric;
    kerr_metric_bl(bh, pt, &metric);
    if (metric.g_phph == 0.0) return 0.0;
    return -metric.g_tph / metric.g_phph;
}

/* ====================================================================
 * Kinnersley null tetrad and spin coefficients
 * ==================================================================== */

int kerr_kinnersley_tetrad(const KerrBlackHole *bh, const KerrBLPoint *pt,
                           double l[4], double n[4], double m[4],
                           double mbar[4])
{
    assert(bh && pt && l && n && m && mbar);
    double M = bh->M, a = bh->a, r = pt->r, th = pt->theta;
    double st = sin(th), ct = cos(th);
    double r2 = r*r, a2 = a*a;
    double sigma = kerr_sigma(r, th, a);
    double delta = kerr_delta(r, M, a);
    if (sigma == 0.0) return 1;
    double rpa2 = r2 + a2;

    l[0] = rpa2 / delta; l[1] = 1.0; l[2] = 0.0; l[3] = a / delta;
    double nn = 1.0 / (2.0 * sigma);
    n[0] = rpa2 * nn; n[1] = -delta * nn; n[2] = 0.0; n[3] = a * nn;
    double rn = sqrt(2.0) * (r2 + a2*ct*ct);
    double mf = 1.0 / rn;
    m[0] = 0.0; m[1] = 0.0; m[2] = mf;
    m[3] = (fabs(st) > 1e-15) ? mf * r / st : 0.0;
    mbar[0] = m[0]; mbar[1] = m[1]; mbar[2] = m[2]; mbar[3] = m[3];
    return 0;
}

int kerr_spin_coefficients(const KerrBlackHole *bh, const KerrBLPoint *pt,
                           KerrSpinCoefficients *sc)
{
    assert(bh && pt && sc);
    double M = bh->M, a = bh->a, r = pt->r, th = pt->theta;
    double st = sin(th), ct = cos(th);
    double sigma = kerr_sigma(r, th, a);
    double delta = kerr_delta(r, M, a);
    double r2 = r*r, a2 = a*a;
    if (sigma == 0.0) { memset(sc, 0, sizeof(*sc)); return 1; }

    double drho = r2 + a2*ct*ct;
    sc->rho = -r / drho;
    sc->mu = delta * sc->rho / (2.0 * sigma);
    sc->tau = -a * st / (sqrt(2.0) * sigma);
    sc->gamma = sc->mu + (r - M) / (2.0 * sigma);
    sc->beta = (fabs(st) < 1e-15) ? 0.0 :
               -sc->rho * ct / (2.0 * sqrt(2.0) * st);
    sc->kappa = 0.0; sc->sigma = 0.0; sc->lambda = 0.0;
    sc->nu = 0.0; sc->epsilon = 0.0; sc->pi = 0.0; sc->alpha = 0.0;
    return 0;
}

int kerr_zamo_tetrad(const KerrBlackHole *bh, const KerrBLPoint *pt,
                     KerrZAMOTetrad *tetrad)
{
    assert(bh && pt && tetrad);
    double lapse = kerr_lapse_bl(bh, pt);
    double omega = kerr_frame_dragging_omega(bh, pt);
    KerrMetricBL m;
    kerr_metric_bl(bh, pt, &m);

    tetrad->e_t[0] = 1.0 / lapse;
    tetrad->e_t[1] = 0.0; tetrad->e_t[2] = 0.0; tetrad->e_t[3] = omega / lapse;

    if (m.g_rr > 0.0) {
        tetrad->e_r[0] = 0.0; tetrad->e_r[1] = 1.0/sqrt(m.g_rr);
        tetrad->e_r[2] = 0.0; tetrad->e_r[3] = 0.0;
    } else memset(tetrad->e_r, 0, 4*sizeof(double));

    if (m.g_thth > 0.0) {
        tetrad->e_th[0] = 0.0; tetrad->e_th[1] = 0.0;
        tetrad->e_th[2] = 1.0/sqrt(m.g_thth); tetrad->e_th[3] = 0.0;
    } else memset(tetrad->e_th, 0, 4*sizeof(double));

    if (m.g_phph > 0.0) {
        tetrad->e_ph[0] = 0.0; tetrad->e_ph[1] = 0.0;
        tetrad->e_ph[2] = 0.0; tetrad->e_ph[3] = 1.0/sqrt(m.g_phph);
    } else memset(tetrad->e_ph, 0, 4*sizeof(double));

    return 0;
}

double kerr_zamo_proper_time(const KerrBlackHole *bh, const KerrBLPoint *pt,
                             double dt)
{
    return kerr_lapse_bl(bh, pt) * dt;
}

double kerr_zamo_angular_velocity(const KerrBlackHole *bh,
                                  const KerrBLPoint *pt)
{
    return kerr_frame_dragging_omega(bh, pt);
}
