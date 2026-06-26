#include "dg_geodesic.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/**
 * RHS of geodesic equation for ODE integration.
 *
 * State y = [x^0, ..., x^{n-1}, v^0, ..., v^{n-1}]
 * dy/dtau = [v^0, ..., v^{n-1}, a^0, ..., a^{n-1}]
 * where a^mu = -Gamma^mu_{alpha beta} v^alpha v^beta
 *
 * Course Ref: Wald Eq.(3.3.5)
 */
void geodesic_rhs(double tau, const double *y, double *dydt, void *params) {
    (void)tau;
    Connection *conn = (Connection*)params;
    int dim = conn->dim;

    /* Positions: dy^mu/dtau = v^mu */
    for (int mu = 0; mu < dim; mu++) {
        dydt[mu] = y[dim + mu];
    }

    /* Velocities: dv^mu/dtau = -Gamma^mu_{alpha beta} v^alpha v^beta */
    for (int mu = 0; mu < dim; mu++) {
        double acc = 0.0;
        for (int alpha = 0; alpha < dim; alpha++) {
            double v_alpha = y[dim + alpha];
            for (int beta = 0; beta < dim; beta++) {
                double v_beta = y[dim + beta];
                acc -= connection_get(conn, mu, alpha, beta)
                     * v_alpha * v_beta;
            }
        }
        dydt[dim + mu] = acc;
    }
}

/**
 * Single RK4 step for geodesic integration.
 *
 * Standard fourth-order Runge-Kutta:
 *   k1 = h * f(t_n, y_n)
 *   k2 = h * f(t_n + h/2, y_n + k1/2)
 *   k3 = h * f(t_n + h/2, y_n + k2/2)
 *   k4 = h * f(t_n + h, y_n + k3)
 *   y_{n+1} = y_n + (k1 + 2*k2 + 2*k3 + k4) / 6
 *
 * Course Ref: Thijssen Sec.7.2
 */
int geodesic_step_rk4(GeodesicPoint *point, const Connection *conn, double dtau) {
    if (!point || !conn) return -1;
    int dim = point->dim;
    int nvars = 2 * dim;
    double y[2 * DG_MAX_DIM] = {0};
    double k1[2 * DG_MAX_DIM] = {0}, k2[2 * DG_MAX_DIM] = {0};
    double k3[2 * DG_MAX_DIM] = {0}, k4[2 * DG_MAX_DIM] = {0};
    double ytmp[2 * DG_MAX_DIM] = {0};

    /* Pack state */
    for (int i = 0; i < dim; i++) {
        y[i]       = point->pos[i];
        y[dim + i] = point->vel[i];
    }

    /* k1 */
    geodesic_rhs(point->tau, y, k1, (void*)conn);
    for (int i = 0; i < nvars; i++) k1[i] *= dtau;

    /* k2 */
    for (int i = 0; i < nvars; i++) ytmp[i] = y[i] + 0.5 * k1[i];
    geodesic_rhs(point->tau + 0.5 * dtau, ytmp, k2, (void*)conn);
    for (int i = 0; i < nvars; i++) k2[i] *= dtau;

    /* k3 */
    for (int i = 0; i < nvars; i++) ytmp[i] = y[i] + 0.5 * k2[i];
    geodesic_rhs(point->tau + 0.5 * dtau, ytmp, k3, (void*)conn);
    for (int i = 0; i < nvars; i++) k3[i] *= dtau;

    /* k4 */
    for (int i = 0; i < nvars; i++) ytmp[i] = y[i] + k3[i];
    geodesic_rhs(point->tau + dtau, ytmp, k4, (void*)conn);
    for (int i = 0; i < nvars; i++) k4[i] *= dtau;

    /* Update */
    for (int i = 0; i < nvars; i++)
        y[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

    /* Unpack */
    for (int i = 0; i < dim; i++) {
        point->pos[i] = y[i];
        point->vel[i] = y[dim + i];
    }
    point->tau += dtau;

    return 0;
}

/**
 * Full geodesic integration using RK4.
 *
 * Integrates from initial point for max_steps steps of size step_size.
 * Records position every record_every steps.
 *
 * Returns array of GeodesicPoint (n_steps long), caller must free via
 * geodesic_points_free().
 */
GeodesicPoint* geodesic_integrate(const GeodesicPoint *initial,
                                   const Connection *conn,
                                   const GeodesicParams *params,
                                   int *n_steps) {
    if (!initial || !conn || !params || !n_steps) return NULL;

    int max_records = params->max_steps / params->record_every + 1;
    GeodesicPoint *traj = (GeodesicPoint*)malloc(
        max_records * sizeof(GeodesicPoint));
    if (!traj) return NULL;

    GeodesicPoint current = *initial;
    traj[0] = current;
    int nrec = 1;

    for (int step = 1; step <= params->max_steps; step++) {
        if (geodesic_step_rk4(&current, conn, params->step_size) != 0) {
            break;
        }
        if (step % params->record_every == 0 && nrec < max_records) {
            traj[nrec++] = current;
        }

        /* Check for NaN/inf */
        int bad = 0;
        for (int i = 0; i < current.dim; i++) {
            if (!isfinite(current.pos[i]) || !isfinite(current.vel[i])) {
                bad = 1;
                break;
            }
        }
        if (bad) break;
    }

    *n_steps = nrec;
    return traj;
}

double geodesic_velocity_norm(const GeodesicPoint *point, const Metric *metric) {
    if (!point || !metric) return 0.0;
    double norm = 0.0;
    int dim = point->dim;
    for (int mu = 0; mu < dim; mu++)
        for (int nu = 0; nu < dim; nu++)
            norm += metric_get(metric, mu, nu)
                  * point->vel[mu] * point->vel[nu];
    return norm;
}

double proper_time_total(const GeodesicPoint *points, int n,
                          const Metric *metric) {
    if (!points || n < 2 || !metric) return 0.0;
    double tau_total = 0.0;
    for (int i = 1; i < n; i++) {
        double dx[DG_MAX_DIM];
        for (int mu = 0; mu < points[i].dim; mu++)
            dx[mu] = points[i].pos[mu] - points[i-1].pos[mu];
        double ds2 = metric_line_element(metric, dx);
        if (ds2 < 0) /* timelike */
            tau_total += sqrt(-ds2);
    }
    return tau_total;
}

/**
 * Parallel transport of vector V along a geodesic trajectory.
 *
 * Solves: dV^mu/dtau + Gamma^mu_{alpha beta} V^alpha v^beta = 0
 * using simple Euler integration along the discrete trajectory.
 *
 * Course Ref: Wald Eq.(3.1.1)
 */
int parallel_transport(const double *V_initial, const GeodesicPoint *points,
                        int n, const Connection *conn, double *V_result) {
    if (!V_initial || !points || n < 1 || !conn || !V_result) return -1;
    int dim = points[0].dim;

    /* Copy initial */
    for (int mu = 0; mu < dim; mu++)
        V_result[mu] = V_initial[mu];

    /* Transport along each segment */
    for (int step = 1; step < n; step++) {
        double dtau = points[step].tau - points[step-1].tau;
        if (dtau <= 0) dtau = 1e-6;

        double V_new[DG_MAX_DIM];
        for (int mu = 0; mu < dim; mu++) {
            double dV = 0.0;
            for (int alpha = 0; alpha < dim; alpha++)
                for (int beta = 0; beta < dim; beta++)
                    dV -= connection_get(conn, mu, alpha, beta)
                        * V_result[alpha] * points[step].vel[beta];
            V_new[mu] = V_result[mu] + dV * dtau;
        }
        for (int mu = 0; mu < dim; mu++)
            V_result[mu] = V_new[mu];
    }
    return 0;
}

/**
 * RHS of geodesic deviation equation:
 *   D^2 xi^mu / dtau^2 = -R^mu_{nu rho sigma} u^nu xi^rho u^sigma
 *
 * Where D/dtau is the absolute derivative.
 *
 * Course Ref: Wald Eq.(3.3.14)
 */
void geodesic_deviation_rhs(const double *xi, const double *dxi_dtau,
                             const GeodesicPoint *point, const RiemannTensor *R,
                             double *ddxi) {
    (void)dxi_dtau;
    if (!xi || !point || !R || !ddxi) return;
    int dim = point->dim;
    for (int mu = 0; mu < dim; mu++) {
        double sum = 0.0;
        for (int nu = 0; nu < dim; nu++)
            for (int rho = 0; rho < dim; rho++)
                for (int sig = 0; sig < dim; sig++)
                    sum -= riemann_get(R, mu, nu, rho, sig)
                         * point->vel[nu] * xi[rho] * point->vel[sig];
        ddxi[mu] = sum;
    }
}

/**
 * Set up initial conditions for a timelike geodesic in Schwarzschild.
 *
 * Schwarzschild metric: ds^2 = -(1-2M/r)dt^2 + dr^2/(1-2M/r)
 *   + r^2(dtheta^2 + sin^2(theta) dphi^2)
 *
 * Initial position: t=0, r=r0, theta=pi/2, phi=0 (equatorial plane)
 * Initial 4-velocity normalized: u^mu u_mu = -1
 *
 * For a particle at rest at infinity with impact parameter b:
 *   E = 1 (energy at infinity), L = b (angular momentum)
 *
 * Course Ref: Wald Sec.6.3, Carroll Sec.5.4
 */
void geodesic_schwarzschild_init(double r0, double b, double M,
                                  GeodesicPoint *point) {
    if (!point) return;
    point->dim = 4;
    point->chart_index = 0;
    point->tau = 0.0;

    /* Position: (t=0, r=r0, theta=pi/2, phi=0) */
    point->pos[0] = 0.0;
    point->pos[1] = r0;
    point->pos[2] = 1.5707963267948966;  /* pi/2 */
    point->pos[3] = 0.0;

    /* Constants of motion for equatorial geodesic:
     * E = -g_{tt} u^t, L = g_{phi phi} u^phi */
    double E = 1.0;
    double L = b;
    double f = 1.0 - 2.0 * M / r0;

    /* u^t = E / f */
    point->vel[0] = E / f;
    /* u^r determined from normalization: g_{mu nu} u^mu u^nu = -1 */
    double ut2 = f * point->vel[0] * point->vel[0]; /* effective potential term for E */
    double uph2 = L * L / (r0 * r0);                 /* effective potential term for L */
    /* -1 = -f (u^t)^2 + (u^r)^2 / f + r^2 (u^phi)^2
     * => (u^r)^2 / f = -1 + f (u^t)^2 - r^2 (u^phi)^2 */
    double ur2 = f * (ut2 - 1.0 - uph2 / (r0 * r0) * r0 * r0 / f);
    /* Actually: -1 = -(1-2M/r)(u^t)^2 + (1/(1-2M/r))(u^r)^2 + r^2(u^phi)^2
     * With u^phi = L/r^2:
     * (u^r)^2 = (u^t)^2 - (1-2M/r)(1 + L^2/r^2)
     */
    ur2 = point->vel[0] * point->vel[0]
        - f * (1.0 + L * L / (r0 * r0));
    point->vel[1] = (ur2 > 0) ? -sqrt(ur2) : 0.0;  /* inward if falling */

    point->vel[2] = 0.0;               /* theta initially constant */
    point->vel[3] = L / (r0 * r0);     /* u^phi = L / r^2 */
}

/**
 * Set up initial conditions for a null geodesic in FRW metric.
 *
 * FRW metric: ds^2 = -dt^2 + a^2(t)[dr^2/(1-kr^2) + r^2 dOmega^2]
 *
 * For a radial light ray (dtheta=dphi=0): ds^2 = 0 => dt = a(t) dr
 *
 * Course Ref: Dodelson Sec.1.2
 */
void geodesic_frw_null_init(double a0, GeodesicPoint *point) {
    if (!point) return;
    point->dim = 4;
    point->chart_index = 0;
    point->tau = 0.0;

    /* Position: t=0, r=0, theta=pi/2, phi=0 */
    point->pos[0] = 0.0;
    point->pos[1] = 0.0;
    point->pos[2] = 1.5707963267948966;
    point->pos[3] = 0.0;

    /* Null condition: g_{mu nu} u^mu u^nu = 0
     * For radial null: -(u^t)^2 + a^2(t) (u^r)^2 = 0
     * => u^r = u^t / a(t)
     * Normalize u^t = 1 (affine parameter scaling) */
    point->vel[0] = 1.0;
    point->vel[1] = 1.0 / a0;
    point->vel[2] = 0.0;
    point->vel[3] = 0.0;
}

/* ---- Utilities ---- */

void geodesic_point_print(const GeodesicPoint *point) {
    if (!point) return;
    printf("tau=%.6f pos=[", point->tau);
    for (int i = 0; i < point->dim; i++)
        printf("%.4f%s", point->pos[i], i < point->dim - 1 ? ", " : "");
    printf("] vel=[");
    for (int i = 0; i < point->dim; i++)
        printf("%.4f%s", point->vel[i], i < point->dim - 1 ? ", " : "");
    printf("]\n");
}

void geodesic_points_free(GeodesicPoint *points) {
    free(points);
}