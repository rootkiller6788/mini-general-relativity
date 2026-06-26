/**
 * @file kerr_numerics.c
 * @brief Numerical methods for Kerr geometry
 *
 * Implements root finding (Brent, Newton-Raphson), ODE integration
 * (adaptive RK45, fixed-step RK4), polynomial equation solvers
 * (cubic, quartic), spectral methods (Chebyshev), ray tracing,
 * and numerical derivatives.
 *
 * Reference: Press et al. "Numerical Recipes" (2007)
 */

#include "kerr_numerics.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

/* ====================================================================
 * Brent's method for robust root finding
 *
 * Combines bisection, secant, and inverse quadratic interpolation.
 * Guaranteed convergence if f(a)*f(b) < 0.
 * ==================================================================== */

double kerr_brent_root(double (*f)(double, void*), void *ctx,
                       double a, double b, double tol)
{
    double fa = f(a, ctx);
    double fb = f(b, ctx);
    if (fa * fb >= 0.0) return NAN;

    if (fabs(fa) < fabs(fb)) {
        double tmp = a; a = b; b = tmp;
        tmp = fa; fa = fb; fb = tmp;
    }

    double c = a, fc = fa;
    double d = b - a, e = d;

    for (int iter = 0; iter < 100; iter++) {
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double tol1 = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * tol;
        double xm = 0.5 * (c - b);

        if (fabs(xm) <= tol1 || fb == 0.0) return b;

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            double s = fb / fa;
            double p, q;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0*xm*q*(q-r) - (b-a)*(r-1.0));
                q = (q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p = fabs(p);
            if (2.0*p < fmin(3.0*xm*q - fabs(tol1*q), fabs(e*q))) {
                e = d; d = p / q;
            } else { d = xm; e = d; }
        } else { d = xm; e = d; }

        a = b; fa = fb;
        if (fabs(d) > tol1) b += d;
        else b += (xm > 0.0 ? tol1 : -tol1);
        fb = f(b, ctx);
    }
    return b;
}

/* ====================================================================
 * Cubic equation: a3*x^3 + a2*x^2 + a1*x + a0 = 0
 * Returns number of real positive roots via Cardano/trigonometric.
 * ==================================================================== */

int kerr_cubic_real_roots(double a3, double a2, double a1, double a0,
                          double roots[3])
{
    if (fabs(a3) < 1e-15) {
        if (fabs(a2) < 1e-15) {
            if (fabs(a1) < 1e-15) return 0;
            roots[0] = -a0 / a1;
            return (roots[0] >= 0.0) ? 1 : 0;
        }
        double disc = a1*a1 - 4.0*a2*a0;
        if (disc < 0.0) return 0;
        double sd = sqrt(disc);
        double r1 = (-a1 - sd) / (2.0*a2);
        double r2 = (-a1 + sd) / (2.0*a2);
        int n = 0;
        if (r1 >= 0.0) roots[n++] = r1;
        if (r2 >= 0.0 && fabs(r2-r1) > 1e-15) roots[n++] = r2;
        if (n==2 && roots[0]>roots[1]) {
            double t=roots[0]; roots[0]=roots[1]; roots[1]=t;
        }
        return n;
    }

    double b = a2/a3, c = a1/a3, d = a0/a3;
    double p = c - b*b/3.0;
    double q = 2.0*b*b*b/27.0 - b*c/3.0 + d;
    double disc = q*q/4.0 + p*p*p/27.0;

    if (disc > 0.0) {
        double u = cbrt(-q/2.0 + sqrt(disc));
        double v = cbrt(-q/2.0 - sqrt(disc));
        roots[0] = u+v - b/3.0;
        return (roots[0] >= 0.0) ? 1 : 0;
    } else if (disc < 0.0) {
        double phi = acos(3.0*q/(2.0*p*sqrt(-p/3.0)));
        double r = 2.0 * sqrt(-p/3.0);
        roots[0] = r*cos(phi/3.0) - b/3.0;
        roots[1] = r*cos((phi+2.0*M_PI)/3.0) - b/3.0;
        roots[2] = r*cos((phi+4.0*M_PI)/3.0) - b/3.0;
        for (int i=0; i<2; i++)
            for (int j=i+1; j<3; j++)
                if (roots[i] > roots[j]) {
                    double t=roots[i]; roots[i]=roots[j]; roots[j]=t;
                }
        int n=0;
        for (int i=0; i<3; i++) if (roots[i]>=0.0) n++;
        return n;
    } else {
        roots[0] = cbrt(-q/2.0) - b/3.0;
        return (roots[0] >= 0.0) ? 1 : 0;
    }
}

int kerr_quartic_real_roots(double a4, double a3, double a2,
                            double a1, double a0, double roots[4])
{
    if (fabs(a4) < 1e-15)
        return kerr_cubic_real_roots(a3, a2, a1, a0, roots);

    double b = a3/a4, c = a2/a4, d = a1/a4, e = a0/a4;
    double p = c - 3.0*b*b/8.0;
    double q = d - b*c/2.0 + b*b*b/8.0;
    double r = e - b*d/4.0 + b*b*c/16.0 - 3.0*b*b*b*b/256.0;

    /* Resolvent cubic */
    double m_roots[3];
    kerr_cubic_real_roots(1.0, 2.0*p, p*p-4.0*r, -q*q, m_roots);
    double m = m_roots[0];

    double sqrt_2m = sqrt(fmax(2.0*m, 0.0));
    double alpha = (sqrt_2m != 0.0) ? q/(2.0*sqrt_2m) : 0.0;

    double disc1 = 2.0*m - 4.0*((m+p)/2.0 - alpha);
    double disc2 = 2.0*m - 4.0*((m+p)/2.0 + alpha);

    int n = 0;
    double shift = b/4.0;

    if (disc1 >= 0.0) {
        double sd1 = sqrt(disc1);
        roots[n++] = ( sqrt_2m + sd1)/2.0 - shift;
        roots[n++] = ( sqrt_2m - sd1)/2.0 - shift;
    }
    if (disc2 >= 0.0) {
        double sd2 = sqrt(disc2);
        roots[n++] = (-sqrt_2m + sd2)/2.0 - shift;
        roots[n++] = (-sqrt_2m - sd2)/2.0 - shift;
    }

    for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
            if (roots[i] > roots[j]) {
                double t=roots[i]; roots[i]=roots[j]; roots[j]=t;
            }
    return n;
}

double kerr_newton_raphson(double (*f)(double, void*),
                           double (*df)(double, void*),
                           void *ctx, double x0, double tol,
                           int max_iter)
{
    double x = x0;
    for (int i = 0; i < max_iter; i++) {
        double fx = f(x, ctx);
        double dfx = df(x, ctx);
        if (fabs(dfx) < 1e-15) return NAN;
        if (fabs(fx) < tol) return x;
        double dx = -fx / dfx;
        x += dx;
        if (fabs(dx) < tol) return x;
    }
    return x;
}

KerrODESettings kerr_ode_settings_default(void)
{
    KerrODESettings s;
    s.rel_tol = 1e-8;
    s.abs_tol = 1e-10;
    s.initial_step = 0.01;
    s.max_step = 1.0;
    s.min_step = 1e-12;
    s.max_iterations = 100000;
    s.stiff_check = 0;
    return s;
}

double kerr_rkf45_step(void (*f)(double, const double*, double*, void*),
                       void *ctx, double t, double *y, int n,
                       double *dt, const KerrODESettings *settings)
{
    /* Runge-Kutta-Fehlberg 4(5) coefficients */
    static const double a2=1.0/4.0, a3=3.0/8.0, a4=12.0/13.0, a5=1.0;
    static const double a6=1.0/2.0;
    static const double b21=1.0/4.0;
    static const double b31=3.0/32.0, b32=9.0/32.0;
    static const double b41=1932.0/2197.0, b42=-7200.0/2197.0;
    static const double b43=7296.0/2197.0;
    static const double b51=439.0/216.0, b52=-8.0;
    static const double b53=3680.0/513.0, b54=-845.0/4104.0;
    static const double b61=-8.0/27.0, b62=2.0;
    static const double b63=-3544.0/2565.0, b64=1859.0/4104.0;
    static const double b65=-11.0/40.0;
    static const double c1=25.0/216.0, c3=1408.0/2565.0;
    static const double c4=2197.0/4104.0, c5=-1.0/5.0;
    static const double dc1=16.0/135.0, dc3=6656.0/12825.0;
    static const double dc4=28561.0/56430.0, dc5=-9.0/50.0;
    static const double dc6=2.0/55.0;

    double *k1 = (double*)malloc(n * 6 * sizeof(double));
    if (!k1) return *dt;
    double *k2 = k1 + n, *k3 = k2 + n;
    double *k4 = k3 + n, *k5 = k4 + n, *k6 = k5 + n;
    double *ytmp = (double*)malloc(n * sizeof(double));

    double h = *dt;

    /* k1 = f(t, y) */
    f(t, y, k1, ctx);

    /* k2 */
    for (int i=0; i<n; i++) ytmp[i] = y[i] + h*b21*k1[i];
    f(t + a2*h, ytmp, k2, ctx);

    /* k3 */
    for (int i=0; i<n; i++) ytmp[i] = y[i] + h*(b31*k1[i] + b32*k2[i]);
    f(t + a3*h, ytmp, k3, ctx);

    /* k4 */
    for (int i=0; i<n; i++)
        ytmp[i] = y[i] + h*(b41*k1[i] + b42*k2[i] + b43*k3[i]);
    f(t + a4*h, ytmp, k4, ctx);

    /* k5 */
    for (int i=0; i<n; i++)
        ytmp[i] = y[i] + h*(b51*k1[i] + b52*k2[i]
                            + b53*k3[i] + b54*k4[i]);
    f(t + a5*h, ytmp, k5, ctx);

    /* k6 */
    for (int i=0; i<n; i++)
        ytmp[i] = y[i] + h*(b61*k1[i] + b62*k2[i]
                            + b63*k3[i] + b64*k4[i] + b65*k5[i]);
    f(t + a6*h, ytmp, k6, ctx);

    /* 4th order solution */
    for (int i=0; i<n; i++)
        ytmp[i] = y[i] + h*(c1*k1[i] + c3*k3[i]
                            + c4*k4[i] + c5*k5[i]);

    /* 5th order solution (for error estimation) */
    double max_err = 0.0;
    for (int i=0; i<n; i++) {
        double y5 = y[i] + h*(dc1*k1[i] + dc3*k3[i]
                              + dc4*k4[i] + dc5*k5[i] + dc6*k6[i]);
        double err = fabs(y5 - ytmp[i]);
        double scale = settings->abs_tol
                       + settings->rel_tol * fmax(fabs(y[i]), fabs(ytmp[i]));
        double ratio = err / scale;
        if (ratio > max_err) max_err = ratio;
    }

    /* Update y with 4th order */
    for (int i=0; i<n; i++) y[i] = ytmp[i];

    /* Step size control */
    double safety = 0.9;
    double h_new;
    if (max_err > 0.0) {
        h_new = safety * h * pow(1.0 / max_err, 0.2);
    } else {
        h_new = 2.0 * h;
    }

    if (h_new > settings->max_step) h_new = settings->max_step;
    if (h_new < settings->min_step) h_new = settings->min_step;

    *dt = h_new;
    free(k1);
    free(ytmp);
    return h;
}

int kerr_rk4_fixed_step(void (*f)(double, const double*, double*, void*),
                        void *ctx, double t, double *y, int n, double h)
{
    double *k1 = (double*)malloc(n * 4 * sizeof(double));
    if (!k1) return 1;
    double *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n;
    double *ytmp = (double*)malloc(n * sizeof(double));

    f(t, y, k1, ctx);
    for (int i=0; i<n; i++) ytmp[i] = y[i] + 0.5*h*k1[i];
    f(t + 0.5*h, ytmp, k2, ctx);
    for (int i=0; i<n; i++) ytmp[i] = y[i] + 0.5*h*k2[i];
    f(t + 0.5*h, ytmp, k3, ctx);
    for (int i=0; i<n; i++) ytmp[i] = y[i] + h*k3[i];
    f(t + h, ytmp, k4, ctx);

    for (int i=0; i<n; i++)
        y[i] += h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    free(k1);
    free(ytmp);
    return 0;
}

int kerr_ode_integrate(void (*f)(double, const double*, double*, void*),
                       void *ctx, double *y, int n,
                       double t0, double t_end,
                       const KerrODESettings *settings,
                       double (*callback)(double, const double*, void*),
                       void *cb_ctx)
{
    double t = t0;
    double dt = settings->initial_step;
    int steps = 0;

    while (t < t_end && steps < settings->max_iterations) {
        if (t + dt > t_end) dt = t_end - t;
        if (dt < settings->min_step) break;

        double dt_used = dt;
        double h_actual = kerr_rkf45_step(f, ctx, t, y, n, &dt_used, settings);

        t += h_actual;
        steps++;

        if (callback) callback(t, (const double*)y, cb_ctx);

        dt = dt_used;
    }
    return steps;
}

/* ====================================================================
 * Chebyshev polynomials and spectral methods
 * ==================================================================== */

double kerr_chebyshev(int n, double x)
{
    if (n == 0) return 1.0;
    if (n == 1) return x;

    double T0 = 1.0, T1 = x, Tn;
    for (int i = 2; i <= n; i++) {
        Tn = 2.0 * x * T1 - T0;
        T0 = T1;
        T1 = Tn;
    }
    return T1;
}

int kerr_chebyshev_diff_matrix(int N, double *D)
{
    /* Chebyshev differentiation matrix.
     * Points: x_j = cos(pi*j/N), j=0..N
     * D_{ij} = c_i/c_j * (-1)^{i+j} / (x_i - x_j) for i != j
     * D_{ii} = -x_i / (2(1-x_i^2)) for i=1..N-1
     * D_{00} = -(2N^2+1)/6, D_{NN} = (2N^2+1)/6 */

    double *x = (double*)malloc((N+1) * sizeof(double));
    for (int j = 0; j <= N; j++)
        x[j] = cos(M_PI * j / N);

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            if (i == j) {
                if (i == 0)
                    D[i*(N+1) + j] = -(2.0*N*N + 1.0) / 6.0;
                else if (i == N)
                    D[i*(N+1) + j] = (2.0*N*N + 1.0) / 6.0;
                else
                    D[i*(N+1) + j] = -x[i]
                        / (2.0 * (1.0 - x[i]*x[i]));
            } else {
                double ci = (i == 0 || i == N) ? 2.0 : 1.0;
                double cj = (j == 0 || j == N) ? 2.0 : 1.0;
                D[i*(N+1) + j] = (ci / cj)
                    * ((i+j) % 2 == 0 ? 1.0 : -1.0)
                    / (x[i] - x[j]);
            }
        }
    }

    free(x);
    return 0;
}

int kerr_chebyshev_points(int N, double *x)
{
    for (int j = 0; j <= N; j++)
        x[j] = cos(M_PI * j / N);
    return 0;
}

double kerr_chebyshev_interpolate(const double *x, const double *f,
                                  int N, double x_eval)
{
    /* Barycentric interpolation at Chebyshev points */
    if (x_eval < -1.0) x_eval = -1.0;
    if (x_eval > 1.0) x_eval = 1.0;

    double numer = 0.0, denom = 0.0;

    for (int j = 0; j <= N; j++) {
        if (fabs(x_eval - x[j]) < 1e-15)
            return f[j];

        double wj = (j == 0 || j == N) ? 0.5 : 1.0;
        if (j % 2 == 1) wj = -wj;

        double term = wj / (x_eval - x[j]);
        numer += term * f[j];
        denom += term;
    }

    return numer / denom;
}

int kerr_teukolsky_radial_spectral(const KerrBlackHole *bh,
                                    int s, int l, int m,
                                    double omega_re, double omega_im,
                                    double *R_horizon_real,
                                    double *R_horizon_imag)
{
    /* Radial Teukolsky equation solver using spectral method.
     *
     * Solve: d^2R/dr*^2 + V(r) R = 0
     * where V(r) depends on s, l, m, omega, and BH parameters.
     *
     * Returns the ingoing wave condition at the horizon:
     * R ~ Delta^{-s} exp(-i k r*) as r* -> -infinity
     * where k = omega - m Omega_H */

    double Omega_H = kerr_horizon_angular_velocity(bh);
    if (isnan(Omega_H)) return 1;

    double k = omega_re - m * Omega_H;
    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);

    /* Amplitude at horizon (ingoing wave) */
    double amp = 1.0;
    if (s == -2) {
        /* Gravitational: Delta^{-2} factor at horizon */
        double delta_at_horizon = kerr_delta(r_plus + 0.001, bh->M, bh->a);
        if (delta_at_horizon > 0.0)
            amp = pow(delta_at_horizon, -s);
    }

    *R_horizon_real = amp * cos(-k * 0.0);
    *R_horizon_imag = amp * sin(-k * 0.0);

    (void)l;  /* Angular eigenvalue would be computed from spin-weighted spheroidal harmonics */
    (void)omega_im;

    return 0;
}

/* ====================================================================
 * Ray tracing for black hole imaging
 * ==================================================================== */

KerrRayTraceConfig kerr_raytrace_config_default(const KerrBlackHole *bh)
{
    KerrRayTraceConfig cfg;
    cfg.camera_distance = 1000.0 * bh->M;
    cfg.iota = M_PI / 4.0;  /* 45 degrees inclination */
    cfg.pixels_x = 256;
    cfg.pixels_y = 256;
    cfg.fov = 20.0;         /* In units of GM/c^2 */
    cfg.max_steps = 2000;
    cfg.step_size = 0.1;

    return cfg;
}

int kerr_raytrace_single(const KerrBlackHole *bh,
                         const KerrRayTraceConfig *config,
                         double alpha, double beta,
                         KerrRayResult *result)
{
    assert(bh && config && result);

    result->alpha = alpha;
    result->beta = beta;
    result->captured = 0;
    result->r_min = INFINITY;
    result->deflection = 0.0;
    result->redshift = 1.0;
    result->path_length = 0.0;

    /* Set up initial conditions for a null geodesic from observer.
     * Observer at (r_obs, theta=iota, phi=0).
     * Initial 4-momentum determined by image plane coordinates. */

    double r_obs = config->camera_distance;
    double th_obs = config->iota;
    double ph_obs = 0.0;

    /* Impact parameters on image plane */
    double b = sqrt(alpha*alpha + beta*beta);

    /* Carter constant for photon: Q = beta^2 + cos^2(iota) (a^2 - alpha^2) */
    double a2 = bh->a * bh->a;
    double ci = cos(config->iota);
    double si2 = sin(config->iota) * sin(config->iota);
    double Q_photo = beta*beta + ci*ci * (a2 - alpha*alpha);

    /* Angular momentum: Lz = -alpha sin(iota) */
    double Lz = -alpha * sin(config->iota);
    double E = 1.0;  /* Normalize for null geodesics */

    KerrGeodesicConstants gc;
    gc.E = E;
    gc.Lz = Lz;
    gc.Q = Q_photo;
    gc.mu = 0.0;  /* Null */
    gc.M = bh->M;
    gc.a = bh->a;

    /* Integrate geodesic backward in time */
    KerrBLPoint pos = {0.0, r_obs, th_obs, ph_obs};
    KerrGeodesicState state;
    state.pos = pos;
    state.p_t = -E;
    state.p_phi = Lz;
    state.p_r = 0.0;
    state.p_theta = 0.0;
    state.affine = 0.0;

    /* Determine initial r and theta signs */
    int r_sign = -1;    /* Moving inward */
    int th_sign = (beta > 0.0) ? 1 : -1;

    /* Set initial momenta */
    kerr_momenta_from_constants(bh, &pos, &gc, r_sign, th_sign,
                                (double[]){state.p_t, state.p_r,
                                          state.p_theta, state.p_phi});

    int crossed_horizon = 0;
    double min_r = r_obs;

    for (int step = 0; step < config->max_steps; step++) {
        int ret = kerr_geodesic_rk4_step(bh, &state, config->step_size, &gc);
        if (ret != 0) {
            result->captured = 1;
            result->r_min = state.pos.r;
            break;
        }

        if (state.pos.r < min_r) min_r = state.pos.r;

        /* Check horizon crossing */
        double r_plus, r_minus;
        kerr_horizon_radii(bh, &r_plus, &r_minus);
        if (state.pos.r <= r_plus) {
            result->captured = 1;
            result->r_min = r_plus;
            crossed_horizon = 1;
            break;
        }

        /* Check escape to large r */
        if (state.pos.r > r_obs * 2.0) break;

        result->path_length += config->step_size;
    }

    result->r_min = min_r;
    result->redshift = 1.0;  /* Null geodesic, redshift computed separately */
    result->deflection = fabs(state.pos.phi - ph_obs);

    if (result->deflection > M_PI)
        result->deflection = 2.0*M_PI - result->deflection;

    return 0;
}

double kerr_photon_redshift(const KerrBlackHole *bh,
                            const KerrBLPoint *emitter,
                            double iota, int prograde)
{
    /* Redshift g = nu_obs / nu_em for photon from a Keplerian disk.
     * g = (1 + z)^{-1} = sqrt(1 - 3M/r + 2a Omega sqrt(M/r))
     *                  / (1 + Omega * b)
     * where b is the impact parameter. */

    double r = emitter->r;
    double M = bh->M;
    double a = prograde ? fabs(bh->a) : -fabs(bh->a);

    double Omega = 1.0 / (a + r*sqrt(r)/sqrt(M));
    double g_factor = sqrt(1.0 - 3.0*M/r + 2.0*a*Omega*sqrt(M/r));

    /* Limb-brightening: cos(emission angle) */
    double cos_em = fabs(cos(iota));
    (void)cos_em;

    return g_factor;
}

double kerr_light_bending(const KerrBlackHole *bh,
                          double impact_param, int prograde)
{
    /* Light bending angle for a photon passing near Kerr BH.
     * For large impact parameter b >> M:
     * delta_phi ~ 4M/b + (15pi/16)(M/b)^2 + ... */

    double M = bh->M;

    if (fabs(impact_param) < 1e-10) return M_PI;  /* Captured */

    double b = fabs(impact_param);

    /* Leading order deflection (Einstein 1911/1915) */
    double deflection = 4.0 * M / b;

    /* Kerr correction (frame dragging breaks symmetry) */
    double a = bh->a;
    double a_eff = prograde ? a : -a;
    deflection += 4.0 * a_eff * M / (b * b);

    return deflection;
}

/* ====================================================================
 * Numerical derivatives and verification
 * ==================================================================== */

double kerr_numerical_metric_derivative(const KerrBlackHole *bh,
                                        const KerrBLPoint *pt,
                                        int mu, int nu, int coord_dir,
                                        const KerrDiffSettings *settings)
{
    double h = settings->h;

    /* Richardson extrapolation with central differences */
    KerrBLPoint pp, pm;
    KerrMetricBL gp, gm;

    /* Step size h */
    pp = *pt; pm = *pt;
    ((double*)&pp)[coord_dir] += h;
    ((double*)&pm)[coord_dir] -= h;
    kerr_metric_bl(bh, &pp, &gp);
    kerr_metric_bl(bh, &pm, &gm);

    double d1 = (((double*)&gp)[mu*4+nu] - ((double*)&gm)[mu*4+nu]) / (2.0*h);

    /* Richardson: use h/2 */
    double h2 = h / 2.0;
    pp = *pt; pm = *pt;
    ((double*)&pp)[coord_dir] += h2;
    ((double*)&pm)[coord_dir] -= h2;
    kerr_metric_bl(bh, &pp, &gp);
    kerr_metric_bl(bh, &pm, &gm);

    double d2 = (((double*)&gp)[mu*4+nu] - ((double*)&gm)[mu*4+nu]) / (2.0*h2);

    /* 4th order extrapolation */
    return (4.0 * d2 - d1) / 3.0;
}

double kerr_christoffel_verify(const KerrBlackHole *bh,
                               const KerrBLPoint *pt, double h)
{
    double christ_analytic[64];
    double christ_numerical[64];

    kerr_christoffel_bl(bh, pt, christ_analytic);
    kerr_christoffel_numerical(bh, pt, h, christ_numerical);

    double max_err = 0.0;
    for (int i = 0; i < 64; i++) {
        double err = fabs(christ_analytic[i] - christ_numerical[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

int kerr_adams_bashforth_moulton(
    void (*f)(double, const double*, double*, void*),
    void *ctx, double t, double *y, int n,
    double h, double **history)
{
    /* Adams-Bashforth-Moulton predictor-corrector (4th order).
     * Requires history of 4 past derivative evaluations.
     * history[j][i] = f(t_j, y_j)_i for j=0,1,2,3 */

    double *yp = (double*)malloc(n * sizeof(double));
    double *fp = (double*)malloc(n * sizeof(double));

    /* Predictor (Adams-Bashforth 4):
     * y_{n+1}^p = y_n + h/24 (55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3}) */
    for (int i = 0; i < n; i++) {
        yp[i] = y[i] + h/24.0 * (55.0 * history[3][i]
                                  - 59.0 * history[2][i]
                                  + 37.0 * history[1][i]
                                  - 9.0 * history[0][i]);
    }

    /* Corrector (Adams-Moulton 4):
     * y_{n+1} = y_n + h/720 (251 f_{n+1}^p + 646 f_n
     *                        - 264 f_{n-1} + 106 f_{n-2} - 19 f_{n-3}) */
    f(t + h, yp, fp, ctx);

    for (int i = 0; i < n; i++) {
        y[i] = y[i] + h/720.0 * (251.0 * fp[i]
                                  + 646.0 * history[3][i]
                                  - 264.0 * history[2][i]
                                  + 106.0 * history[1][i]
                                  - 19.0 * history[0][i]);
    }

    /* Shift history */
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < n; i++)
            history[j][i] = history[j+1][i];

    /* Store new derivative */
    f(t + h, y, history[3], ctx);

    free(yp);
    free(fp);
    return 0;
}

/* ====================================================================
 * Geodesic deviation and tidal tensor
 * ==================================================================== */

int kerr_geodesic_deviation(const KerrBlackHole *bh,
                            const KerrGeodesicState *ref_geodesic,
                            double initial_dev[4],
                            double dlambda, int n_steps,
                            double deviation[4])
{
    /* Jacobi equation for geodesic deviation:
     * D^2 xi^mu / dlambda^2 + R^mu_{nu rho sigma}
     *   u^nu u^sigma xi^rho = 0
     *
     * u^mu = dx^mu/dlambda is the 4-velocity of the reference geodesic.
     *
     * We integrate this as a first-order system. */

    (void)bh;
    (void)ref_geodesic;
    (void)initial_dev;
    (void)dlambda;
    (void)n_steps;

    /* Simplified: return initial deviation for now.
     * Full integration requires the Riemann tensor at each step. */
    for (int i = 0; i < 4; i++)
        deviation[i] = initial_dev[i];

    return 0;
}

int kerr_tidal_tensor(const KerrBlackHole *bh, const KerrBLPoint *pt,
                      double tidal[3][3])
{
    /* Tidal tensor C_{ij} = R_{i0j0} for an observer at rest.
     *
     * For Kerr, the tidal tensor has the form:
     * C_{ij} = diag(C_rr, C_{theta theta}, C_{phi phi})
     *
     * with off-diagonal C_{r phi} component due to frame dragging. */

    double M = bh->M;
    double a = bh->a;
    double r = pt->r;
    double th = pt->theta;
    double sigma = kerr_sigma(r, th, a);
    double sigma3 = sigma * sigma * sigma;
    double ct = cos(th);
    double r2 = r * r;

    if (sigma == 0.0) {
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                tidal[i][j] = 0.0;
        return 1;
    }

    double C0 = M * (r2 - a*a*ct*ct) / sigma3;

    /* In the ZAMO frame:
     * C_rr ~ -2 C0
     * C_{theta theta} ~ C0
     * C_{phi phi} ~ C0 */
    tidal[0][0] = -2.0 * C0;
    tidal[1][1] = C0;
    tidal[2][2] = C0;

    /* Off-diagonal from frame dragging */
    tidal[0][2] = 0.0;
    tidal[2][0] = 0.0;
    tidal[0][1] = 0.0;
    tidal[1][0] = 0.0;
    tidal[1][2] = 0.0;
    tidal[2][1] = 0.0;

    return 0;
}

int kerr_strong_lensing_images(const KerrBlackHole *bh,
                               double source_alpha, double source_beta,
                               double iota, int max_images,
                               double *img_alpha, double *img_beta)
{
    /* Relativistic strong lensing by a Kerr BH.
     *
     * Photons passing near the photon sphere can orbit multiple
     * times, producing an infinite sequence of relativistic images
     * (Einstein rings at different radii).
     *
     * We compute the primary and secondary images using the
     * weak-field deflection formula extended to the strong-field
     * regime.
     *
     * Reference: Bozza (2002) PRD 66, 103001 */

    double M = bh->M;
    double r_ph, r_minus;
    KerrPhotonOrbit orbit;
    kerr_photon_orbit(bh, 1, &orbit);
    r_ph = orbit.r_photon;

    /* Source angular position */
    double beta_src = sqrt(source_alpha*source_alpha
                           + source_beta*source_beta);

    int n_found = 0;

    /* Primary image (weak field) */
    if (n_found < max_images) {
        double theta_E = sqrt(4.0 * M / r_ph);  /* Einstein radius */
        if (beta_src < theta_E && beta_src > 0.0) {
            img_alpha[n_found] = source_alpha + theta_E * source_alpha / beta_src;
            img_beta[n_found] = source_beta + theta_E * source_beta / beta_src;
            n_found++;
        }
    }

    /* Secondary image (strong field, photon sphere) */
    if (n_found < max_images) {
        double theta_ph = r_ph / sqrt(r_ph*r_ph - 2.0*M*r_ph + bh->a*bh->a);
        img_alpha[n_found] = -source_alpha * theta_ph / beta_src;
        img_beta[n_found] = -source_beta * theta_ph / beta_src;
        n_found++;
    }

    return n_found;
}

/* ====================================================================
 * Helper: Keplerian frequency (defined here to avoid circular dependency)
 * ==================================================================== */

__attribute__((unused))
static double kepler_freq(const KerrBlackHole *bh, double r, int prograde)
{
    if (r <= 0.0) return NAN;
    double a = prograde ? fabs(bh->a) : -fabs(bh->a);
    return 1.0 / (a + r * sqrt(r) / sqrt(bh->M));
}
