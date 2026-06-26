/**
 * @file coordinate.c
 * @brief Coordinate systems and metric realizations in GR.
 *
 * Reference: Wald (1984), Ch.6; Carroll (2004), Ch.5
 *            Chandrasekhar "The Mathematical Theory of Black Holes" (1983)
 *            Misner, Thorne, Wheeler "Gravitation" Ch.31-33
 *
 * Exact solutions to the Einstein equations in various coordinate systems.
 */

#include "coordinate.h"
#include <math.h>
#include <string.h>

/* ========================================================================
 * L3: Coordinate transformation utilities
 * ========================================================================*/

void transform_tensor2_covariant(const Tensor2 J_inv, const Tensor2 T,
                                  Tensor2 T_prime)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    sum += J_inv[mu][alpha]
                         * J_inv[nu][beta]
                         * T[alpha][beta];
                }
            }
            T_prime[mu][nu] = sum;
        }
    }
}

void transform_christoffel(const Tensor2 J, const Tensor2 J_inv,
                            const Tensor3 H, const Tensor3 Gamma,
                            Tensor3 Gamma_p)
{
    for (int lambda = 0; lambda < 4; lambda++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                double sum = 0.0;

                /* Homogeneous part: J^lambda_rho * J_inv^sigma_mu * J_inv^tau_nu * Gamma^rho_{sigma,tau} */
                for (int rho = 0; rho < 4; rho++) {
                    for (int sigma = 0; sigma < 4; sigma++) {
                        for (int tau = 0; tau < 4; tau++) {
                            sum += J[lambda][rho]
                                 * J_inv[sigma][mu]
                                 * J_inv[tau][nu]
                                 * Gamma[rho][sigma][tau];
                        }
                    }
                }

                /* Inhomogeneous part: J^lambda_rho * H^rho_{mu,nu} */
                for (int rho = 0; rho < 4; rho++) {
                    sum += J[lambda][rho] * H[rho][mu][nu];
                }

                Gamma_p[lambda][mu][nu] = sum;
            }
        }
    }
}

void compute_jacobian(const double x[4], CoordinateTransform transform,
                       double h, Tensor2 J)
{
    double xp_center[4];
    transform(x, xp_center);

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double x_plus[4], x_minus[4], xp_plus[4], xp_minus[4];
            for (int i = 0; i < 4; i++) {
                x_plus[i] = x[i];
                x_minus[i] = x[i];
            }
            x_plus[nu] += h;
            x_minus[nu] -= h;

            transform(x_plus, xp_plus);
            transform(x_minus, xp_minus);

            /* Central difference */
            J[mu][nu] = (xp_plus[mu] - xp_minus[mu]) / (2.0 * h);
        }
    }
}

/* ========================================================================
 * L6: Schwarzschild metric
 * ========================================================================*/

/**
 * Schwarzschild metric (exterior solution, r > 2M).
 *
 * This is the unique static, spherically symmetric vacuum solution
 * (Birkhoff's theorem 1923).
 *
 * The metric has two singularities:
 *   r = 2M: coordinate singularity (event horizon)
 *   r = 0:  physical singularity (curvature diverges, K ~ 48M²/r⁶)
 *
 * Knowledge point: The Schwarzschild solution (1916) was the first exact
 * solution to Einstein's equations, discovered only months after the
 * theory was completed. It describes black holes and the exterior of
 * any spherically symmetric body.
 *
 * Birkhoff's theorem: The Schwarzschild metric is the unique spherically
 * symmetric vacuum solution — even time-dependent spherical sources cannot
 * emit gravitational waves (they produce no time-varying quadrupole moment).
 */
void schwarzschild_metric(double M, double r, double theta, Metric *m)
{
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            m->g[mu][nu] = 0.0;

    double f = 1.0 - 2.0 * M / r;

    m->g[0][0] = -f;           /* g_{tt} */
    m->g[1][1] = 1.0 / f;      /* g_{rr} */
    m->g[2][2] = r * r;        /* g_{θθ} */
    m->g[3][3] = r * r * sin(theta) * sin(theta); /* g_{φφ} */

    metric_compute_inverse(m);
}

/**
 * Eddington-Finkelstein (ingoing): regular at horizon.
 *
 * Transform t → v = t + r_* where r_* = r + 2M ln|r/(2M) - 1|.
 * The metric becomes:
 *   ds² = -(1 - 2M/r) dv² + 2 dv dr + r² dΩ²
 *
 * Knowledge point: Eddington-Finkelstein coordinates extend the
 * Schwarzschild solution across the event horizon, revealing that
 * r = 2M is a regular null surface (not a true singularity).
 * Ingoing coordinates cover the black hole interior (region II).
 */
void eddington_finkelstein_ingoing(double M, double r, double theta, Metric *m)
{
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            m->g[mu][nu] = 0.0;

    double f = 1.0 - 2.0 * M / r;

    /* Coordinates: (v, r, θ, φ) */
    m->g[0][0] = -f;       /* g_{vv} */
    m->g[0][1] =  1.0;     /* g_{vr} */
    m->g[1][0] =  1.0;     /* g_{rv} */
    m->g[1][1] =  0.0;     /* g_{rr} */
    m->g[2][2] =  r * r;
    m->g[3][3] =  r * r * sin(theta) * sin(theta);

    metric_compute_inverse(m);
}

void eddington_finkelstein_outgoing(double M, double r, double theta, Metric *m)
{
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            m->g[mu][nu] = 0.0;

    double f = 1.0 - 2.0 * M / r;

    /* Coordinates: (u, r, θ, φ) — u = t - r_* (retarded time) */
    m->g[0][0] = -f;       /* g_{uu} */
    m->g[0][1] = -1.0;     /* g_{ur} */
    m->g[1][0] = -1.0;     /* g_{ru} */
    m->g[1][1] =  0.0;     /* g_{rr} */
    m->g[2][2] =  r * r;
    m->g[3][3] =  r * r * sin(theta) * sin(theta);

    metric_compute_inverse(m);
}

/**
 * Kruskal-Szekeres maximal extension.
 *
 * Transform to Kruskal coordinates (T, X):
 *   T = sqrt(r/(2M) - 1) exp(r/(4M)) sinh(t/(4M))  (exterior)
 *   X = sqrt(r/(2M) - 1) exp(r/(4M)) cosh(t/(4M))
 *
 * The line element:
 *   ds² = (32M³/r) exp(-r/(2M)) (-dT² + dX²) + r² dΩ²
 *
 * Knowledge point: Kruskal-Szekeres coordinates (1960) give the maximal
 * analytic extension of Schwarzschild, revealing four regions:
 *   I:  exterior (our universe)
 *   II: black hole interior (future trapped)
 *   III: white hole interior (past trapped)
 *   IV: parallel exterior (Einstein-Rosen bridge)
 */
void kruskal_szekeres_metric(double M, double T, double X, double theta,
                              Metric *m)
{
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            m->g[mu][nu] = 0.0;

    /* Need r implicitly from T² - X² = (1 - r/(2M)) exp(r/(2M)).
     * Since there's no closed-form inverse, the caller must provide
     * consistent (T, X) such that r can be determined.
     * For numerical simplicity, we assume r is passed indirectly.
     * Here we use a Newton iteration. */

    /* Initial guess for r: */
    double val = T*T - X*X;
    double r;
    if (val > 0) {
        /* Exterior: val = (r/(2M) - 1) exp(r/(2M)) — wait.
         * Correct: 1 - r/(2M) = val * exp(-r/(2M))
         * So: r/(2M) = 1 - val * exp(-r/(2M)) — implicit. */
        r = 2.0 * M * (1.0 + val); /* rough guess for exterior */
    } else {
        /* Interior */
        r = 2.0 * M * 0.5;
    }

    /* Simple Newton iteration for r:
     * f(r) = (1 - r/(2M)) exp(r/(2M)) - (T² - X²) = 0
     * Actually the standard relation: T² - X² = (1 - r/(2M)) exp(r/(2M)) */
    for (int iter = 0; iter < 50; iter++) {
        double r_over_2M = r / (2.0 * M);
        double exp_term = exp(r_over_2M);
        double f = (1.0 - r_over_2M) * exp_term - (T*T - X*X);
        double df = -1.0/(2.0*M) * exp_term + (1.0 - r_over_2M) * exp_term / (2.0*M);
        double dr = f / df;
        r -= dr;
        if (fabs(dr) < 1e-12) break;
    }
    if (r < 1e-10) r = 2.0 * M; /* floor */

    double conformal = (32.0 * M*M*M / r) * exp(-r / (2.0 * M));

    m->g[0][0] = -conformal;  /* g_{TT} */
    m->g[1][1] =  conformal;  /* g_{XX} */
    m->g[2][2] =  r * r;
    m->g[3][3] =  r * r * sin(theta) * sin(theta);

    metric_compute_inverse(m);
}

/**
 * Tortoise coordinate: r_* = r + 2M ln|r/(2M) - 1|
 *
 * Maps r in (2M, infinity) to r_* in (-infinity, infinity).
 *
 * Knowledge point: The tortoise coordinate stretches the radial coordinate
 * so that radial null geodesics satisfy d(r_*)/dt = +/- 1 (ingoing/outgoing).
 * This makes the causal structure transparent with 45-degree light cones.
 */
double tortoise_coordinate(double r, double M)
{
    if (r <= 2.0 * M) return -1e100; /* inside horizon */
    return r + 2.0 * M * log(fabs(r / (2.0 * M) - 1.0));
}

double inverse_tortoise_coordinate(double rstar, double M)
{
    /* Solve r + 2M ln(r/(2M) - 1) = rstar via Newton */
    if (rstar > 100.0 * M) return rstar; /* asymptotic: r ≈ rstar */

    double r = rstar;
    if (r < 2.001 * M) r = 2.001 * M;

    for (int i = 0; i < 100; i++) {
        double f = r + 2.0*M * log(r/(2.0*M) - 1.0) - rstar;
        double df = 1.0 + 2.0*M / (r - 2.0*M);
        double dr = f / df;
        r -= dr;
        if (fabs(dr) < 1e-12) break;
    }
    return r;
}

/* ========================================================================
 * L6: Kerr metric (Boyer-Lindquist coordinates)
 * ========================================================================*/

/**
 * Kerr metric — rotating black hole solution (Roy Kerr, 1963).
 *
 * Parameters:
 *   M: mass
 *   a: spin parameter (J/M), 0 ≤ a ≤ M
 *
 * For a=0: reduces to Schwarzschild
 * For a=M: extremal Kerr (maximal spin)
 *
 * The Kerr solution is axisymmetric but not spherically symmetric.
 *
 * Key features:
 *   - Event horizon: r_+ = M + √(M² - a²)
 *   - Cauchy horizon: r_- = M - √(M² - a²)
 *   - Ergosphere: g_{tt} = 0 at r = M + √(M² - a² cos²θ)
 *   - Ring singularity at r=0, θ=π/2
 *
 * Knowledge point: The Kerr metric (1963) is the unique stationary,
 * axisymmetric vacuum black hole solution (Kerr-Newman uniqueness
 * theorem: black holes are characterized only by M, J, Q).
 * Rotating black holes are astrophysically relevant — most black holes
 * in nature are expected to be Kerr.
 */
void kerr_metric(double M, double a, double r, double theta, Metric *m)
{
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            m->g[mu][nu] = 0.0;

    double st = sin(theta);
    double ct = cos(theta);

    double rho2 = r*r + a*a * ct*ct;       /* ρ² = r² + a² cos²θ */
    double Delta = r*r - 2.0*M*r + a*a;    /* Δ = r² - 2Mr + a² */

    /* g_{tt} = -(1 - 2Mr/ρ²) */
    m->g[0][0] = -(1.0 - 2.0 * M * r / rho2);

    /* g_{tφ} = g_{φt} = -2Mar sin²θ / ρ² */
    m->g[0][3] = -2.0 * M * a * r * st * st / rho2;
    m->g[3][0] = m->g[0][3];

    /* g_{rr} = ρ² / Δ */
    m->g[1][1] = rho2 / Delta;

    /* g_{θθ} = ρ² */
    m->g[2][2] = rho2;

    /* g_{φφ} = (r² + a² + 2Ma²r sin²θ/ρ²) sin²θ */
    m->g[3][3] = (r*r + a*a + 2.0*M*a*a*r*st*st/rho2) * st * st;

    metric_compute_inverse(m);
}

double kerr_ergosphere_radius(double M, double a, double theta)
{
    double ct = cos(theta);
    return M + sqrt(M*M - a*a * ct*ct);
}

double kerr_horizon_radius(double M, double a)
{
    return M + sqrt(M*M - a*a);
}

double kerr_cauchy_horizon_radius(double M, double a)
{
    return M - sqrt(M*M - a*a);
}

/* ========================================================================
 * L7: Schwarzschild in isotropic coordinates
 * ========================================================================*/

/**
 * Schwarzschild in isotropic coordinates (r̄, θ, φ):
 *   ds² = -[(1-M/(2r̄))/(1+M/(2r̄))]² dt² + [1+M/(2r̄)]⁴ (dr̄² + r̄² dΩ²)
 *
 * The spatial part is conformally flat, which is useful for initial
 * data construction in numerical relativity (Bowen-York puncture method).
 *
 * Relation to Schwarzschild r: r = r̄ [1 + M/(2r̄)]²
 *
 * Knowledge point: Isotropic coordinates make the spatial metric manifestly
 * conformally flat, simplifying the constraint equations in the initial
 * value problem of numerical relativity.
 */
void schwarzschild_isotropic_metric(double M, double r_bar, double theta,
                                     Metric *m)
{
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            m->g[mu][nu] = 0.0;

    double M_over_2r = M / (2.0 * r_bar);
    double psi = 1.0 + M_over_2r;          /* conformal factor */
    double psi2 = psi * psi;
    double psi4 = psi2 * psi2;

    double lapse = (1.0 - M_over_2r) / psi; /* sqrt of g_{tt} prefactor */
    m->g[0][0] = -lapse * lapse;

    double spatial_factor = psi4;
    m->g[1][1] = spatial_factor;
    m->g[2][2] = spatial_factor * r_bar * r_bar;
    m->g[3][3] = spatial_factor * r_bar * r_bar * sin(theta) * sin(theta);

    metric_compute_inverse(m);
}
