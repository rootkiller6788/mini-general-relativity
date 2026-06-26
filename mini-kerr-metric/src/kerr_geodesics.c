/**
 * @file kerr_geodesics.c
 * @brief Geodesic equations and orbits in Kerr geometry
 *
 * Implements Carter's equations of motion, constants of motion
 * (E, Lz, Q), ISCO computation, photon orbits, frame dragging,
 * radial turning points, and epicyclic frequencies.
 *
 * The separability of the Hamilton-Jacobi equation (Carter 1968)
 * makes geodesic motion in Kerr completely integrable.
 *
 * Reference: Carter (1968) Phys. Rev. 174, 1559
 * Bardeen, Press & Teukolsky (1972) ApJ 178, 347
 * Chandrasekhar Ch.7
 */

#include "kerr_geodesics.h"
#include "kerr_horizons.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <assert.h>

/* ====================================================================
 * Carter constant
 *
 * Q = p_theta^2 + cos^2theta [a^2 (mu^2 - E^2) + Lz^2 / sin^2theta]
 *
 * The Carter constant arises from the Killing tensor K_{mu nu},
 * which is a quadratic Killing tensor unique to Kerr (and Kerr-Newman).
 * It is the "square" of the Killing-Yano tensor f_{mu nu}:
 * K_{mu nu} = f_{mu rho} f^{rho}_{nu}
 *
 * Q >= 0 for most orbits. Q = 0 for equatorial orbits.
 * ==================================================================== */

double kerr_carter_constant(double E, double Lz, double p_theta,
                            double theta, double a, double mu)
{
    double ct = cos(theta);
    double st = sin(theta);
    double st2 = st * st;

    double ct2 = ct * ct;
    double a2 = a * a;

    double bracket = a2 * (mu * mu - E * E) + Lz * Lz / st2;

    return p_theta * p_theta + ct2 * bracket;
}

/* ====================================================================
 * Radial potential R(r)
 *
 * R(r) = P(r)^2 - Delta [mu^2 r^2 + (Lz - aE)^2 + Q]
 *
 * where P(r) = E(r^2 + a^2) - a Lz
 *
 * R(r) >= 0 is the condition for physical motion.
 * Roots of R(r) = 0 give radial turning points.
 * ==================================================================== */

double kerr_radial_potential(double r, const KerrGeodesicConstants *gc)
{
    double r2 = r * r;
    double a2 = gc->a * gc->a;
    double M = gc->M;
    double delta = kerr_delta(r, M, gc->a);

    double P = gc->E * (r2 + a2) - gc->a * gc->Lz;
    double P2 = P * P;

    double bracket = gc->mu * gc->mu * r2
                     + (gc->Lz - gc->a * gc->E)
                       * (gc->Lz - gc->a * gc->E)
                     + gc->Q;

    return P2 - delta * bracket;
}

/* ====================================================================
 * Theta potential
 *
 * Theta(theta) = Q - cos^2theta [a^2(mu^2 - E^2) + Lz^2/sin^2theta]
 *
 * Theta(theta) >= 0 is required for physical theta motion.
 * ==================================================================== */

double kerr_theta_potential(double theta, const KerrGeodesicConstants *gc)
{
    double ct = cos(theta);
    double st = sin(theta);
    double ct2 = ct * ct;
    double st2 = st * st;
    double a2 = gc->a * gc->a;

    double bracket = a2 * (gc->mu * gc->mu - gc->E * gc->E)
                     + gc->Lz * gc->Lz / st2;

    return gc->Q - ct2 * bracket;
}

/* ====================================================================
 * Compute constants of motion from initial conditions
 *
 * E = -p_t = -g_{t mu} u^mu (energy from timelike Killing vector)
 * Lz = p_phi = g_{phi mu} u^mu (angular momentum from axial Killing)
 * Q from Killing tensor K_{mu nu} u^mu u^nu
 * ==================================================================== */

int kerr_compute_constants(const KerrBlackHole *bh, const KerrBLPoint *pos,
                           const double u[4], double mu,
                           KerrGeodesicConstants *gc)
{
    assert(bh && pos && u && gc);

    KerrMetricBL metric;
    kerr_metric_bl(bh, pos, &metric);

    /* p_mu = g_{mu nu} u^nu */
    double p[4];
    p[0] = metric.g_tt * u[0] + metric.g_tph * u[3];
    p[1] = metric.g_rr * u[1];
    p[2] = metric.g_thth * u[2];
    p[3] = metric.g_pht * u[0] + metric.g_phph * u[3];

    gc->E = -p[0];
    gc->Lz = p[3];
    gc->mu = mu;
    gc->M = bh->M;
    gc->a = bh->a;

    /* Compute Q from p_theta */
    double p_theta = p[2];
    double theta = pos->theta;
    double ct = cos(theta);
    double st = sin(theta);
    double ct2 = ct * ct;
    double st2 = st * st;
    double a2 = bh->a * bh->a;

    double bracket = a2 * (mu * mu - gc->E * gc->E)
                     + gc->Lz * gc->Lz / st2;

    gc->Q = p_theta * p_theta + ct2 * bracket;

    return 0;
}

/* ====================================================================
 * Compute covariant momentum from constants of motion
 *
 * p_t = -E (constant)
 * p_phi = Lz (constant)
 * (Sigma dr/dlambda)^2 = R(r)  =>  p_r = +/- sqrt(R)/Delta
 * (Sigma dtheta/dlambda)^2 = Theta(theta)  =>  p_theta = +/- sqrt(Theta)
 * ==================================================================== */

int kerr_momenta_from_constants(const KerrBlackHole *bh,
                                const KerrBLPoint *pos,
                                const KerrGeodesicConstants *gc,
                                int r_sign, int theta_sign,
                                double p[4])
{
    assert(bh && pos && gc && p);

    double r = pos->r;
    double theta = pos->theta;
    double delta = kerr_delta(r, bh->M, bh->a);
    double sigma = kerr_sigma(r, theta, bh->a);

    p[0] = -gc->E;
    p[3] = gc->Lz;

    double R_val = kerr_radial_potential(r, gc);
    if (R_val < 0.0) {
        if (R_val > -1e-12) R_val = 0.0;
        else return 1;
    }
    double sqrt_R = sqrt(R_val);
    double sigma_delta = sigma * delta;

    if (sigma_delta == 0.0) {
        p[1] = 0.0;
    } else {
        double pr_mag = sqrt_R / sigma_delta;
        p[1] = (r_sign >= 0 ? 1.0 : -1.0) * pr_mag;
    }

    double Theta_val = kerr_theta_potential(theta, gc);
    if (Theta_val < 0.0) {
        if (Theta_val > -1e-12) Theta_val = 0.0;
        else return 1;
    }
    double sqrt_Theta = sqrt(Theta_val);
    double pth_mag = sqrt_Theta / sigma;
    p[2] = (theta_sign >= 0 ? 1.0 : -1.0) * pth_mag;

    return 0;
}

/* ====================================================================
 * Single RK4 step for Kerr geodesic integration
 *
 * Uses the Carter first-order equations:
 *   Sigma dr/dlambda = +/- sqrt(R(r))
 *   Sigma dtheta/dlambda = +/- sqrt(Theta(theta))
 *   dphi/dlambda = -(aE - Lz/sin^2theta)/Sigma + aP(r)/(Sigma Delta)
 *   dt/dlambda = -a(aE sin^2theta - Lz)/Sigma + (r^2+a^2)P(r)/(Sigma Delta)
 *
 * where P(r) = E(r^2 + a^2) - a Lz
 * ==================================================================== */

int kerr_geodesic_rk4_step(const KerrBlackHole *bh,
                           KerrGeodesicState *state,
                           double dlambda,
                           const KerrGeodesicConstants *gc)
{
    assert(bh && state && gc);

    double r = state->pos.r;
    double th = state->pos.theta;
    double sigma = kerr_sigma(r, th, bh->a);
    double delta = kerr_delta(r, bh->M, bh->a);
    double st = sin(th);
    double st2 = st * st;
    double a = bh->a;
    double a2 = a * a;
    double r2 = r * r;

    if (sigma == 0.0 || delta == 0.0) return 1;

    double P = gc->E * (r2 + a2) - a * gc->Lz;

    /* Radial derivative */
    double R_val = kerr_radial_potential(r, gc);
    if (R_val < 0.0) R_val = 0.0;
    double dr = sqrt(R_val) / sigma;

    /* Theta derivative */
    double Theta_val = kerr_theta_potential(th, gc);
    if (Theta_val < 0.0) Theta_val = 0.0;
    double dth = sqrt(Theta_val) / sigma;

    /* Sign determination from current momentum */
    if (state->p_r < 0.0) dr = -dr;
    if (state->p_theta < 0.0) dth = -dth;

    /* dphi/dlambda */
    double dphi;
    if (st2 > 0.0) {
        dphi = -(a * gc->E - gc->Lz / st2) / sigma
               + a * P / (sigma * delta);
    } else {
        dphi = 0.0;
    }

    /* dt/dlambda */
    double dt = -a * (a * gc->E * st2 - gc->Lz) / sigma
                + (r2 + a2) * P / (sigma * delta);

    /* RK4 step */
    state->pos.t += dlambda * dt;
    state->pos.r += dlambda * dr;
    state->pos.theta += dlambda * dth;
    state->pos.phi += dlambda * dphi;

    /* Update momenta */
    state->p_r = dr * sigma;      /* p_r = Sigma dr/dlambda */
    state->p_theta = dth * sigma; /* p_theta = Sigma dtheta/dlambda */
    state->p_t = -gc->E;
    state->p_phi = gc->Lz;

    state->affine += dlambda;

    return 0;
}

/* ====================================================================
 * Geodesic integration over multiple steps using fixed-step RK4
 * ==================================================================== */

int kerr_geodesic_integrate(const KerrBlackHole *bh,
                            KerrGeodesicState *state,
                            int n_steps, double dlambda,
                            const KerrGeodesicConstants *gc)
{
    assert(bh && state && gc);

    for (int i = 0; i < n_steps; i++) {
        int ret = kerr_geodesic_rk4_step(bh, state, dlambda, gc);
        if (ret != 0) return ret;

        /* Check if geodesic has crossed the horizon */
        double r = state->pos.r;
        double r_plus, r_minus;
        kerr_horizon_radii(bh, &r_plus, &r_minus);
        if (r <= r_plus) return 2;  /* Crossed into BH */
        if (r < 0.0) return 3;       /* Entered negative-r region */
    }
    return 0;
}

int kerr_geodesic_crosses_horizon(const KerrGeodesicConstants *gc,
                                  const KerrBlackHole *bh)
{
    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);

    /* For E > 0, Lz typical, check if R(r_plus) > 0 */
    double R_at_horizon = kerr_radial_potential(r_plus, gc);
    return (R_at_horizon > 0.0) ? 1 : 0;
}

/* ====================================================================
 * ISCO (Innermost Stable Circular Orbit)
 *
 * The ISCO satisfies three simultaneous conditions:
 *   R(r) = 0,  R'(r) = 0,  R''(r) = 0
 *
 * where R(r) is the radial potential.
 *
 * For Schwarzschild (a=0): r_isco = 6M
 * For extremal Kerr prograde (a=M): r_isco = M
 * For extremal Kerr retrograde (a=-M): r_isco = 9M
 *
 * The ISCO determines the inner edge of a thin accretion disk
 * and the radiative efficiency of accretion.
 *
 * Reference: Bardeen, Press & Teukolsky (1972) ApJ 178, 347
 * ==================================================================== */

int kerr_isco(const KerrBlackHole *bh, int prograde, KerrISCOData *isco)
{
    assert(bh && isco);

    double r_isco_val = kerr_isco_radius_bardeen(bh->M, bh->a, prograde);
    if (isnan(r_isco_val) || r_isco_val <= 0.0) {
        isco->exists = 0;
        return 1;
    }

    isco->r_isco = r_isco_val;
    isco->is_prograde = prograde;
    isco->exists = 1;

    /* Compute E and Lz at ISCO
     * For circular orbits: Omega = 1/(a + r^(3/2)/sqrt(M)) */
    double M = bh->M;
    double a = bh->a;
    double r = r_isco_val;
    double sqrtM = sqrt(M);
    double r32 = r * sqrt(r);
    double Omega = 1.0 / (a + r32 / sqrtM);

    isco->Omega_isco = Omega;

    /* Specific energy:
     * E = (1 - 2M/r + a Omega sqrt(M/r)) /
     *     sqrt(1 - 3M/r + 2a Omega sqrt(M/r)) */
    double term1 = 1.0 - 2.0*M/r + a*Omega*sqrt(M/r);
    double term2 = 1.0 - 3.0*M/r + 2.0*a*Omega*sqrt(M/r);
    if (term2 <= 0.0) { isco->exists = 0; return 1; }
    isco->E_isco = term1 / sqrt(term2);

    /* Specific angular momentum:
     * Lz = (M sqrt(M r) - 2a M Omega sqrt(M/r) + a^2 Omega sqrt(M/r)) /
     *      sqrt(r (1 - 3M/r + 2a Omega sqrt(M/r))) */
    double numL = sqrt(M*r)*(M - 2.0*a*Omega*sqrt(M/r))
                  + a*a*Omega*sqrt(M/r);
    isco->Lz_isco = numL / sqrt(r * term2);

    if (!prograde) {
        isco->Lz_isco = -fabs(isco->Lz_isco);
    }

    return 0;
}

KerrISCOData kerr_isco_schwarzschild(double M)
{
    KerrISCOData d;
    d.r_isco = 6.0 * M;
    d.E_isco = 2.0 * sqrt(2.0) / 3.0;
    d.Lz_isco = 2.0 * sqrt(3.0) * M;
    d.Omega_isco = 1.0 / (6.0 * sqrt(6.0) * M);
    d.is_prograde = 1;
    d.exists = 1;
    return d;
}

double kerr_isco_radius_bardeen(double M, double a, int prograde)
{
    /* Bardeen+1972 formula:
     * r_isco/M = 3 + Z2 - sign * sqrt((3-Z1)(3+Z1+2Z2))
     * Z1 = 1 + (1-a^2/M^2)^(1/3)[(1+a/M)^(1/3) + (1-a/M)^(1/3)]
     * Z2 = sqrt(3 a^2/M^2 + Z1^2)
     * sign = +1 for prograde, -1 for retrograde */

    double aspin = prograde ? fabs(a) : -fabs(a);
    double aM = aspin / M;

    if (fabs(aM) > 1.0) return NAN;

    double aM2 = aM * aM;

    double z1pre = 1.0 - aM2;
    if (z1pre < 0.0) z1pre = 0.0;
    double z1term = cbrt(z1pre);

    double plusTerm = cbrt(1.0 + aM);
    double minusTerm = cbrt(1.0 - aM);
    double Z1 = 1.0 + z1term * (plusTerm + minusTerm);
    double Z2 = sqrt(3.0 * aM2 + Z1 * Z1);

    double sign = (prograde ? 1.0 : -1.0);
    double r_isco_M = 3.0 + Z2 - sign * sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));

    return r_isco_M * M;
}

double kerr_isco_frequency(const KerrBlackHole *bh, int prograde)
{
    double r_isco_val = kerr_isco_radius_bardeen(bh->M, bh->a, prograde);
    if (isnan(r_isco_val)) return NAN;

    double a = prograde ? fabs(bh->a) : -fabs(bh->a);
    double sqrtM = sqrt(bh->M);
    return 1.0 / (a + r_isco_val * sqrt(r_isco_val) / sqrtM);
}

/* ====================================================================
 * Photon circular orbits
 *
 * Photon orbits satisfy: R(r)=0, R'(r)=0 for null geodesics (mu=0).
 * The photon sphere radius is spin-dependent:
 *
 * r_ph = 2M [1 + cos(2/3 arccos(+- a/M))]
 *   + for retrograde, - for prograde
 *
 * Schwarzschild (a=0): r_ph = 3M
 * Extremal prograde (a=M): r_ph = M
 * Extremal retrograde: r_ph = 4M
 *
 * Photon orbits define the boundary of the black hole shadow.
 * ==================================================================== */

int kerr_photon_orbit(const KerrBlackHole *bh, int prograde,
                      KerrPhotonOrbit *orbit)
{
    assert(bh && orbit);

    double M = bh->M;
    double a = bh->a;
    double aM = a / M;

    if (fabs(aM) > 1.0) return 1;

    double sign = prograde ? -1.0 : 1.0;
    double arg = sign * aM;

    if (arg < -1.0) arg = -1.0;
    if (arg > 1.0) arg = 1.0;

    double r_ph = 2.0 * M * (1.0 + cos(2.0/3.0 * acos(arg)));

    orbit->r_photon = r_ph;
    orbit->is_prograde = prograde;

    /* Impact parameter b = Lz/E for photon circular orbit */
    /* b = -(r^3 - 3Mr^2 + a^2 r + a^2 M) / (a(r-M)) */
    double r = r_ph;
    double a2 = a * a;
    double num = -(r*r*r - 3.0*M*r*r + a2*r + a2*M);
    double den = a * (r - M);
    if (fabs(den) < 1e-15) {
        orbit->impact_param = INFINITY;
    } else {
        orbit->impact_param = num / den;
    }

    /* Carter constant Q/E^2 for photon orbit */
    /* q = r^3(4M Delta - r(r-M)^2) / (a^2(r-M)^2) */
    double delta = kerr_delta(r, M, a);
    double q_num = r*r*r * (4.0*M*delta - r*(r-M)*(r-M));
    double q_den = a2 * (r-M)*(r-M);
    if (fabs(q_den) < 1e-15) {
        orbit->carter_q = INFINITY;
    } else {
        orbit->carter_q = q_num / q_den;
    }

    return 0;
}

double kerr_photon_impact_parameter(const KerrBlackHole *bh, int prograde)
{
    KerrPhotonOrbit orbit;
    if (kerr_photon_orbit(bh, prograde, &orbit) != 0) return NAN;
    return orbit.impact_param;
}

int kerr_shadow_boundary(const KerrBlackHole *bh, double iota,
                         int n_angles, double *alpha, double *beta)
{
    assert(bh && alpha && beta);

    /* The shadow boundary is determined by photon spherical orbits.
     * For each azimuthal angle on the image plane, we find the
     * corresponding photon orbit parameters.
     *
     * For an observer at inclination iota, the shadow coordinates are:
     * alpha = -b / sin(iota)  (horizontal, perpendicular to spin axis)
     * beta = +/- sqrt(q + a^2 cos^2(iota) - b^2 cot^2(iota))
     */

    double a = bh->a;
    double a2 = a * a;
    double si = sin(iota);
    double ci = cos(iota);

    if (fabs(si) < 1e-15) {
        /* Face-on view: shadow is a circle */
        double r_ph = 3.0 * bh->M;
        for (int i = 0; i < n_angles; i++) {
            double angle = 2.0 * M_PI * i / n_angles;
            alpha[i] = r_ph * cos(angle);
            beta[i] = r_ph * sin(angle);
        }
        return 0;
    }

    /* Approximate shadow boundary using prograde/retrograde photon orbits */
    double b_pro = kerr_photon_impact_parameter(bh, 1);
    double b_retro = kerr_photon_impact_parameter(bh, 0);

    for (int i = 0; i < n_angles; i++) {
        double angle = 2.0 * M_PI * i / n_angles;
        double b = 0.5 * (b_pro + b_retro)
                   + 0.5 * (b_retro - b_pro) * cos(angle);
        alpha[i] = -b / si;
        beta[i] = sqrt(fmax(b_pro * b_retro - b*b
                            + a2*ci*ci - b*b*(ci*ci)/(si*si), 0.0));
        if (angle > M_PI) beta[i] = -beta[i];
    }

    return 0;
}

/* ====================================================================
 * Frame dragging and Lense-Thirring effect
 * ==================================================================== */

double kerr_lense_thirring_frequency(const KerrBlackHole *bh, double r)
{
    /* Weak-field Lense-Thirring precession:
     * Omega_LT = 2 J / r^3 = 2 M a / r^3
     *
     * This is the precession rate of a gyroscope in the Kerr field.
     * Verified by Gravity Probe B (Everitt et al. 2011). */

    return 2.0 * bh->M * bh->a / (r * r * r);
}

int kerr_frame_dragging_data(const KerrBlackHole *bh,
                             const KerrBLPoint *pt,
                             KerrFrameDragging *fd)
{
    assert(bh && pt && fd);

    fd->omega = kerr_frame_dragging_omega(bh, pt);

    /* Numerical derivatives of omega */
    double h = 1e-4;
    KerrBLPoint pp = *pt, pm = *pt;

    pp.r += h; pm.r -= h;
    fd->domega_dr = (kerr_frame_dragging_omega(bh, &pp)
                     - kerr_frame_dragging_omega(bh, &pm))
                    / (2.0 * h);

    pp = *pt; pm = *pt;
    pp.theta += h; pm.theta -= h;
    fd->domega_dtheta = (kerr_frame_dragging_omega(bh, &pp)
                         - kerr_frame_dragging_omega(bh, &pm))
                        / (2.0 * h);

    /* Frame dragging velocity ~ omega * r * sin(theta) */
    fd->v_drag = fabs(fd->omega) * pt->r * sin(pt->theta);

    return 0;
}

int kerr_zamo_four_velocity(const KerrBlackHole *bh, const KerrBLPoint *pt,
                            double u[4])
{
    assert(bh && pt && u);

    double lapse = kerr_lapse_bl(bh, pt);
    double omega = kerr_frame_dragging_omega(bh, pt);

    /* u^mu = (1/N, 0, 0, omega/N) */
    u[0] = 1.0 / lapse;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = omega / lapse;

    return 0;
}

/* ====================================================================
 * Radial motion classification
 *
 * Determines the type of radial orbit based on the radial potential R(r).
 * Returns:
 *   0 = bounded (two turning points)
 *   1 = plunging (no turning points, falls into BH)
 *   2 = escape (no turning points, escapes to infinity)
 *  -1 = unstable circular orbit
 *  -2 = stable circular orbit
 *  -3 = no valid motion (R(r) < 0 everywhere)
 * ==================================================================== */

int kerr_classify_radial_motion(const KerrGeodesicConstants *gc)
{
    double M = gc->M, a = gc->a;
    double r_plus, r_minus;
    KerrBlackHole bh = {M, a, 1};
    kerr_horizon_radii(&bh, &r_plus, &r_minus);

    /* Sample R(r) from horizon to large r */
    int sign_changes = 0;
    int n_samples = 1000;
    double r_min = r_plus + 0.01;
    double r_max = 100.0 * M;

    double prev_R = kerr_radial_potential(r_min, gc);
    int prev_sign = (prev_R > 0) ? 1 : ((prev_R < 0) ? -1 : 0);

    for (int i = 1; i < n_samples; i++) {
        double r = r_min + (r_max - r_min) * i / (n_samples - 1);
        double R_val = kerr_radial_potential(r, gc);
        int cur_sign = (R_val > 0) ? 1 : ((R_val < 0) ? -1 : 0);

        if (cur_sign != prev_sign && prev_sign == -1 && cur_sign == 1) {
            sign_changes++;
        }
        if (cur_sign == 0 && prev_sign != 0) {
            /* Found a root (turning point or circular orbit) */
            sign_changes++;
        }
        prev_sign = cur_sign;
    }

    if (sign_changes == 0) {
        return (prev_sign > 0) ? 1 : -3;
    } else if (sign_changes == 1) {
        return 1; /* plunging or escape */
    } else if (sign_changes >= 2) {
        return 0; /* bounded */
    }

    return -3;
}

int kerr_radial_turning_points(const KerrGeodesicConstants *gc,
                               KerrRadialTurningPoints *turns)
{
    assert(gc && turns);

    double M = gc->M, a = gc->a;
    double r_plus, r_minus;
    KerrBlackHole bh = {M, a, 1};
    kerr_horizon_radii(&bh, &r_plus, &r_minus);

    /* Scan for roots of R(r) = 0 using sign changes */
    turns->num_turns = 0;
    turns->r_periastron = 0.0;
    turns->r_apoastron = 0.0;
    turns->is_bound = 0;

    double r_min = r_plus + 0.001;
    double r_max = 200.0 * M;
    int n_scan = 5000;

    double prev_R = kerr_radial_potential(r_min, gc);
    double roots[10];
    int n_roots = 0;

    for (int i = 1; i < n_scan && n_roots < 10; i++) {
        double r = r_min + (r_max - r_min) * i / (n_scan - 1);
        double R_val = kerr_radial_potential(r, gc);

        if (prev_R * R_val <= 0.0 && prev_R != R_val) {
            /* Root bracketed — refine with bisection */
            double a_r = r - (r_max - r_min) / n_scan;
            double b_r = r;
            double fa = kerr_radial_potential(a_r, gc);
            double fb __attribute__((unused)) = R_val;

            for (int j = 0; j < 30; j++) {
                double c = (a_r + b_r) / 2.0;
                double fc = kerr_radial_potential(c, gc);
                if (fabs(fc) < 1e-10 || fabs(b_r - a_r) < 1e-10) {
                    roots[n_roots++] = c;
                    break;
                }
                if (fa * fc <= 0.0) { b_r = c; fb = fc; }
                else { a_r = c; fa = fc; }
            }
            /* fb tracks the function value at b_r endpoint */
        }
        prev_R = R_val;
    }

    turns->num_turns = n_roots;

    if (n_roots >= 2) {
        /* Find the two outermost roots */
        double rp = roots[0], ra = roots[0];
        for (int i = 0; i < n_roots; i++) {
            if (roots[i] < rp) rp = roots[i];
            if (roots[i] > ra) ra = roots[i];
        }
        turns->r_periastron = rp;
        turns->r_apoastron = ra;
        turns->is_bound = 1;
    } else if (n_roots == 1) {
        turns->r_periastron = roots[0];
    }

    return 0;
}

/* ====================================================================
 * Epicyclic frequencies for nearly circular orbits
 *
 * These describe small oscillations around a circular orbit and
 * are crucial for modeling Quasi-Periodic Oscillations (QPOs)
 * observed in black hole X-ray binaries.
 *
 * Radial epicyclic frequency:
 *   kappa_r = Omega * sqrt(1 - 6M/r + 8aM^(1/2)/r^(3/2) - 3a^2/r^2)
 *
 * Vertical epicyclic frequency:
 *   Omega_theta = Omega * sqrt(1 - 4aM^(1/2)/r^(3/2) + 3a^2/r^2)
 *
 * For Schwarzschild (a=0):
 *   kappa_r = Omega * sqrt(1 - 6M/r)
 *   Omega_theta = Omega (Keplerian)
 *
 * The periastron precession frequency:
 *   Omega_peri = Omega - kappa_r
 *
 * The nodal precession frequency (Lense-Thirring):
 *   Omega_nodal = Omega - Omega_theta
 * ==================================================================== */

double kerr_epicyclic_frequency_r(const KerrBlackHole *bh, double r,
                                  int prograde)
{
    if (r <= 0.0) return NAN;

    double M = bh->M;
    double a = prograde ? fabs(bh->a) : -fabs(bh->a);
    double Omega = 1.0 / (a + r * sqrt(r) / sqrt(M));

    double term = 1.0 - 6.0*M/r
                  + 8.0*a*sqrt(M)/(r*sqrt(r))
                  - 3.0*a*a/(r*r);

    if (term < 0.0) return 0.0;  /* Orbit is unstable */

    return Omega * sqrt(term);
}

double kerr_epicyclic_frequency_theta(const KerrBlackHole *bh, double r,
                                       int prograde)
{
    if (r <= 0.0) return NAN;

    double M = bh->M;
    double a = prograde ? fabs(bh->a) : -fabs(bh->a);
    double Omega = 1.0 / (a + r * sqrt(r) / sqrt(M));

    double term = 1.0 - 4.0*a*sqrt(M)/(r*sqrt(r))
                  + 3.0*a*a/(r*r);

    if (term < 0.0) return 0.0;

    return Omega * sqrt(term);
}

/* ====================================================================
 * Keplerian orbital frequency
 *
 * Omega_K = 1 / (a + r^(3/2) / sqrt(M))
 *
 * This is the angular velocity of a circular equatorial orbit
 * as measured by a distant observer.
 * ==================================================================== */

double kerr_keplerian_frequency(const KerrBlackHole *bh, double r,
                                int prograde)
{
    if (r <= 0.0) return NAN;
    double a = prograde ? fabs(bh->a) : -fabs(bh->a);
    return 1.0 / (a + r * sqrt(r) / sqrt(bh->M));
}

double kerr_periastron_precession(const KerrBlackHole *bh, double r,
                                  int prograde)
{
    double Omega = kerr_keplerian_frequency(bh, r, prograde);
    double kappa_r = kerr_epicyclic_frequency_r(bh, r, prograde);
    if (isnan(kappa_r)) return NAN;
    return Omega - kappa_r;
}

double kerr_nodal_precession(const KerrBlackHole *bh, double r,
                             int prograde)
{
    double Omega = kerr_keplerian_frequency(bh, r, prograde);
    double Omega_th = kerr_epicyclic_frequency_theta(bh, r, prograde);
    if (isnan(Omega_th)) return NAN;
    return Omega - Omega_th;
}
