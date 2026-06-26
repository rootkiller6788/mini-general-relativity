/**
 * schwarzschild_geodesics.c — Geodesic equations in Schwarzschild spacetime
 *
 * Reference: Wald (1984) Ch.6; Carroll (2004) Ch.5; MTW (1973) Ch.25
 * MIT 8.962 — General Relativity
 *
 * Knowledge coverage:
 *   L2 — Core Concepts: constants of motion, effective potential
 *   L4 — Fundamental Laws: geodesic equation
 *   L6 — Canonical Systems: circular orbits, ISCO, photon sphere, bound orbits
 */

#include "schwarzschild_geodesics.h"
#include "schwarzschild_defs.h"
#include <math.h>
#include <float.h>

/* ==========================================================================
 * Effective potential
 * ========================================================================== */

double schwarzschild_V_eff_massive(double r, double rs, double L_sq) {
    if (r <= 0.0) return INFINITY;
    double f = 1.0 - rs / r;
    return f * (1.0 + L_sq / (r * r));
}

double schwarzschild_V_eff_massless(double r, double rs, double L_sq) {
    if (r <= 0.0) return INFINITY;
    double f = 1.0 - rs / r;
    return f * (L_sq / (r * r));
}

double schwarzschild_dV_eff_dr_massive(double r, double rs, double L_sq) {
    if (r <= 0.0) return 0.0;
    double r2 = r * r;
    double r3 = r2 * r;
    double f = 1.0 - rs / r;
    /*
     * V_eff = (1 - rs/r) * (1 + L^2/r^2)
     * dV/dr = (rs/r^2)*(1 + L^2/r^2) - (1 - rs/r)*(2*L^2/r^3)
     *       = rs/r^2 + rs*L^2/r^4 - 2*L^2/r^3 + 2*rs*L^2/r^4
     *       = rs/r^2 - 2*L^2/r^3 + 3*rs*L^2/r^4
     */
    return rs / r2 - 2.0 * L_sq / r3 + 3.0 * rs * L_sq / (r3 * r);
}

double schwarzschild_dV_eff_dr_massless(double r, double rs, double L_sq) {
    if (r <= 0.0) return 0.0;
    double r2 = r * r;
    double r3 = r2 * r;
    /*
     * V_eff = (1 - rs/r) * L^2/r^2
     * dV/dr = (rs*L^2/r^4) - 2*L^2*(1-rs/r)/r^3
     *       = L^2*rs/r^4 - 2*L^2/r^3 + 2*rs*L^2/r^4
     *       = -2*L^2/r^3 + 3*rs*L^2/r^4
     */
    return -2.0 * L_sq / r3 + 3.0 * rs * L_sq / (r3 * r);
}

double schwarzschild_d2V_eff_dr2_massive(double r, double rs, double L_sq) {
    if (r <= 0.0) return 0.0;
    double r3 = r * r * r;
    double r4 = r3 * r;
    double r5 = r4 * r;
    /*
     * d^2V/dr^2 = -2*rs/r^3 + 6*L^2/r^4 - 12*rs*L^2/r^5
     */
    return -2.0 * rs / r3 + 6.0 * L_sq / r4 - 12.0 * rs * L_sq / r5;
}

void schwarzschild_effective_potential_massive(double r, double rs, double L_sq,
                                                SchwarzschildEffectivePot *pot) {
    if (!pot) return;
    pot->r = r;
    pot->V_eff = schwarzschild_V_eff_massive(r, rs, L_sq);
    pot->dV_dr = schwarzschild_dV_eff_dr_massive(r, rs, L_sq);
    pot->d2V_dr2 = schwarzschild_d2V_eff_dr2_massive(r, rs, L_sq);
}

void schwarzschild_effective_potential_massless(double r, double rs, double L_sq,
                                                 SchwarzschildEffectivePot *pot) {
    if (!pot) return;
    pot->r = r;
    pot->V_eff = schwarzschild_V_eff_massless(r, rs, L_sq);
    pot->dV_dr = schwarzschild_dV_eff_dr_massless(r, rs, L_sq);
    pot->d2V_dr2 = 0.0;  /* not implemented for massless yet */
}

/* ==========================================================================
 * Circular orbits
 * ========================================================================== */

double schwarzschild_L_for_circular_orbit(double r, double rs) {
    /*
     * From dV_eff/dr = 0:
     *   rs/r^2 - 2*L^2/r^3 + 3*rs*L^2/r^4 = 0
     * => L^2 * (3*rs - 2*r) / r^4 = -rs/r^2
     * => L^2 = rs*r^2 / (r - 3*rs/2) ... wait let me do this properly.
     *
     * dV/dr = rs/r^2 - 2*L^2/r^3 + 3*rs*L^2/r^4 = 0
     * rs*r^2 - 2*L^2*r + 3*rs*L^2 = 0
     * L^2*(3*rs - 2*r) = -rs*r^2
     * L^2 = rs*r^2 / (2*r - 3*rs)
     *
     * Actually, let me re-derive more carefully.
     * V_eff = (1 - rs/r)(1 + L^2/r^2) = (1 - 2m/r)(1 + L^2/r^2) in geometric units
     *
     * dV_dr = (2m/r^2)(1 + L^2/r^2) - (1 - 2m/r)(2L^2/r^3)
     *       = 2m/r^2 + 2mL^2/r^4 - 2L^2/r^3 + 4mL^2/r^4
     *       = 2m/r^2 - 2L^2/r^3 + 6mL^2/r^4
     *
     * Set dV_dr = 0:
     * 2m/r^2 - 2L^2/r^3 + 6mL^2/r^4 = 0
     * Multiply by r^4: 2m*r^2 - 2L^2*r + 6m*L^2 = 0
     * Divide by 2: m*r^2 - L^2*r + 3m*L^2 = 0
     * L^2*(r - 3m) = m*r^2
     * L^2 = m*r^2 / (r - 3m)
     *
     * In terms of rs = 2m: L^2 = (rs/2)*r^2 / (r - 3*rs/2) = rs*r^2 / (2*r - 3*rs)
     */
    double m = rs / 2.0;
    if (r <= 3.0 * m) return INFINITY;
    return m * r * r / (r - 3.0 * m);
}

double schwarzschild_E_for_circular_orbit(double r, double rs) {
    /*
     * E^2 = V_eff(r_circ) = (1 - 2m/r)(1 + L^2/r^2)
     * With L^2 = m*r^2/(r - 3m):
     * E^2 = (r - 2m)^2 / (r*(r - 3m))
     */
    double m = rs / 2.0;
    if (r <= 3.0 * m) return 0.0;
    double num = (r - 2.0 * m) * (r - 2.0 * m);
    double den = r * (r - 3.0 * m);
    return sqrt(num / den);
}

int schwarzschild_is_stable_circular_orbit(double r, double rs) {
    /*
     * Stable if d^2V/dr^2 > 0 at the circular orbit.
     * d^2V/dr^2 = -2*rs/r^3 + 6*L^2/r^4 - 12*rs*L^2/r^5
     * With L^2 from the circular condition, this simplifies to:
     *   d^2V/dr^2 = 2*m*(r - 6*m) / (r^3*(r - 3*m))
     * So stable for r > 6*m.
     */
    double m = rs / 2.0;
    if (r <= 3.0 * m) return 0;
    return (r > 6.0 * m) ? 1 : 0;
}

double schwarzschild_photon_sphere_radius(double rs) {
    double m = rs / 2.0;
    return 3.0 * m;
}

double schwarzschild_isco_radius(double rs) {
    double m = rs / 2.0;
    return 6.0 * m;
}

double schwarzschild_marginally_bound_radius(double rs) {
    double m = rs / 2.0;
    return 4.0 * m;
}

double schwarzschild_L_isco(double rs) {
    /* L_isco = sqrt(12) * m */
    double m = rs / 2.0;
    return sqrt(12.0) * m;
}

double schwarzschild_E_isco(void) {
    /* E_isco = sqrt(8/9) */
    return sqrt(8.0 / 9.0);
}

double schwarzschild_Omega_isco(double mass) {
    /*
     * Omega_isco = sqrt(G*M / r_isco^3) in Newtonian,
     * but the relativistic value is:
     * Omega_isco = c^3 / (G*M * 6*sqrt(6))
     *
     * Actually: Omega = dphi/dt for circular orbit.
     * From the geodesic equation:
     * Omega = sqrt(m/r^3) (Kepler's 3rd in geometric units)
     * For r = 6m: Omega = sqrt(m/(216*m^3)) = 1/(6*sqrt(6)*m) in geometric.
     * In SI: Omega = c^3/(6*sqrt(6)*G*M)
     */
    if (mass <= 0.0) return 0.0;
    double m_geo = geometric_mass(mass);
    return SCHW_C / (6.0 * sqrt(6.0) * m_geo);
}

/* ==========================================================================
 * Orbital periods
 * ========================================================================== */

double schwarzschild_orbital_period_keplerian(double r, double mass) {
    if (mass <= 0.0 || r <= 0.0) return 0.0;
    return 2.0 * M_PI * sqrt(r * r * r / (SCHW_G_N * mass));
}

double schwarzschild_orbital_period_relativistic(double r, double mass) {
    if (mass <= 0.0 || r <= 0.0) return 0.0;
    double m = geometric_mass(mass);
    double T_kep = schwarzschild_orbital_period_keplerian(r, mass);
    double f1 = 1.0 - 3.0 * m / r;
    double f2 = 1.0 - 2.0 * m / r;
    if (f1 <= 0.0 || f2 <= 0.0) return 0.0;
    return T_kep * sqrt(f1 / f2);
}

/* ==========================================================================
 * Radial acceleration
 * ========================================================================== */

double schwarzschild_radial_accel_massive(double r, double rs, double L_sq) {
    /*
     * From the geodesic equation:
     * d^2r/dtau^2 = -m/r^2 + L^2/r^3 - 3*m*L^2/r^4
     * (Geometric units, m = rs/2)
     *
     * In SI: d^2r/dtau^2 = -GM/r^2 + L^2/r^3 - 3*(GM/c^2)*L^2/r^4
     * where L is angular momentum per unit mass.
     */
    if (r <= 0.0) return 0.0;
    double m = rs / 2.0;
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r3 * r;
    return -m / r2 + L_sq / r3 - 3.0 * m * L_sq / r4;
}

double schwarzschild_radial_accel_massless(double r, double rs, double L_sq) {
    /*
     * For null geodesics:
     * d^2r/dlambda^2 = L^2/r^3 - 3*m*L^2/r^4
     * (No -m/r^2 term since massless particles have no rest mass term)
     */
    if (r <= 0.0) return 0.0;
    double m = rs / 2.0;
    double r3 = r * r * r;
    double r4 = r3 * r;
    return L_sq / r3 - 3.0 * m * L_sq / r4;
}

/* ==========================================================================
 * Orbit classification
 * ========================================================================== */

SchwarzschildOrbitType schwarzschild_classify_orbit(double r, double E_sq,
                                                      double L_sq, double rs) {
    double V = schwarzschild_V_eff_massive(r, rs, L_sq);
    double dV = schwarzschild_dV_eff_dr_massive(r, rs, L_sq);
    double d2V = schwarzschild_d2V_eff_dr2_massive(r, rs, L_sq);

    if (r <= rs) return ORBIT_PLUNGE;
    if (E_sq > V) {
        /* Particle has enough energy to be at this radius */
        if (fabs(dV) < 1e-12) {
            /* At an extremum of V_eff */
            if (d2V > 0.0) return ORBIT_CIRCULAR;
            else return ORBIT_UNSTABLE_CIRCULAR;
        }
        /* Check if this radius is between turning points */
        if (E_sq < 1.0) return ORBIT_BOUND;
        else if (E_sq > 1.0) return ORBIT_UNBOUND;
        else return ORBIT_MARGINAL;
    }
    return ORBIT_PLUNGE;
}

int schwarzschild_radial_turning_points(double E_sq, double L_sq, double rs,
                                         double roots[4]) {
    /*
     * Solve E^2 = (1 - rs/r)(1 + L^2/r^2) for r.
     * => E^2 = (1 - rs/r) + L^2/r^2 - rs*L^2/r^3
     * => 0 = (E^2 - 1)*r^3 + rs*r^2 - L^2*r + rs*L^2
     *
     * This is a cubic in r. For bound orbits (E^2 < 1):
     * there are typically 2 positive roots (periastron, apastron)
     * and 1 negative root (unphysical).
     *
     * Use Cardano's formula.
     */
    double a = E_sq - 1.0;
    double b = rs;
    double c = -L_sq;
    double d = rs * L_sq;

    /*
     * Normalize: r^3 + (b/a)*r^2 + (c/a)*r + (d/a) = 0
     */
    if (fabs(a) < 1e-15) {
        /* E^2 = 1: marginal case, quadratic */
        if (fabs(b) < 1e-15) return 0;
        double disc = c * c - 4.0 * b * d;
        if (disc < 0) return 0;
        roots[0] = (-c + sqrt(disc)) / (2.0 * b);
        roots[1] = (-c - sqrt(disc)) / (2.0 * b);
        roots[2] = 0.0; roots[3] = 0.0;
        return (roots[0] > rs) + (roots[1] > rs);
    }

    double bn = b / a;
    double cn = c / a;
    double dn = d / a;

    /* Cardano: let r = y - bn/3 */
    double p = cn - bn * bn / 3.0;
    double q = 2.0 * bn * bn * bn / 27.0 - bn * cn / 3.0 + dn;

    double disc = q * q / 4.0 + p * p * p / 27.0;

    double y1, y2, y3;
    if (disc > 0) {
        double u = cbrt(-q / 2.0 + sqrt(disc));
        double v = cbrt(-q / 2.0 - sqrt(disc));
        y1 = u + v;
        y2 = -(u + v) / 2.0;
        y3 = y2;
    } else if (fabs(disc) < 1e-15) {
        y1 = 2.0 * cbrt(-q / 2.0);
        y2 = -cbrt(-q / 2.0);
        y3 = y2;
    } else {
        double phi = acos(-q / (2.0 * sqrt(-p * p * p / 27.0)));
        double r_mag = 2.0 * sqrt(-p / 3.0);
        y1 = r_mag * cos(phi / 3.0);
        y2 = r_mag * cos((phi + 2.0 * M_PI) / 3.0);
        y3 = r_mag * cos((phi + 4.0 * M_PI) / 3.0);
    }

    double shift = bn / 3.0;
    double r1 = y1 - shift;
    double r2 = y2 - shift;
    double r3 = y3 - shift;

    /* Collect positive roots > rs */
    int count = 0;
    if (r1 > rs) roots[count++] = r1;
    if (r2 > rs && fabs(r2 - r1) > 1e-10) roots[count++] = r2;
    if (r3 > rs && fabs(r3 - r1) > 1e-10 && fabs(r3 - r2) > 1e-10)
        roots[count++] = r3;
    roots[count] = 0.0;

    /* Sort ascending */
    for (int i = 0; i < count - 1; i++)
        for (int j = i + 1; j < count; j++)
            if (roots[i] > roots[j]) {
                double tmp = roots[i]; roots[i] = roots[j]; roots[j] = tmp;
            }

    return count;
}

/* ==========================================================================
 * Angular equations of motion
 * ========================================================================== */

double schwarzschild_dphi_dtau_massive(double r, double L) {
    if (r <= 0.0) return 0.0;
    return L / (r * r);
}

double schwarzschild_dtheta_dtau(double r, double theta, double Q, double L_sq) {
    if (r <= 0.0) return 0.0;
    double tan_theta = tan(theta);
    double term = (tan_theta != 0.0) ? L_sq / (tan_theta * tan_theta) : 0.0;
    double arg = Q - term;
    if (arg < 0.0) return 0.0;
    return sqrt(arg) / (r * r);
}

double schwarzschild_dt_dtau_massive(double E, double r, double rs) {
    if (r <= rs) return INFINITY;
    return E / (1.0 - rs / r);
}

/* ==========================================================================
 * Proper time and coordinate time for radial infall
 * ========================================================================== */

double schwarzschild_proper_time_radial_infall(double r_start, double r_end,
                                                double rs) {
    /*
     * For radial infall from rest at infinity (E = 1, L = 0):
     *
     * (dr/dtau)^2 = E^2 - V_eff = 1 - (1 - rs/r) = rs/r
     *
     * So dr/dtau = -sqrt(rs/r) (infalling)
     *
     * tau = integral_{r_start}^{r_end} dr / sqrt(rs/r)
     *     = integral_{r_start}^{r_end} sqrt(r/rs) dr
     *     = (2/3) * (r_start^{3/2} - r_end^{3/2}) / sqrt(rs)
     */
    if (r_start <= rs || r_end <= rs || rs <= 0.0) return 0.0;
    if (r_start < r_end) { double t = r_start; r_start = r_end; r_end = t; }
    double sqrt_rs = sqrt(rs);
    return (2.0 / 3.0) * (pow(r_start, 1.5) - pow(r_end, 1.5)) / sqrt_rs;
}

double schwarzschild_coordinate_time_radial_infall(double r_start, double r_end,
                                                    double rs) {
    /*
     * Coordinate time for radial infall (E = 1):
     *
     * dt/dr = (dt/dtau) / (dr/dtau) = [1/(1-rs/r)] / [-sqrt(rs/r)]
     *       = -sqrt(r^3/rs) / (r - rs)
     *
     * t = integral dr * sqrt(r^3/rs) / (r - rs)
     *
     * This integral diverges as r -> rs (infinite coordinate time to cross horizon).
     *
     * The analytic antiderivative is:
     * t = -(2/3)*sqrt(r^3/rs) - 2*sqrt(r*rs) + rs*ln|(sqrt(r)+sqrt(rs))/(sqrt(r)-sqrt(rs))|
     *
     * evaluated at r_start and r_end.
     */
    if (r_start <= rs || r_end <= rs || rs <= 0.0) return 0.0;
    if (r_start < r_end) { double tmp = r_start; r_start = r_end; r_end = tmp; }
    double t_start = -(2.0/3.0)*sqrt(r_start*r_start*r_start/rs)
                     - 2.0*sqrt(r_start*rs)
                     + rs*log((sqrt(r_start)+sqrt(rs))/(sqrt(r_start)-sqrt(rs)));
    double t_end = -(2.0/3.0)*sqrt(r_end*r_end*r_end/rs)
                   - 2.0*sqrt(r_end*rs)
                   + rs*log((sqrt(r_end)+sqrt(rs))/(sqrt(r_end)-sqrt(rs)));
    return t_end - t_start;
}

/* ==========================================================================
 * 8D geodesic system (Hamiltonian form)
 * ========================================================================== */

void schwarzschild_geodesic_derivs(double lambda, const double y[8],
                                    double dydlambda[8], double rs) {
    (void)lambda;
    double r = y[1];
    double theta = y[2];
    double pt = y[4];
    double pr = y[5];
    double ptheta = y[6];
    double pphi = y[7];

    if (r <= 0.0) {
        for (int i = 0; i < 8; i++) dydlambda[i] = 0.0;
        return;
    }

    double st = sin(theta);
    double f = 1.0 - rs / r;

    /*
     * Hamiltonian: H = (1/2) g^{mu nu} p_mu p_nu
     * For timelike: H = -1/2
     * For null: H = 0
     *
     * dx^mu/dlambda = dH/dp_mu = g^{mu nu} p_nu
     * dp_mu/dlambda = -dH/dx^mu = -(1/2) p_alpha p_beta d(g^{alpha beta})/dx^mu
     */

    /* dt/dlambda = -pt / f */
    dydlambda[0] = (f != 0.0) ? -pt / f : 0.0;

    /* dr/dlambda = f * pr */
    dydlambda[1] = f * pr;

    /* dtheta/dlambda = ptheta / r^2 */
    dydlambda[2] = ptheta / (r * r);

    /* dphi/dlambda = pphi / (r^2 * sin^2(theta)) */
    if (st != 0.0)
        dydlambda[3] = pphi / (r * r * st * st);
    else
        dydlambda[3] = 0.0;

    /* dp_t/dlambda = 0 (stationary spacetime, t is cyclic) */
    dydlambda[4] = 0.0;

    /*
     * dp_r/dlambda = -(1/2) * [pt^2 * dg^{tt}/dr + pr^2 * dg^{rr}/dr
     *                          + ptheta^2 * dg^{theta theta}/dr
     *                          + pphi^2 * dg^{phi phi}/dr]
     *
     * g^{tt} = -1/f  =>  dg^{tt}/dr = -rs/(r^2*f^2)
     * g^{rr} = f     =>  dg^{rr}/dr = rs/r^2
     * g^{theta theta} = 1/r^2  =>  d/dr = -2/r^3
     * g^{phi phi} = 1/(r^2*st^2)  =>  d/dr = -2/(r^3*st^2)
     */
    double r2 = r * r;
    double r3 = r2 * r;
    dydlambda[5] = 0.0;
    if (f != 0.0) {
        dydlambda[5] -= 0.5 * pt * pt * (-rs / (r2 * f * f));
    }
    dydlambda[5] -= 0.5 * pr * pr * (rs / r2);
    dydlambda[5] -= 0.5 * ptheta * ptheta * (-2.0 / r3);
    if (st != 0.0)
        dydlambda[5] -= 0.5 * pphi * pphi * (-2.0 / (r3 * st * st));

    /*
     * dp_theta/dlambda = -(1/2) * [pphi^2 * dg^{phi phi}/dtheta]
     * dg^{phi phi}/dtheta = -2*cos(theta)/(r^2*sin^3(theta))
     *
     * = -(1/2) * pphi^2 * (-2*ct) / (r^2*st^3)
     * = pphi^2 * ct / (r^2*st^3)
     */
    if (st != 0.0)
        dydlambda[6] = pphi * pphi * cos(theta) / (r2 * st * st * st);
    else
        dydlambda[6] = 0.0;

    /* dp_phi/dlambda = 0 (axisymmetry) */
    dydlambda[7] = 0.0;
}