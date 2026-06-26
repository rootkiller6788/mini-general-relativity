/**
 * schwarzschild_observables.c — Observable effects of Schwarzschild spacetime
 *
 * Reference: Wald (1984) Ch.6; Carroll (2004) Ch.5; Will (2018)
 * MIT 8.962 — General Relativity
 *
 * Knowledge coverage:
 *   L2 — Core Concepts: gravitational time dilation, redshift
 *   L6 — Canonical Systems: light deflection, perihelion precession, Shapiro delay
 *   L7 — Applications: GPS corrections, black hole shadow, accretion disk
 */

#include "schwarzschild_observables.h"
#include "schwarzschild_defs.h"
#include "schwarzschild_geodesics.h"
#include <math.h>
#include <float.h>

/* ==========================================================================
 * Gravitational Time Dilation
 * ========================================================================== */

double schwarzschild_time_dilation_factor(double r, double rs) {
    if (r <= rs) return 0.0;
    return sqrt(1.0 - rs / r);
}

double schwarzschild_time_dilation_ratio(double r1, double r2, double rs) {
    double f1 = schwarzschild_time_dilation_factor(r1, rs);
    double f2 = schwarzschild_time_dilation_factor(r2, rs);
    if (f2 <= 0.0) return INFINITY;
    return f1 / f2;
}

double schwarzschild_gps_time_dilation(double r_satellite, double r_earth,
                                        double mass_earth, double v_satellite) {
    /*
     * Combined GR + SR time dilation for GPS satellites.
     *
     * GR contribution (clocks run FASTER in weaker gravity):
     *   (Delta f / f)_GR = GM/(c^2) * (1/r_earth - 1/r_satellite) > 0
     *   => Satellite clock faster by ~45.7 us/day
     *
     * SR contribution (clocks run SLOWER due to velocity):
     *   (Delta f / f)_SR = -v^2/(2*c^2) < 0
     *   => Satellite clock slower by ~7.2 us/day
     *
     * Total: (Delta f / f)_total = GM/c^2 * (1/r_earth - 1/r_satellite) - v^2/(2*c^2)
     *
     * For GPS: net ~ +38.5 us/day
     *
     * Reference: Ashby (2003), Living Reviews in Relativity
     */
    if (r_satellite <= 0.0 || r_earth <= 0.0 || mass_earth <= 0.0) return 0.0;
    double gr = SCHW_G_N * mass_earth / (SCHW_C * SCHW_C)
                * (1.0 / r_earth - 1.0 / r_satellite);
    double sr = -0.5 * v_satellite * v_satellite / (SCHW_C * SCHW_C);
    return gr + sr;
}

/* ==========================================================================
 * Gravitational Redshift
 * ========================================================================== */

double schwarzschild_redshift_at_infinity(double r_emitter, double rs) {
    if (r_emitter <= rs) return INFINITY;
    return 1.0 / sqrt(1.0 - rs / r_emitter) - 1.0;
}

double schwarzschild_redshift_finite(double r_emitter, double r_observer, double rs) {
    if (r_emitter <= rs || r_observer <= rs) return 0.0;
    double num = 1.0 - rs / r_observer;
    double den = 1.0 - rs / r_emitter;
    if (num < 0.0 || den <= 0.0) return 0.0;
    return sqrt(num / den) - 1.0;
}

double schwarzschild_redshift_isco_at_infinity(double rs) {
    double r_isco = schwarzschild_isco_radius(rs);
    return schwarzschild_redshift_at_infinity(r_isco, rs);
}

double schwarzschild_doppler_circular_orbit(double r, double rs,
                                              double impact_param, double incl_angle) {
    /*
     * g = E_obs / E_em = 1 / [u^t * (1 - Omega * b * sin(incl))]
     *
     * For circular orbit at radius r:
     * u^t = 1 / sqrt(1 - 3*m/r)
     * Omega = sqrt(m/r^3)
     * b = impact parameter (at infinity)
     */
    double m = rs / 2.0;
    if (r <= 3.0 * m) return 0.0;
    double ut = 1.0 / sqrt(1.0 - 3.0 * m / r);
    double Omega = sqrt(m / (r * r * r));
    double denom = ut * (1.0 - Omega * impact_param * sin(incl_angle));
    if (denom <= 0.0) return 0.0;
    return 1.0 / denom;
}

/* ==========================================================================
 * Light Deflection
 * ========================================================================== */

double schwarzschild_light_deflection_angle(double impact_param, double mass) {
    /*
     * alpha = 4*G*M / (b * c^2) = 2*rs / b
     *
     * For the Sun: M_sun = 1.98847e30 kg, R_sun = 6.96e8 m
     * b = R_sun => alpha = 4*6.674e-11*1.988e30/(6.96e8*9e16) = 8.48e-6 rad = 1.75 arcsec
     */
    if (impact_param <= 0.0 || mass <= 0.0) return 0.0;
    double rs = schwarzschild_radius(mass);
    return 2.0 * rs / impact_param;
}

double schwarzschild_light_deflection_exact(double impact_param, double rs,
                                             double tolerance) {
    /*
     * Exact deflection via numerical integration:
     * alpha = 2 * integral_{r0}^{infty} dr / (r * sqrt(r^2/impact_param^2 - (1 - rs/r))) - pi
     *
     * where r0 is the turning point: r0^3 - impact_param^2*(r0 - rs) = 0
     *
     * Simplified: use u = impact_param/r substitution:
     * du/dphi = sqrt(1 - impact_param^2 * u^2 * (1 - rs*u))
     *
     * Actually: for a null geodesic, the orbit equation is:
     * (du/dphi)^2 + u^2 = (rs/b^2) * u^3 + u^2 / b^2 ... hmm.
     *
     * Let me use the standard approach:
     * d^2u/dphi^2 + u = (3/2)*rs*u^2
     *
     * For small deflection, integrate perturbatively.
     *
     * For the exact formula, use elliptic integrals. The proper expression is:
     * alpha = 4/sqrt(r0*(r0-rs)*(3*rs/r0 - rs^2/r0^2)) * K(k) - pi
     *
     * where K is the complete elliptic integral of the first kind.
     *
     * We'll use numerical integration with adaptive Simpson for robustness.
     */
    if (impact_param <= 0.0 || rs <= 0.0) return 0.0;

    /* Find turning point r0: r0^3 - impact_param^2*(r0 - rs) = 0 */
    double r0 = impact_param;
    if (impact_param < 3.0 * sqrt(3.0) * rs / 2.0) {
        return INFINITY;  /* photon captured */
    }
    /* Newton-Raphson for r0 */
    for (int i = 0; i < 50; i++) {
        double f = r0 * r0 * r0 - impact_param * impact_param * (r0 - rs);
        double df = 3.0 * r0 * r0 - impact_param * impact_param;
        if (fabs(df) < 1e-15) break;
        double dr = -f / df;
        r0 += dr;
        if (fabs(dr) < tolerance * r0) break;
    }

    /* Integrate: alpha = 2 * integral_{r0}^{R_max} dr/(r*sqrt(r^2/b^2 - 1 + rs/r)) - pi */
    double R_max = 1000.0 * fmax(impact_param, rs);
    double integral = 0.0;

    /* Use u = 1/r substitution: integral_{u_min}^{0} du / sqrt(1/b^2 - u^2 + rs*u^3) */
    /* But this has a root at u=0. Better to integrate in r. */
    int n = 2000;
    double log_r0 = log(r0);
    double log_Rmax = log(R_max);
    double dlogr = (log_Rmax - log_r0) / n;

    for (int i = 0; i < n; i++) {
        double log_r1 = log_r0 + i * dlogr;
        double log_r2 = log_r0 + (i + 1) * dlogr;
        double r_a = exp(log_r1);
        double r_b = exp(log_r2);
        double r_mid = exp((log_r1 + log_r2) / 2.0);

        double integrand(double rr) {
            double arg = rr * rr / (impact_param * impact_param) - (1.0 - rs / rr);
            if (arg <= 0.0) return 0.0;
            return 1.0 / (rr * sqrt(arg));
        }

        double fa = integrand(r_a);
        double fb = integrand(r_b);
        double fm = integrand(r_mid);
        integral += (r_b - r_a) * (fa + 4.0 * fm + fb) / 6.0;
    }

    return 2.0 * integral - M_PI;
}

double schwarzschild_einstein_ring_radius(double mass, double D_l,
                                            double D_s, double D_ls) {
    /*
     * theta_E = sqrt(4*G*M/c^2 * D_ls/(D_l*D_s))
     *
     * For a galaxy lens M ~ 1e12 Msun at D_l ~ 1 Gpc, source at D_s ~ 2 Gpc:
     * theta_E ~ 1 arcsecond
     */
    if (D_l <= 0.0 || D_s <= 0.0 || D_ls <= 0.0 || mass <= 0.0) return 0.0;
    double rs = schwarzschild_radius(mass);
    return sqrt(2.0 * rs * D_ls / (D_l * D_s));
}

double schwarzschild_lens_magnification(double u) {
    /*
     * mu = (u^2 + 2) / (u * sqrt(u^2 + 4))
     *
     * u = beta / theta_E is the angular separation in Einstein radius units.
     * For u << 1 (near-perfect alignment): mu ~ 1/u (diverges)
     * For u >> 1 (large separation): mu ~ 1 + 2/u^4
     */
    if (fabs(u) < 1e-15) return INFINITY;
    return (u * u + 2.0) / (u * sqrt(u * u + 4.0));
}

int schwarzschild_lens_image_positions(double beta, double theta_E,
                                         double *theta_plus, double *theta_minus) {
    /*
     * theta_{1,2} = (beta +/- sqrt(beta^2 + 4*theta_E^2)) / 2
     */
    if (!theta_plus || !theta_minus || theta_E <= 0.0) return 0;
    double disc = beta * beta + 4.0 * theta_E * theta_E;
    double sqrt_disc = sqrt(disc);
    *theta_plus = 0.5 * (beta + sqrt_disc);
    *theta_minus = 0.5 * (beta - sqrt_disc);
    return 2;
}

/* ==========================================================================
 * Shapiro Time Delay
 * ========================================================================== */

double schwarzschild_shapiro_log_term(double r_emitter, double r_reflector,
                                        double impact_param, double geometric_mass) {
    if (r_emitter <= 0.0 || r_reflector <= 0.0 || impact_param <= 0.0) return 0.0;
    /*
     * Delta_t = 2*m/c * ln(4*r_e*r_r / b^2)
     * m = geometric_mass = GM/c^2
     */
    double arg = 4.0 * r_emitter * r_reflector / (impact_param * impact_param);
    if (arg <= 0.0) return 0.0;
    return 2.0 * geometric_mass / SCHW_C * log(arg);
}

double schwarzschild_shapiro_time_delay(double r_emitter, double r_reflector,
                                          double impact_param, double mass) {
    double m_geo = geometric_mass(mass);
    return schwarzschild_shapiro_log_term(r_emitter, r_reflector,
                                            impact_param, m_geo);
}

double schwarzschild_shapiro_earth_venus(double earth_orbit_radius,
                                           double venus_orbit_radius,
                                           double solar_mass) {
    /*
     * For Earth-Venus-Sun superior conjunction:
     * r_earth ~ 1 AU, r_venus ~ 0.723 AU
     * impact_param ~ R_sun ~ 6.96e8 m
     *
     * Expected: Delta_t ~ 240 microseconds (round trip)
     *
     * The delay is maximum at superior conjunction when the radar
     * signal passes closest to the Sun.
     */
    double r_sun = 6.957e8; /* solar radius [m] */
    double b = r_sun;
    double m_geo = geometric_mass(solar_mass);
    /* Round trip: 2x the one-way delay */
    return 2.0 * schwarzschild_shapiro_log_term(earth_orbit_radius,
                                                  venus_orbit_radius, b, m_geo);
}

/* ==========================================================================
 * Perihelion Precession
 * ========================================================================== */

double schwarzschild_perihelion_precession_per_orbit(double semi_major_axis,
                                                       double eccentricity,
                                                       double mass) {
    /*
     * Delta_phi = 6*pi*G*M / (a * c^2 * (1 - e^2))
     *
     * For Mercury: a = 5.79e10 m, e = 0.2056, M = M_sun
     * Delta_phi = 6*pi*6.6743e-11*1.9885e30 / (5.79e10 * 9e16 * (1-0.2056^2))
     *            = 5.01e-7 rad/orbit = 0.1035 arcsec/orbit
     * Period = 0.2409 yr => 0.1035/0.2409 = 0.429 arcsec/yr = 42.9 arcsec/century
     */
    if (semi_major_axis <= 0.0 || eccentricity < 0.0 || eccentricity >= 1.0)
        return 0.0;
    double num = 6.0 * M_PI * SCHW_G_N * mass;
    double den = semi_major_axis * SCHW_C * SCHW_C * (1.0 - eccentricity * eccentricity);
    if (den <= 0.0) return 0.0;
    return num / den;
}

double schwarzschild_perihelion_precession_rate(double semi_major_axis,
                                                  double eccentricity,
                                                  double mass) {
    double per_orbit = schwarzschild_perihelion_precession_per_orbit(
        semi_major_axis, eccentricity, mass);
    double T = schwarzschild_orbital_period_keplerian(semi_major_axis, mass);
    if (T <= 0.0) return 0.0;
    return per_orbit / T;
}

double schwarzschild_hulse_taylor_precession(double m1, double m2,
                                               double semi_major_axis,
                                               double eccentricity) {
    /*
     * For a binary system, the periastron advance is:
     * domega/dt = 3 * (G*(m1+m2)/c^3)^(2/3) * (P_b/(2*pi))^(-5/3) * 1/(1-e^2)
     *
     * where P_b is the orbital period.
     *
     * Simplified: Delta_omega = 6*pi*G*(m1+m2) / (a*c^2*(1-e^2)) per orbit
     * Rate: domega/dt = Delta_omega / P_b
     *
     * For PSR B1913+16: m1 ~ 1.44 Msun, m2 ~ 1.39 Msun,
     * a ~ 1.95e9 m, e ~ 0.617
     * => domega/dt ~ 4.226 degrees/year
     */
    double M_total = m1 + m2;
    double per_orbit = schwarzschild_perihelion_precession_per_orbit(
        semi_major_axis, eccentricity, M_total);
    double P_b = 2.0 * M_PI * sqrt(semi_major_axis * semi_major_axis * semi_major_axis
                                    / (SCHW_G_N * M_total));
    if (P_b <= 0.0) return 0.0;
    return per_orbit / P_b;
}

/* ==========================================================================
 * Black Hole Shadow
 * ========================================================================== */

double schwarzschild_shadow_angular_radius(double mass, double distance) {
    /*
     * theta_shadow = sqrt(27) * m / D
     *
     * For M87*: M ~ 6.5e9 Msun, D ~ 16.8 Mpc
     * m = GM/c^2 ~ 6.5e9 * 1.477 km ~ 9.6e12 m
     * D ~ 16.8 * 3.086e22 ~ 5.18e23 m
     * theta ~ sqrt(27)*9.6e12/5.18e23 ~ 9.66e-11 rad ~ 20 uas
     * (plus factor for photon ring vs shadow boundary, actual ~42 uas)
     */
    if (distance <= 0.0 || mass <= 0.0) return 0.0;
    double m = geometric_mass(mass);
    return sqrt(27.0) * m / distance;
}

double schwarzschild_photon_capture_cross_section(double mass) {
    double m = geometric_mass(mass);
    return 27.0 * M_PI * m * m;
}

double schwarzschild_critical_impact_parameter(double rs) {
    return 3.0 * sqrt(3.0) * rs / 2.0;
}

/* ==========================================================================
 * Accretion Disk Observables
 * ========================================================================== */

double schwarzschild_accretion_efficiency(void) {
    return 1.0 - sqrt(8.0 / 9.0);
}

double schwarzschild_disk_inner_edge(double rs) {
    return schwarzschild_isco_radius(rs);
}

double schwarzschild_disk_flux(double r, double mass, double mdot) {
    /*
     * F(r) = (3*G*M*Mdot)/(8*pi*r^3) * (1 - sqrt(r_isco/r))
     *
     * Standard Shakura-Sunyaev flux with GR ISCO boundary condition.
     * The factor (1 - sqrt(r_isco/r)) enforces zero stress at ISCO.
     */
    if (r <= 0.0 || mass <= 0.0) return 0.0;
    double r_isco = schwarzschild_isco_radius(schwarzschild_radius(mass));
    if (r < r_isco) return 0.0;
    double factor = 1.0 - sqrt(r_isco / r);
    if (factor <= 0.0) return 0.0;
    return (3.0 * SCHW_G_N * mass * mdot) / (8.0 * M_PI * r * r * r) * factor;
}

double schwarzschild_disk_temperature(double r, double mass, double mdot) {
    /*
     * T(r) = (F(r) / sigma)^{1/4}
     * sigma = Stefan-Boltzmann constant = 5.670367e-8 W/(m^2*K^4)
     */
    double sigma = 5.670367e-8;
    double F = schwarzschild_disk_flux(r, mass, mdot);
    if (F <= 0.0) return 0.0;
    return pow(F / sigma, 0.25);
}