/**
 * schwarzschild_geodesics.h — Geodesic equations in Schwarzschild spacetime
 *
 * Reference: Wald (1984) Ch.6; Carroll (2004) Ch.5; MTW (1973) Ch.25
 * MIT 8.962 — General Relativity
 *
 * Covers:
 *   L2 — Core Concepts: geodesic motion, constants of motion
 *   L4 — Fundamental Laws: geodesic equation from metric
 *   L6 — Canonical Systems: Keplerian orbits, ISCO, photon sphere
 */

#ifndef SCHWARZSCHILD_GEODESICS_H
#define SCHWARZSCHILD_GEODESICS_H

#include "schwarzschild_defs.h"

/* Orbit types */
typedef enum {
    ORBIT_BOUND,
    ORBIT_UNBOUND,
    ORBIT_MARGINAL,
    ORBIT_PLUNGE,
    ORBIT_CIRCULAR,
    ORBIT_UNSTABLE_CIRCULAR
} SchwarzschildOrbitType;

/* Constants of motion */
typedef struct {
    double E;
    double L;
    double Q;
} SchwarzschildConstants;

/* Initial conditions for a geodesic */
typedef struct {
    SchwarzschildPoint initial_point;
    double dt_dlambda;
    double dr_dlambda;
    double dtheta_dlambda;
    double dphi_dlambda;
} SchwarzschildGeodesicInit;

/* Effective potential parameters */
typedef struct {
    double r;
    double V_eff;
    double dV_dr;
    double d2V_dr2;
} SchwarzschildEffectivePot;

/* ===== Effective Potential ===== */

/* V_eff(r) = (1 - rs/r)(1 + L^2/r^2) for massive particles (per unit mass) */
double schwarzschild_V_eff_massive(double r, double rs, double L_sq);

/* V_eff(r) = (1 - rs/r) * L^2/r^2 for massless particles */
double schwarzschild_V_eff_massless(double r, double rs, double L_sq);

/* dV_eff/dr for massive particles */
double schwarzschild_dV_eff_dr_massive(double r, double rs, double L_sq);

/* dV_eff/dr for massless particles */
double schwarzschild_dV_eff_dr_massless(double r, double rs, double L_sq);

/* d^2V_eff/dr^2 for massive particles */
double schwarzschild_d2V_eff_dr2_massive(double r, double rs, double L_sq);

/* Compute full effective potential info at radius r */
void schwarzschild_effective_potential_massive(double r, double rs, double L_sq,
                                                SchwarzschildEffectivePot *pot);

void schwarzschild_effective_potential_massless(double r, double rs, double L_sq,
                                                 SchwarzschildEffectivePot *pot);

/* ===== Circular Orbits ===== */

/* Circular orbit radius from angular momentum L:
   r_circ solves dV_eff/dr = 0 => r_circ = (L^2/(2*m))*(1 + sqrt(1 - 12*m^2/L^2))
   For given L, returns the stable circular orbit radius. */
double schwarzschild_circular_orbit_r_from_L(double L, double rs);

/* Angular momentum for circular orbit at radius r:
   L_circ^2 = m * r^2 / (r - 3*m) */
double schwarzschild_L_for_circular_orbit(double r, double rs);

/* Energy for circular orbit at radius r:
   E_circ^2 = (r - 2*m)^2 / (r*(r - 3*m)) */
double schwarzschild_E_for_circular_orbit(double r, double rs);

/* Check if circular orbit at r is stable: d^2V_eff/dr^2 > 0 */
int schwarzschild_is_stable_circular_orbit(double r, double rs);

/* ===== Special Radii ===== */

/* Photon sphere radius: r_ph = 3*m = 1.5*rs */
double schwarzschild_photon_sphere_radius(double rs);

/* ISCO radius: r_isco = 6*m = 3*rs */
double schwarzschild_isco_radius(double rs);

/* Marginally bound orbit radius: r_mb = 4*m = 2*rs */
double schwarzschild_marginally_bound_radius(double rs);

/* ISCO angular momentum per unit rest mass: L_isco = sqrt(12)*m */
double schwarzschild_L_isco(double rs);

/* ISCO energy per unit rest mass: E_isco = sqrt(8/9) */
double schwarzschild_E_isco(void);

/* ISCO angular velocity (Keplerian): Omega_isco = c^3/(6^(3/2)*G*M) */
double schwarzschild_Omega_isco(double mass);

/* ===== Orbital period ===== */

/* Keplerian orbital period for circular orbit: T = 2*pi*sqrt(r^3/(G*M)) */
double schwarzschild_orbital_period_keplerian(double r, double mass);

/* Relativistic orbital period (including GR corrections):
   T_GR = T_Kepler * sqrt(1 - 3*m/r) / sqrt(1 - 2*m/r) */
double schwarzschild_orbital_period_relativistic(double r, double mass);

/* ===== Radial geodesic equation (d^2r/dtau^2) ===== */

/* Radial acceleration for massive particle:
   d^2r/dtau^2 = -m/r^2 + L^2/r^3 - 3*m*L^2/r^4 */
double schwarzschild_radial_accel_massive(double r, double rs, double L_sq);

/* Radial acceleration for massless particle */
double schwarzschild_radial_accel_massless(double r, double rs, double L_sq);

/* ===== Radial motion classification ===== */

/* Classify the type of orbit based on E, L and position */
SchwarzschildOrbitType schwarzschild_classify_orbit(double r, double E_sq,
                                                      double L_sq, double rs);

/* Radial turning points: solve E^2 = V_eff(r) for r.
   Returns number of positive real roots (up to 4). */
int schwarzschild_radial_turning_points(double E_sq, double L_sq, double rs,
                                         double roots[4]);

/* ===== Angular motion (theta equation) ===== */

/* Azimuthal equation of motion: dphi/dtau = L / r^2 */
double schwarzschild_dphi_dtau_massive(double r, double L);

/* Theta equation for equatorial plane (theta = pi/2): dtheta/dtau = 0 */
/* Non-equatorial: dtheta/dtau = sqrt(Q - L^2/tan^2(theta)) / r^2 */
double schwarzschild_dtheta_dtau(double r, double theta, double Q, double L_sq);

/* ===== Proper time vs coordinate time ===== */

/* dt/dtau for massive particle: dt/dtau = E / (1 - rs/r) */
double schwarzschild_dt_dtau_massive(double E, double r, double rs);

/* Total proper time for radial infall from r0 to r1:
   tau = integral_{r0}^{r1} dr / sqrt(E^2 - (1-rs/r)(1+L^2/r^2)) */
double schwarzschild_proper_time_radial_infall(double r_start, double r_end,
                                                double rs);

/* Coordinate time for radial infall (diverges at horizon) */
double schwarzschild_coordinate_time_radial_infall(double r_start, double r_end,
                                                    double rs);

/* ===== Complete first-order geodesic system ===== */

/* Get the 8-component state vector derivative dy/dlambda for a geodesic.
   State: (t, r, theta, phi, p_t, p_r, p_theta, p_phi)
   where p_mu is the covariant 4-momentum (lower index).
   This is the Hamiltonian form which avoids the affine parameter issues
   for null geodesics.
   
   Reference: Levin & Perez-Giz (2008), "A Periodic Table for Black Hole Orbits" */
void schwarzschild_geodesic_derivs(double lambda, const double y[8],
                                    double dydlambda[8], double rs);

/* ===== Affine parameter for null geodesics ===== */

/* Normalize affine parameter so that E = 1 (useful for ray tracing) */
void schwarzschild_normalize_null_geodesic(double *E, double *L, double impact_param);

#endif /* SCHWARZSCHILD_GEODESICS_H */