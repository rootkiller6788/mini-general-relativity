/**
 * @file kerr_geodesics.h
 * @brief Geodesic equations, constants of motion, and orbits in Kerr geometry
 *
 * The Kerr metric admits four constants of motion:
 * 1. mu = -1, 0, +1 for timelike, null, spacelike geodesics
 * 2. E = -p_t (energy at infinity, from stationarity)
 * 3. L_z = p_phi (azimuthal angular momentum, from axisymmetry)
 * 4. Q (Carter constant, from Killing tensor K_munu)
 *
 * The Carter constant Q makes geodesic motion in Kerr completely integrable
 * (Liouville integrable). This is a hidden symmetry from the Killing-Yano tensor.
 *
 * Reference: Carter (1968) Phys. Rev. 174, 1559; Chandrasekhar Ch.7
 * MIT 8.962, Wald Ch.7.4, Teukolsky (2015)
 */

#ifndef KERR_GEODESICS_H
#define KERR_GEODESICS_H

#include "kerr_metric.h"

/* ==========================================================================
 * L1 -- Structures
 * ========================================================================== */

/** Constants of motion for a Kerr geodesic */
typedef struct {
    double E;          /* Energy at infinity (per unit rest mass for massive) */
    double Lz;         /* Angular momentum about the symmetry axis */
    double Q;          /* Carter constant */
    double mu;         /* -1 = timelike, 0 = null, 1 = spacelike */
    double M;          /* Mass of the black hole */
    double a;          /* Spin parameter */
} KerrGeodesicConstants;

/** State vector for a Kerr geodesic */
typedef struct {
    KerrBLPoint pos;   /* Current position (t, r, theta, phi) */
    double p_t;        /* Covariant momentum components */
    double p_r;
    double p_theta;
    double p_phi;
    double affine;     /* Affine parameter (proper time for massive) */
} KerrGeodesicState;

/** ISCO (Innermost Stable Circular Orbit) data */
typedef struct {
    double r_isco;     /* ISCO radius */
    double E_isco;     /* Specific energy at ISCO */
    double Lz_isco;    /* Specific angular momentum at ISCO */
    double Omega_isco; /* Orbital angular velocity at ISCO */
    int is_prograde;   /* 1 for prograde, 0 for retrograde */
    int exists;        /* 1 if ISCO exists (a not too large for retrograde) */
} KerrISCOData;

/** Photon orbit data */
typedef struct {
    double r_photon;   /* Photon circular orbit radius (prograde/retrograde) */
    double impact_param; /* Impact parameter b = Lz/E for photon orbit */
    double carter_q;   /* Carter constant Q/E^2 for photon orbit */
    int is_prograde;   /* 1 for prograde photon orbit, 0 for retrograde */
} KerrPhotonOrbit;

/** Radial turning points for bound orbits */
typedef struct {
    double r_periastron; /* Innermost radial turning point */
    double r_apoastron;  /* Outermost radial turning point */
    int is_bound;        /* 1 if both turning points exist and are finite */
    int num_turns;       /* Number of real roots of radial potential */
} KerrRadialTurningPoints;

/** Frame dragging velocity at a point */
typedef struct {
    double omega;        /* Angular velocity of ZAMO frame */
    double v_drag;       /* Magnitude of frame-dragging effect */
    double domega_dr;    /* Radial gradient of omega */
    double domega_dtheta;/* Angular gradient of omega */
} KerrFrameDragging;

/* ==========================================================================
 * L2 -- Constants of Motion
 * ========================================================================== */

/**
 * Compute the Carter constant Q from initial conditions.
 * Q = p_theta^2 + cos^2(theta) [a^2 (mu^2 - E^2) + Lz^2 / sin^2(theta)]
 *
 * For equatorial orbits (theta = pi/2): Q = 0.
 */
double kerr_carter_constant(double E, double Lz, double p_theta,
                            double theta, double a, double mu);

/**
 * Radial potential R(r) for Kerr geodesics.
 *
 * R(r) = [E(r^2 + a^2) - a Lz]^2 - Delta [mu^2 r^2 + (Lz - aE)^2 + Q]
 *
 * The roots of R(r) = 0 determine radial turning points.
 * R(r) >= 0 is required for physical motion.
 */
double kerr_radial_potential(double r, const KerrGeodesicConstants *gc);

/**
 * Theta potential Theta(theta) for Kerr geodesics.
 *
 * Theta(theta) = Q - cos^2(theta) [a^2(mu^2 - E^2) + Lz^2/sin^2(theta)]
 *
 * Theta(theta) >= 0 is required for physical motion in theta.
 */
double kerr_theta_potential(double theta, const KerrGeodesicConstants *gc);

/**
 * Compute all 4 constants of motion from initial position and 4-velocity.
 * Uses the Killing vectors and Killing tensor.
 */
int kerr_compute_constants(const KerrBlackHole *bh, const KerrBLPoint *pos,
                           const double u[4], double mu,
                           KerrGeodesicConstants *gc);

/**
 * Compute the 4 covariant momentum components from constants of motion.
 * p_t = -E, p_phi = Lz
 * p_r^2 = R(r)/Delta^2, p_theta^2 = Theta(theta)
 */
int kerr_momenta_from_constants(const KerrBlackHole *bh,
                                const KerrBLPoint *pos,
                                const KerrGeodesicConstants *gc,
                                int r_sign, int theta_sign,
                                double p[4]);

/* ==========================================================================
 * L3 -- Geodesic Integration
 * ========================================================================== */

/**
 * Single RK4 step for Kerr geodesic integration using Carter equations.
 *
 * The first-order form separates the r and theta equations:
 *   Sigma dr/dlambda = +/- sqrt(R(r))
 *   Sigma dtheta/dlambda = +/- sqrt(Theta(theta))
 *   Sigma dphi/dlambda = -(aE - Lz/sin^2theta) + (a/Delta) P(r)
 *   Sigma dt/dlambda = -a(aE sin^2theta - Lz) + (r^2+a^2)/Delta P(r)
 * where P(r) = E(r^2 + a^2) - a Lz
 */
int kerr_geodesic_rk4_step(const KerrBlackHole *bh,
                           KerrGeodesicState *state,
                           double dlambda,
                           const KerrGeodesicConstants *gc);

/**
 * Integrate a geodesic from start position over affine parameter range.
 * Uses adaptive RK4(5) Cash-Karp method.
 *
 * @param bh      Black hole parameters
 * @param state   Initial state (updated in place)
 * @param n_steps Number of integration steps
 * @param dlambda Step size in affine parameter
 * @param gc      Constants of motion
 * @return 0 on success, 1 if geodesic enters singularity
 */
int kerr_geodesic_integrate(const KerrBlackHole *bh,
                            KerrGeodesicState *state,
                            int n_steps, double dlambda,
                            const KerrGeodesicConstants *gc);

/**
 * Check whether a geodesic will cross the event horizon.
 * Determined by examining the radial potential at the horizon.
 */
int kerr_geodesic_crosses_horizon(const KerrGeodesicConstants *gc,
                                  const KerrBlackHole *bh);

/* ==========================================================================
 * L4 -- ISCO and Circular Orbits
 * ========================================================================== */

/**
 * Compute the ISCO radius for Kerr.
 *
 * For Schwarzschild limit (a=0): r_isco = 6M
 * For extremal Kerr prograde (a=M): r_isco = M (degenerate with horizon)
 * For extremal Kerr retrograde (a=-M): r_isco = 9M
 *
 * The ISCO satisfies: R(r)=0, R'(r)=0, R''(r)=0 simultaneously.
 *
 * Reference: Bardeen, Press & Teukolsky (1972), ApJ 178, 347
 */
int kerr_isco(const KerrBlackHole *bh, int prograde, KerrISCOData *isco);

/**
 * Compute ISCO properties analytically for the Schwarzschild limit a=0.
 */
KerrISCOData kerr_isco_schwarzschild(double M);

/**
 * ISCO radius as a function of spin, fitting formula (Bardeen+1972).
 * r_isco/M = 3 + Z2 - sgn(a) * sqrt((3-Z1)(3+Z1+2Z2))
 * where Z1 = 1 + (1-a^2/M^2)^(1/3)[(1+a/M)^(1/3) + (1-a/M)^(1/3)]
 *       Z2 = sqrt(3 a^2/M^2 + Z1^2)
 */
double kerr_isco_radius_bardeen(double M, double a, int prograde);

/**
 * ISCO orbital frequency Omega_isco for a Kerr black hole.
 * Omega = 1 / (a + r^(3/2)/sqrt(M)) in geometric units.
 */
double kerr_isco_frequency(const KerrBlackHole *bh, int prograde);

/* ==========================================================================
 * L5 -- Photon Orbits
 * ========================================================================== */

/**
 * Compute photon circular orbit radii in Kerr.
 *
 * The photon orbit satisfies:
 * r_ph = 2M [1 + cos(2/3 arccos(-a/M))]  for prograde
 * r_ph = 2M [1 + cos(2/3 arccos(+a/M))]  for retrograde
 *
 * Schwarzschild (a=0): r_ph = 3M (photon sphere)
 * Extremal prograde (a=M): r_ph = M
 * Extremal retrograde: r_ph = 4M
 */
int kerr_photon_orbit(const KerrBlackHole *bh, int prograde,
                      KerrPhotonOrbit *orbit);

/**
 * Photon orbit impact parameter b = Lz/E.
 * For Schwarzschild: b = 3*sqrt(3) M.
 * The shadow of a black hole is determined by photon orbits.
 */
double kerr_photon_impact_parameter(const KerrBlackHole *bh, int prograde);

/**
 * Compute the black hole shadow shape as seen by a distant observer.
 * The shadow boundary is determined by photon spherical orbits.
 *
 * Returns the cartesian shadow coordinates (alpha, beta) on the observer''s
 * image plane for viewing inclination iota.
 */
int kerr_shadow_boundary(const KerrBlackHole *bh, double iota,
                         int n_angles, double *alpha, double *beta);

/* ==========================================================================
 * L6 -- Frame Dragging
 * ========================================================================== */

/**
 * Frame-dragging angular velocity: omega = -g_tphi / g_phiphi.
 *
 * For ZAMOs (Zero Angular Momentum Observers), this is the angular
 * velocity they acquire due to frame dragging.
 */
double kerr_frame_dragging_omega(const KerrBlackHole *bh,
                                 const KerrBLPoint *pt);

/**
 * Lense-Thirring precession frequency at large r (weak-field limit):
 * Omega_LT = 2 J / r^3 = 2 M a / r^3.
 *
 * This is the rate at which gyroscopes precess in the Kerr field.
 * Key prediction of GR verified by Gravity Probe B.
 */
double kerr_lense_thirring_frequency(const KerrBlackHole *bh, double r);

/**
 * Compute full frame-dragging data at a point.
 */
int kerr_frame_dragging_data(const KerrBlackHole *bh,
                             const KerrBLPoint *pt,
                             KerrFrameDragging *fd);

/**
 * Compute the angular velocity of the "locally non-rotating frame"
 * (LNRO/ZAMO frame). Same as frame_dragging_omega but returns
 * the full 4-velocity components.
 */
int kerr_zamo_four_velocity(const KerrBlackHole *bh, const KerrBLPoint *pt,
                            double u[4]);

/* ==========================================================================
 * Radial and Polar Motion Classification
 * ========================================================================== */

/**
 * Classify the type of radial motion:
 * Returns: 0=bounded, 1=plunging, 2=escape, -1=circular (unstable),
 *          -2=circular (stable), -3=no valid motion
 */
int kerr_classify_radial_motion(const KerrGeodesicConstants *gc);

/**
 * Determine radial turning points for given constants of motion.
 * Computes all real positive roots of the radial potential R(r).
 */
int kerr_radial_turning_points(const KerrGeodesicConstants *gc,
                               KerrRadialTurningPoints *turns);

/**
 * Compute the radial epicyclic frequency kappa_r for nearly circular orbits.
 * For Kerr: kappa_r^2 = Omega^2 (1 - 6M/r + 8a M^(1/2)/r^(3/2) - 3a^2/r^2)
 *
 * Important for modeling QPOs (Quasi-Periodic Oscillations) in BH systems.
 */
double kerr_epicyclic_frequency_r(const KerrBlackHole *bh, double r,
                                  int prograde);

/**
 * Compute the vertical epicyclic frequency Omega_theta for nearly
 * circular orbits.
 *
 * Omega_theta^2 = Omega^2 (1 - 4a M^(1/2)/r^(3/2) + 3a^2/r^2)
 */
double kerr_epicyclic_frequency_theta(const KerrBlackHole *bh, double r,
                                       int prograde);

/**
 * Keplerian orbital frequency for circular equatorial orbits.
 * Omega_K = 1 / (a + r^(3/2) / sqrt(M))
 */
double kerr_keplerian_frequency(const KerrBlackHole *bh, double r,
                                int prograde);

/**
 * Periastron precession frequency: Omega_peri = Omega - kappa_r
 */
double kerr_periastron_precession(const KerrBlackHole *bh, double r,
                                  int prograde);

/**
 * Nodal precession frequency (Lense-Thirring): Omega_nodal = Omega - Omega_theta
 */
double kerr_nodal_precession(const KerrBlackHole *bh, double r,
                             int prograde);

#endif /* KERR_GEODESICS_H */