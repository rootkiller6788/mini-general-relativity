/**
 * schwarzschild_observables.h — Observable effects of Schwarzschild spacetime
 *
 * Reference: Wald (1984) Ch.6; Carroll (2004) Ch.5; Will (2018) "Gravity"
 * MIT 8.962 — General Relativity
 *
 * Covers:
 *   L2 — Core Concepts: time dilation, gravitational redshift
 *   L6 — Canonical Systems: light deflection, perihelion precession, Shapiro delay
 *   L7 — Applications: GPS corrections, astrophysical observations
 */

#ifndef SCHWARZSCHILD_OBSERVABLES_H
#define SCHWARZSCHILD_OBSERVABLES_H

#include "schwarzschild_defs.h"

/* ==========================================================================
 * Gravitational Time Dilation
 * ========================================================================== */

/**
 * Gravitational time dilation factor between a clock at radius r and infinity.
 *
 * dtau = dt * sqrt(1 - rs/r)
 *
 * Clocks deeper in the gravitational well run slower.
 * For r >> rs, this reduces to: dtau/dt ~ 1 - rs/(2r) = 1 - GM/(c^2 r)
 *
 * Reference: Einstein (1907), "On the Relativity Principle..."
 */
double schwarzschild_time_dilation_factor(double r, double rs);

/**
 * Time dilation between two clocks at radii r1 and r2:
 * tau1 / tau2 = sqrt((1 - rs/r1) / (1 - rs/r2))
 */
double schwarzschild_time_dilation_ratio(double r1, double r2, double rs);

/**
 * Time dilation experienced by a GPS satellite at Earth altitude.
 *
 * Combines special relativistic (velocity) and general relativistic
 * (gravitational) effects:
 *   Delta_t/t = -GM/(c^2 * R_E) * (1 - R_E/r) - (3/2) * GM/(c^2 * R_E) * (1 - R_E/r) + ...
 *
 * For GPS (r ~ 26600 km from Earth center, R_E ~ 6371 km):
 *   GR: +45.7 us/day (satellite clocks run faster)
 *   SR: -7.2 us/day (satellite clocks run slower due to velocity)
 *   Net: +38.5 us/day
 *
 * Without correction, GPS error accumulates ~11 km/day.
 *
 * Reference: Ashby (2003), "Relativity in the Global Positioning System"
 */
double schwarzschild_gps_time_dilation(double r_satellite, double r_earth,
                                        double mass_earth, double v_satellite);

/* ==========================================================================
 * Gravitational Redshift
 * ========================================================================== */

/**
 * Gravitational redshift: light emitted at r_e observed at infinity.
 *
 * z = lambda_obs / lambda_em - 1 = 1/sqrt(1 - rs/r_e) - 1
 *
 * For weak field (r_e >> rs): z ~ GM/(c^2 r_e)
 */
double schwarzschild_redshift_at_infinity(double r_emitter, double rs);

/**
 * Redshift between emitter at r_e and observer at r_o:
 * 1 + z = sqrt((1 - rs/r_o) / (1 - rs/r_e))
 *
 * If r_o > r_e (observer farther from source), z > 0 (redshift).
 * If r_o < r_e (observer closer), z < 0 (blueshift).
 */
double schwarzschild_redshift_finite(double r_emitter, double r_observer, double rs);

/**
 * Redshift of light emitted from ISCO and observed at infinity.
 * This is a key observational signature of accretion disks around
 * black holes — the "red wing" of iron K-alpha lines.
 *
 * Reference: Fabian et al. (1989), "X-ray fluorescence from the inner disc..."
 */
double schwarzschild_redshift_isco_at_infinity(double rs);

/**
 * Doppler shift for a source in circular orbit at radius r,
 * observed at infinity along the line of sight.
 *
 * g = E_obs / E_em = 1 / (u^t * (1 - Omega * b))
 *
 * where u^t = 1/sqrt(1 - 3*m/r), Omega = sqrt(m/r^3), b is impact parameter.
 */
double schwarzschild_doppler_circular_orbit(double r, double rs,
                                              double impact_param, double incl_angle);

/* ==========================================================================
 * Light Deflection (Gravitational Lensing)
 * ========================================================================== */

/**
 * Einstein deflection angle for light passing near a mass M with
 * impact parameter b:
 *
 * alpha = 4*G*M / (b * c^2) = 2*rs / b   [radians]
 *
 * This is the first-order (weak field) approximation.
 * For the Sun: alpha_sun = 1.75 arcseconds for a ray grazing the solar limb.
 *
 * The factor 2 difference from Newtonian expectation (2*GM/(b*c^2))
 * was the first experimental test of GR (Eddington 1919).
 *
 * Reference: Einstein (1915), Dyson, Eddington & Davidson (1920)
 */
double schwarzschild_light_deflection_angle(double impact_param, double mass);

/**
 * Exact light deflection angle (strong field, elliptic integral):
 *
 * alpha = 2 * integral_{r0}^{infinity} dr / (r * sqrt(r^2/b^2 - (1 - rs/r)))
 *          - pi
 *
 * where r0 is the turning point: r0^3 - b^2*(r0 - rs) = 0
 */
double schwarzschild_light_deflection_exact(double impact_param, double rs,
                                             double tolerance);

/**
 * Einstein ring angular radius for a point mass lens:
 *
 * theta_E = sqrt(4*G*M/c^2 * D_ls/(D_l * D_s))
 *
 * where D_l, D_s, D_ls are distances to lens, source, and lens-source.
 *
 * Reference: Schneider, Ehlers & Falco (1992), "Gravitational Lenses"
 */
double schwarzschild_einstein_ring_radius(double mass, double D_l,
                                            double D_s, double D_ls);

/**
 * Magnification factor for a point mass gravitational lens:
 *
 * mu = (u^2 + 2) / (u * sqrt(u^2 + 4))
 *
 * where u = beta/theta_E is the angular separation in units of Einstein radius.
 * For perfect alignment (u -> 0): mu -> infinity (Einstein ring).
 */
double schwarzschild_lens_magnification(double angular_separation_thetaE_units);

/**
 * Image positions for a point mass lens (lens equation solutions):
 *
 * theta_{1,2} = (beta +/- sqrt(beta^2 + 4*theta_E^2)) / 2
 *
 * Two images: one outside, one inside the Einstein ring.
 * Returns the number of images (always 2 for point mass, except beta=0 yields ring).
 */
int schwarzschild_lens_image_positions(double beta, double theta_E,
                                         double *theta_plus, double *theta_minus);

/* ==========================================================================
 * Shapiro Time Delay
 * ========================================================================== */

/**
 * Shapiro time delay for a radar signal passing near mass M.
 *
 * The coordinate time delay (relative to flat spacetime) is:
 * Delta_t = 2*GM/c^3 * ln(4 * r_earth * r_target / b^2)
 *
 * For a round trip Earth -> target -> Earth, with closest approach b.
 *
 * Verified by Shapiro et al. (1971) using radar ranging to Venus
 * and Mercury, confirming GR at the 5% level.
 *
 * Reference: Shapiro (1964), "Fourth Test of General Relativity"
 */
double schwarzschild_shapiro_time_delay(double r_emitter, double r_reflector,
                                          double impact_param, double mass);

/**
 * Shapiro delay for a planetary radar experiment (Earth-Venus).
 * Returns the excess round-trip time in seconds.
 */
double schwarzschild_shapiro_earth_venus(double earth_orbit_radius,
                                           double venus_orbit_radius,
                                           double solar_mass);

/**
 * Shapiro delay logarithmic term (the part beyond the geometric path):
 * Delta_t = 2*m/c * ln(4*r_e*r_r / b^2)
 * where m = GM/c^2 (geometric mass).
 */
double schwarzschild_shapiro_log_term(double r_emitter, double r_reflector,
                                        double impact_param, double geometric_mass);

/* ==========================================================================
 * Perihelion/Periastron Precession
 * ========================================================================== */

/**
 * Perihelion precession per orbit for a test particle on a bound
 * geodesic (eccentric orbit):
 *
 * Delta_phi = 6*pi*G*M / (a * c^2 * (1 - e^2))   [radians per orbit]
 *
 * where a = semi-major axis, e = eccentricity.
 *
 * Mercury: Delta_phi ~ 43 arcseconds per century (GR contribution,
 * after accounting for Newtonian perturbations from other planets).
 *
 * Reference: Einstein (1915), "Erklarung der Perihelbewegung des Merkur..."
 */
double schwarzschild_perihelion_precession_per_orbit(double semi_major_axis,
                                                       double eccentricity,
                                                       double mass);

/**
 * Perihelion precession rate: dphi/dt = Omega * Delta_phi / (2*pi)
 * where Omega = sqrt(G*M / a^3) is the mean motion.
 *
 * Returns precession in radians per second.
 */
double schwarzschild_perihelion_precession_rate(double semi_major_axis,
                                                  double eccentricity,
                                                  double mass);

/**
 * Precession of the longitude of periastron for the Hulse-Taylor
 * binary pulsar PSR B1913+16.
 *
 * Observed: ~4.2 degrees per year
 * GR prediction: ~4.226 degrees per year
 *
 * This was one of the earliest confirmations of GR in strong-field
 * binary systems, and the orbital decay due to gravitational radiation
 * provided the first indirect detection of gravitational waves
 * (Nobel Prize 1993: Hulse & Taylor).
 *
 * Reference: Taylor & Weisberg (1982)
 */
double schwarzschild_hulse_taylor_precession(double m1, double m2,
                                               double semi_major_axis,
                                               double eccentricity);

/* ==========================================================================
 * Black Hole Shadow
 * ========================================================================== */

/**
 * Angular radius of the black hole shadow for a Schwarzschild BH:
 *
 * theta_shadow = sqrt(27) * m / D = 3*sqrt(3) * m / D
 *
 * where D is the distance to the black hole.
 *
 * For M87* (M ~ 6.5e9 Msun, D ~ 16.8 Mpc):
 * theta_shadow ~ 42 microarcseconds
 *
 * This was directly imaged by the Event Horizon Telescope (EHT) in 2019.
 *
 * Reference: Event Horizon Telescope Collaboration (2019), ApJL 875, L1
 */
double schwarzschild_shadow_angular_radius(double mass, double distance);

/**
 * Photon capture cross section: sigma = 27*pi*m^2
 * (photons with impact parameter b < 3*sqrt(3)*m are captured)
 *
 * Reference: Synge (1966), "The Escape of Photons from Gravitationally Intense Stars"
 */
double schwarzschild_photon_capture_cross_section(double mass);

/**
 * Critical impact parameter for photon capture:
 * b_crit = 3*sqrt(3) * m = (3*sqrt(3)/2) * rs
 *
 * Photons with b < b_crit fall into the black hole.
 * Photons with b > b_crit are deflected.
 * Exactly at b = b_crit, photons orbit on the photon sphere (unstable).
 */
double schwarzschild_critical_impact_parameter(double rs);

/* ==========================================================================
 * Accretion Disk Observables
 * ========================================================================== */

/**
 * Radiative efficiency of a thin accretion disk:
 * eta = 1 - E_isco
 *
 * For Schwarzschild: eta = 1 - sqrt(8/9) ~ 0.057 (5.7%)
 * For maximally rotating Kerr: eta ~ 0.42 (42%)
 *
 * The energy released as matter spirals from infinity to ISCO is
 * converted to radiation, heating the accretion disk.
 *
 * Reference: Novikov & Thorne (1973), "Astrophysics of Black Holes"
 */
double schwarzschild_accretion_efficiency(void);

/**
 * Temperature profile of a thin accretion disk (Shakura-Sunyaev model):
 *
 * T(r) = [3*G*M*Mdot / (8*pi*r^3*sigma) * (1 - sqrt(r_isco/r))]^{1/4}
 *
 * where Mdot = accretion rate, sigma = Stefan-Boltzmann constant.
 *
 * Reference: Shakura & Sunyaev (1973), A&A 24, 337
 */
double schwarzschild_disk_temperature(double r, double mass, double mdot);

/**
 * Flux from thin accretion disk at radius r (Newtonian approximation with
 * GR boundary condition at ISCO):
 *
 * F(r) = (3*G*M*Mdot) / (8*pi*r^3) * (1 - sqrt(r_isco/r))
 */
double schwarzschild_disk_flux(double r, double mass, double mdot);

/**
 * Inner edge of a Schwarzschild accretion disk: r_in = r_isco = 6*m.
 * Matter inside ISCO plunges rapidly into the black hole,
 * contributing negligible radiation.
 */
double schwarzschild_disk_inner_edge(double rs);

#endif /* SCHWARZSCHILD_OBSERVABLES_H */