/**
 * schwarzschild_extensions.h — Extensions of Schwarzschild: Kerr, RN, and beyond
 *
 * Reference: Wald (1984) Ch.6-7, Ch.12; Carroll (2004) Ch.5-6
 * MIT 8.962 — General Relativity
 *
 * Covers:
 *   L8 — Advanced Topics: Kerr metric, Reissner-Nordström, black hole thermodynamics
 */

#ifndef SCHWARZSCHILD_EXTENSIONS_H
#define SCHWARZSCHILD_EXTENSIONS_H

#include "schwarzschild_defs.h"

/* ==========================================================================
 * Reissner-Nordström (charged black hole)
 * ========================================================================== */

/**
 * Reissner-Nordström metric parameters.
 *
 * ds^2 = -Delta/r^2 * c^2*dt^2 + r^2/Delta * dr^2 + r^2*dOmega^2
 * where Delta = r^2 - 2*m*r + Q^2
 * m = GM/c^2, Q^2 = G*q^2/(4*pi*epsilon_0*c^4)
 */
typedef struct {
    double mass;
    double charge;
    double geometric_mass;
    double geometric_charge;
    double r_plus;
    double r_minus;
    int is_extremal;
    int is_naked_singularity;
} ReissnerNordstromBH;

void rn_black_hole_init(double mass, double charge, ReissnerNordstromBH *rn);
double rn_delta(double r, const ReissnerNordstromBH *rn);
double rn_g_tt(double r, const ReissnerNordstromBH *rn);
double rn_g_rr(double r, const ReissnerNordstromBH *rn);
double rn_horizon_outer(const ReissnerNordstromBH *rn);
double rn_horizon_inner(const ReissnerNordstromBH *rn);
double rn_extremal_charge(double mass);
double rn_hawking_temperature(const ReissnerNordstromBH *rn);
double rn_electric_field(double r, const ReissnerNordstromBH *rn);
void rn_metric_covariant(const SchwarzschildPoint *p, const ReissnerNordstromBH *rn,
                           SchwarzschildSymmetric4x4 *metric);
double rn_kretschmann(double r, const ReissnerNordstromBH *rn);

/* ==========================================================================
 * Kerr Metric (rotating black hole)
 * ========================================================================== */

/**
 * Kerr black hole parameters.
 *
 * ds^2 = -(1 - 2*m*r/Sigma)*c^2*dt^2 - (4*m*a*r*sin^2(theta)/Sigma)*c*dt*dphi
 *        + (Sigma/Delta)*dr^2 + Sigma*dtheta^2
 *        + (r^2 + a^2 + 2*m*a^2*r*sin^2(theta)/Sigma)*sin^2(theta)*dphi^2
 *
 * where Sigma = r^2 + a^2*cos^2(theta), Delta = r^2 - 2*m*r + a^2
 * a = J/(M*c) is the spin parameter [m]
 */
typedef struct {
    double mass;
    double spin;
    double geometric_mass;
    double geometric_spin;
    double r_plus;
    double r_minus;
    double r_ergo_equator;
    double isco_prograde;
    double isco_retrograde;
} KerrBH;

void kerr_black_hole_init(double mass, double spin, KerrBH *kerr);
double kerr_delta(double r, const KerrBH *kerr);
double kerr_sigma(double r, double theta, const KerrBH *kerr);
double kerr_g_tt(double r, double theta, const KerrBH *kerr);
double kerr_g_tphi(double r, double theta, const KerrBH *kerr);
double kerr_g_rr(double r, double theta, const KerrBH *kerr);
double kerr_g_theta_theta(double r, double theta, const KerrBH *kerr);
double kerr_g_phi_phi(double r, double theta, const KerrBH *kerr);
double kerr_horizon_outer(const KerrBH *kerr);
double kerr_horizon_inner(const KerrBH *kerr);
double kerr_ergosphere_radius(double theta, const KerrBH *kerr);
double kerr_isco_prograde(const KerrBH *kerr);
double kerr_isco_retrograde(const KerrBH *kerr);
double kerr_frame_dragging_omega(double r, double theta, const KerrBH *kerr);
double kerr_photon_sphere_prograde(const KerrBH *kerr);
double kerr_photon_sphere_retrograde(const KerrBH *kerr);
void kerr_metric_covariant(const SchwarzschildPoint *p, const KerrBH *kerr,
                             SchwarzschildSymmetric4x4 *metric);
void kerr_metric_contravariant(const SchwarzschildPoint *p, const KerrBH *kerr,
                                 SchwarzschildSymmetric4x4 *inv_metric);

/* ==========================================================================
 * Black Hole Thermodynamics
 * ========================================================================== */

/**
 * Bekenstein-Hawking entropy: S = k_B * c^3 * A / (4 * G * hbar)
 *
 * where A = event horizon area = 16*pi*m^2 for Schwarzschild.
 * This identifies the horizon area / 4 (in Planck units) as entropy,
 * linking GR (geometry) with thermodynamics (entropy).
 *
 * Reference: Bekenstein (1973), Hawking (1975)
 */
double schwarzschild_bekenstein_hawking_entropy(double mass);

/**
 * Hawking temperature: T_H = hbar * c^3 / (8 * pi * G * M * k_B)
 *
 * For a solar mass black hole: T_H ~ 6.17e-8 K (extremely cold)
 * For a primordial BH of mass 1e12 kg: T_H ~ 1.23e11 K (extremely hot,
 * would have evaporated by now)
 *
 * Reference: Hawking (1975), "Particle Creation by Black Holes"
 */
double schwarzschild_hawking_temperature(double mass);

/**
 * Black hole luminosity from Hawking radiation (Stefan-Boltzmann):
 *
 * L = sigma * A * T_H^4 = (hbar*c^6) / (15360*pi*G^2*M^2)
 *
 * Reference: Page (1976), "Particle emission rates from a black hole"
 */
double schwarzschild_hawking_luminosity(double mass);

/**
 * Black hole evaporation time (lifetime) for a Schwarzschild BH:
 *
 * t_evap ~ (5120 * pi * G^2 * M^3) / (hbar * c^4)
 *
 * For M = 1 M_sun: t_evap ~ 2.1e67 years (>> age of universe)
 * For M = 1e12 kg (primordial): t_evap ~ 2.7e10 years (~ age of universe)
 *
 * Reference: Page (1976)
 */
double schwarzschild_evaporation_time(double mass);

/**
 * First law of black hole thermodynamics:
 * dM = (kappa/(8*pi))*dA + Omega*dJ + Phi*dQ
 *
 * For Schwarzschild (dJ=0, dQ=0):
 * dM = T_H * dS = (c^2/(8*pi*G))*kappa*dA
 *
 * This function verifies the relation numerically.
 */
double schwarzschild_first_law_dM_dA(double mass);

/**
 * Second law (area theorem): the total horizon area A never decreases
 * in classical GR.
 *
 * Hawking (1971): "The event horizon area of a black hole never decreases."
 * Analogous to the second law of thermodynamics.
 *
 * This function returns the area increase when mass Delta_M is added
 * to a Schwarzschild black hole.
 */
double schwarzschild_area_increase(double mass_initial, double mass_added);

/* ==========================================================================
 * Penrose Process and Energy Extraction
 * ========================================================================== */

/**
 * Maximum energy extractable via the Penrose process from a Kerr black hole.
 *
 * A particle in the ergosphere can have negative energy (as measured at
 * infinity). If it splits into two particles, one falling in with E < 0,
 * the other escaping with E > E_initial, energy is extracted from the BH.
 *
 * Maximum efficiency for extremal Kerr (a = m): eta = 1 - 1/sqrt(2) ~ 29.3%
 *
 * Reference: Penrose (1969), "Gravitational Collapse..."
 */
double kerr_penrose_max_efficiency(double spin_param);

/**
 * Blandford-Znajek mechanism power estimate:
 *
 * P_BZ ~ (1/8) * B^2 * r_h^2 * (a/m)^2 * c
 *
 * Electromagnetic extraction of rotational energy from a Kerr black hole
 * via magnetic fields threading the ergosphere. Believed to power AGN jets.
 *
 * Reference: Blandford & Znajek (1977), MNRAS 179, 433
 */
double kerr_blandford_znajek_power(double mass, double spin,
                                     double magnetic_field);

/* ==========================================================================
 * Gravitational Wave Ringdown (Quasinormal Modes)
 * ========================================================================== */

/**
 * Fundamental quasinormal mode frequency for Schwarzschild (l=2, n=0):
 *
 * omega = (1/M) * (0.37367 - i*0.08896) in geometric units
 *
 * The real part gives the ringdown frequency, the imaginary part
 * gives the damping time.
 *
 * Reference: Chandrasekhar & Detweiler (1975), Leaver (1985)
 */
void schwarzschild_qnm_fundamental_l2_n0(double mass,
                                            double *freq_Hz,
                                            double *damping_time_s);

/**
 * Ringdown frequency for the dominant l=2, m=2 mode.
 *
 * f_lmn ~ (1/(2*pi*M)) * [1 - 0.63*(1-a)^{3/10}]
 *
 * Fitting formula from Berti, Cardoso & Will (2006), PRD 73, 064030.
 * Valid for all spin values a/m in [0, 1).
 *
 * Returns frequency in Hz.
 */
double kerr_qnm_ringdown_frequency(double mass, double spin);

/**
 * Quality factor Q = pi * f * tau for the ringdown.
 *
 * Measures how many oscillations occur during the damping time.
 * Higher Q means the BH "rings" for more cycles.
 */
double schwarzschild_qnm_quality_factor(double mass);

/**
 * Gravitational wave strain amplitude from ringdown:
 *
 * h ~ sqrt(epsilon * M / D) * exp(-t/tau) * cos(2*pi*f*t)
 *
 * where epsilon is the fraction of mass radiated (~1-3% for BBH merger).
 *
 * Returns peak strain for given parameters.
 */
double schwarzschild_ringdown_strain(double mass, double distance,
                                       double radiated_fraction);

#endif /* SCHWARZSCHILD_EXTENSIONS_H */