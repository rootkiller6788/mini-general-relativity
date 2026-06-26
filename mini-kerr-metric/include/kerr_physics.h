/**
 * @file kerr_physics.h
 * @brief Physical processes in Kerr geometry: energy extraction,
 *        superradiance, Blandford-Znajek mechanism, and accretion
 *
 * The Kerr geometry enables several unique astrophysical processes:
 * - Penrose process: extract rotational energy from inside the ergosphere
 * - Superradiance: wave amplification by scattering off rotating BH
 * - Blandford-Znajek: electromagnetic energy extraction via magnetic fields
 * - Accretion disk physics: thin disk model (Novikov-Thorne)
 *
 * Reference: Penrose & Floyd (1971), Blandford & Znajek (1977)
 * MIT 8.962, Cambridge Part III Astro, Princeton PHY 535
 */

#ifndef KERR_PHYSICS_H
#define KERR_PHYSICS_H

#include "kerr_metric.h"
#include "kerr_horizons.h"
#include "kerr_geodesics.h"

/* ==========================================================================
 * L1 -- Energy Extraction Structures
 * ========================================================================== */

/** Penrose process configuration */
typedef struct {
    double r_decay;       /* Radial coordinate of particle decay */
    double theta_decay;   /* Polar angle of decay */
    double E_in;          /* Energy of incoming particle */
    double Lz_in;         /* Angular momentum of incoming particle */
    double E_out;         /* Energy of escaping particle */
    double E_negative;    /* Energy of particle falling into BH */
    double efficiency;    /* eta = (E_out - E_in) / E_in */
    int feasible;         /* 1 if energy extraction is possible */
} PenroseProcess;

/** Blandford-Znajek process configuration */
typedef struct {
    double B_field;       /* Poloidal magnetic field strength at horizon */
    double Omega_F;       /* Field line angular velocity (~ Omega_H / 2) */
    double power;         /* Extracted power P_BZ */
    double efficiency;    /* Efficiency relative to B^2 M^2 */
    double jet_power;     /* Observable jet power estimate */
} BlandfordZnajek;

/** Superradiance condition data */
typedef struct {
    double omega;         /* Wave frequency */
    double m;             /* Azimuthal mode number */
    double omega_H;       /* Horizon angular frequency */
    double threshold;     /* omega < m * omega_H */
    int is_superradiant;  /* 1 if superradiance condition holds */
    double amplification; /* Amplification factor (if computable) */
} SuperradianceCondition;

/** Thin accretion disk model (Novikov-Thorne) */
typedef struct {
    double M_dot;         /* Mass accretion rate */
    double r_in;          /* Inner disk edge (ISCO) */
    double r_out;         /* Outer disk edge */
    double L_disk;        /* Total disk luminosity */
    double T_max;         /* Maximum disk temperature */
    double eta;           /* Radiative efficiency eta = 1 - E_isco */
    double r_max_flux;    /* Radius of maximum flux */
} AccretionDisk;

/** Black hole shadow observables */
typedef struct {
    double shadow_radius;    /* Characteristic shadow size */
    double shadow_deviation; /* Deviation from Schwarzschild shadow */
    double asymmetry;        /* Shadow asymmetry parameter */
    double spin_estimate;    /* Spin estimate from shadow shape */
} KerrShadow;

/* ==========================================================================
 * L2 -- Penrose Process
 * ========================================================================== */

/**
 * Maximum efficiency of the Penrose process.
 *
 * For a maximally rotating (a=M) Kerr BH:
 * eta_max = 1 - 1/sqrt(2) ~ 0.2929 (about 29% of rest mass energy)
 *
 * The process requires a particle to decay inside the ergosphere,
 * with one fragment falling into the BH with negative energy.
 */
double kerr_penrose_max_efficiency(const KerrBlackHole *bh);

/**
 * Check if the Penrose process is feasible at a given point.
 * Requires: inside ergosphere, particle decay can produce
 * one negative-energy fragment (E < 0 relative to infinity).
 */
int kerr_penrose_feasibility(const KerrBlackHole *bh,
                             const KerrBLPoint *pt,
                             double E_particle, double Lz_particle);

/**
 * Compute Penrose process outcome given decay parameters.
 * Simulates a particle of rest mass mu decaying into two fragments
 * at the specified location.
 */
int kerr_penrose_process(const KerrBlackHole *bh,
                         const KerrBLPoint *decay_pt,
                         double E_in, double Lz_in, double mu_in,
                         double mu_frag1, double mu_frag2,
                         PenroseProcess *result);

/**
 * Collisional Penrose process: two particles collide in ergosphere.
 * Maximum CM energy can be arbitrarily large near the horizon
 * for extremal Kerr (Banados, Silk & West 2009).
 *
 * Returns: CM energy squared in units of rest mass.
 */
double kerr_collisional_cm_energy(const KerrBlackHole *bh,
                                  const KerrBLPoint *pt,
                                  const KerrGeodesicConstants *gc1,
                                  const KerrGeodesicConstants *gc2);

/* ==========================================================================
 * L3 -- Superradiance
 * ========================================================================== */

/**
 * Superradiance condition for a scalar field with frequency omega
 * and azimuthal quantum number m:
 *
 * omega < m * Omega_H  (where Omega_H = a / (r_plus^2 + a^2))
 *
 * When satisfied, the reflected wave has larger amplitude than
 * the incident wave, extracting rotational energy from the BH.
 */
int kerr_superradiance_condition(const KerrBlackHole *bh,
                                 double omega, int m,
                                 SuperradianceCondition *cond);

/**
 * Scalar wave superradiance amplification factor in the
 * low-frequency limit (omega M << 1).
 *
 * Z_abs ~ 4 (omega - m Omega_H) omega (r_plus^2 + a^2)
 *
 * Reference: Starobinsky (1973), Teukolsky & Press (1974)
 */
double kerr_superradiance_amplification_scalar(const KerrBlackHole *bh,
                                                double omega, int m,
                                                int l);

/**
 * Compute the maximum superradiant amplification for scalar waves.
 * Scans over frequency omega to find peak amplification.
 */
double kerr_superradiance_max_gain_scalar(const KerrBlackHole *bh,
                                          int l, int m);

/**
 * Superradiant instability timescale for a massive scalar field
 * forming a "gravitational atom" around the BH.
 *
 * tau_inst ~ 10^7 M exp(1.84 M mu) years for a/M ~ 0.99
 *
 * Reference: Dolan (2007), Brito, Cardoso & Pani (2015)
 */
double kerr_superradiant_instability_timescale(const KerrBlackHole *bh,
                                                double mu_scalar);

/* ==========================================================================
 * L4 -- Blandford-Znajek Mechanism
 * ========================================================================== */

/**
 * Blandford-Znajek electromagnetic power:
 *
 * P_BZ ~ (kappa / 4pi c) Phi_BH^2 Omega_H (Omega_F) (1 - k Omega_F/Omega_H)
 *
 * where Omega_F ~ Omega_H/2 (optimal), k ~ 1 for split monopole.
 *
 * Simplified: P_BZ = (1/6pi) B^2 M^2 a^2 / M^4  (geometric units)
 * In SI: P_BZ ~ 10^45 (B/10^4 G)^2 (M/10^9 Msol)^2 (a/M)^2 erg/s
 */
int kerr_blandford_znajek(const KerrBlackHole *bh, double B_poloidal,
                          BlandfordZnajek *bz);

/**
 * Maximum BZ power (for optimal field configuration, Omega_F = Omega_H/2):
 * P_max = (1/6) B^2 r_plus^2 a^2 Omega_H
 */
double kerr_bz_max_power(const KerrBlackHole *bh, double B_poloidal);

/**
 * BZ jet power in physical units given BH mass in solar masses
 * and B field in Gauss.
 *
 * Typical AGN: M=10^9 Msol, B=10^4 G => P ~ 10^45 erg/s
 */
double kerr_bz_power_physical(double M_solar, double a_spin,
                               double B_gauss);

/**
 * Estimate magnetic field at the horizon from observed jet power
 * (inverse BZ problem).
 */
double kerr_bz_estimate_B_field(double M_solar, double a_spin,
                                 double jet_power_erg_s);

/* ==========================================================================
 * L5 -- Accretion Disk Physics
 * ========================================================================== */

/**
 * Novikov-Thorne thin accretion disk model for Kerr.
 *
 * Radiative efficiency: eta = 1 - E_isco
 * For Schwarzschild: eta ~ 0.057
 * For extremal Kerr prograde: eta ~ 0.423
 */
int kerr_accretion_disk(const KerrBlackHole *bh, double M_dot,
                        AccretionDisk *disk);

/**
 * Disk flux as a function of radius (Novikov-Thorne formula):
 * F(r) = (M_dot / 4pi) * f(r, a) / M^2
 *
 * The function f(r,a) is given by the relativistic correction factor.
 */
double kerr_disk_flux_radius(const KerrBlackHole *bh, double M_dot, double r);

/**
 * Disk temperature at radius r (assuming blackbody):
 * T(r) = (F(r) / sigma_SB)^(1/4)
 */
double kerr_disk_temperature(const KerrBlackHole *bh, double M_dot, double r);

/**
 * Total disk spectrum: compute F_nu by integrating over disk radii
 * with relativistic Doppler and gravitational redshift corrections.
 */
double kerr_disk_spectrum(const KerrBlackHole *bh, double M_dot,
                          double nu, double iota);

/**
 * Iron K-alpha line profile from a relativistic accretion disk.
 * The line is broadened and skewed by Doppler shifts and
 * gravitational redshift (Fabian+1989, Laor 1991).
 *
 * Key observable for measuring black hole spin.
 */
int kerr_iron_line_profile(const KerrBlackHole *bh, double iota,
                           double r_in, double r_out,
                           int n_bins, double *energy_bins,
                           double *flux);

/* ==========================================================================
 * L6 -- Black Hole Shadow
 * ========================================================================== */

/**
 * Compute the shadow shape and characteristic parameters.
 * The shadow is the photon capture cross-section as seen by
 * a distant observer at inclination angle iota.
 */
int kerr_shadow_observables(const KerrBlackHole *bh, double iota,
                            KerrShadow *shadow);

/**
 * Shadow radius approximation (EHT collaboration):
 * R_shadow ~ 5 M for Schwarzschild
 * For Kerr: depends on spin and inclination
 */
double kerr_shadow_radius_approx(const KerrBlackHole *bh, double iota);

/**
 * Spin measurement from shadow shape using the EHT method.
 * The asymmetry of the shadow crescent encodes the spin.
 *
 * Reference: Event Horizon Telescope Collab. (2019), ApJ 875, L1
 */
double kerr_spin_from_shadow_deviation(double deviation, double iota);

/* ==========================================================================
 * L7 -- Gravitational Wave Ringdown
 * ========================================================================== */

/**
 * Quasi-normal mode (QNM) frequency for Kerr l=m=2 mode.
 *
 * omega_QNM = omega_R + i omega_I
 *
 * The real part gives the ringdown frequency, the imaginary part
 * gives the damping time.
 *
 * Reference: Berti, Cardoso & Will (2006), PRD 73, 064030
 */
double kerr_qnm_frequency_real(const KerrBlackHole *bh, int l, int m, int n);

/** QNM damping time: tau = 1 / |omega_I| */
double kerr_qnm_damping_time(const KerrBlackHole *bh, int l, int m, int n);

/**
 * QNM frequency using the fitting formula by Berti+2006.
 * omega = [f1 + f2 (1-a)^f3] / M  (approximate)
 */
double kerr_qnm_frequency_fit(const KerrBlackHole *bh, int l, int m, int n);

/** Quality factor Q = omega_R / (2 omega_I) */
double kerr_qnm_quality_factor(const KerrBlackHole *bh, int l, int m, int n);

#endif /* KERR_PHYSICS_H */