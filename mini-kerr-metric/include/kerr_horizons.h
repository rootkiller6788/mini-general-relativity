/**
 * @file kerr_horizons.h
 * @brief Kerr black hole horizons, ergosphere, and ring singularity
 *
 * The Kerr geometry has:
 * - Two event horizons at r_plus/minus = M +/- sqrt(M^2 - a^2)
 * - An ergosphere bounded by r_ergo(theta) = M + sqrt(M^2 - a^2 cos^2(theta))
 * - A ring singularity at r=0, theta=pi/2 of radius a in Kerr-Schild coords
 *
 * Reference: Wald Ch.7.5, MTW Ch.33, Chandrasekhar Ch.6
 * MIT 8.962 / Cambridge Part III GR / Princeton PHY 535
 */

#ifndef KERR_HORIZONS_H
#define KERR_HORIZONS_H

#include "kerr_metric.h"

/* ==========================================================================
 * L1 -- Horizon Structures
 * ========================================================================== */

/** Kerr horizon geometry data */
typedef struct {
    double r_plus;       /* Outer event horizon radius */
    double r_minus;      /* Inner (Cauchy) horizon radius */
    double surface_area; /* Area = 4pi(r_plus^2 + a^2) */
    double kappa;        /* Surface gravity at outer horizon */
    double omega_H;      /* Horizon angular velocity */
    double temperature;  /* Hawking temperature T = kappa hbar / (2pi k_B c) */
    double entropy;      /* Bekenstein-Hawking entropy S = A k_B c^3 / (4 G hbar) */
    double irreducible_mass; /* M_irr = sqrt(A/16pi) */
} KerrHorizonData;

/** Ergosphere boundary data at a given polar angle */
typedef struct {
    double r_ergo;       /* Ergosphere radius r_ergo(theta) */
    int inside;          /* 1 if point is inside ergosphere */
    double g_tt;         /* g_tt component (positive inside ergosphere) */
} KerrErgosphereData;

/** Ring singularity characterization */
typedef struct {
    double ring_radius;  /* a in BL coordinates = ring radius in KS */
    double r_bl;         /* r=0 in BL (entire disk r=0) */
    double theta_bl;     /* theta=pi/2 for the singular ring */
    int is_singular;     /* 1 if point is on the singular ring */
} KerrRingSingularity;

/** Black hole type classification */
typedef enum {
    BH_SCHWARZSCHILD = 0,  /* a = 0, non-rotating */
    BH_KERR_SLOW     = 1,  /* 0 < a < M, generic rotating */
    BH_KERR_EXTREMAL = 2,  /* a = M, maximal spin */
    BH_NAKED_SING    = 3,  /* a > M, cosmic censorship violation */
    BH_INVALID       = -1  /* Negative mass or other invalid params */
} KerrBHType;

/* ==========================================================================
 * L2 -- Horizon Properties
 * ========================================================================== */

/**
 * Compute both horizon radii: r_plus = M + sqrt(M^2 - a^2),
 * r_minus = M - sqrt(M^2 - a^2).
 * For a=M (extremal): r_plus = r_minus = M.
 * For a>M (naked): returns NaN.
 *
 * Complexity: O(1)
 * Theorem: Roots of Delta(r) (Wald Thm 7.1)
 */
int kerr_horizon_radii(const KerrBlackHole *bh, double *r_plus,
                       double *r_minus);

/** Classify the black hole type based on M and a */
KerrBHType kerr_bh_classify(const KerrBlackHole *bh);

/**
 * Outer horizon surface area: A = 4pi(r_plus^2 + a^2).
 * The area of a Kerr black hole is always less than Schwarzschild
 * with the same mass: A_Kerr <= 16pi M^2.
 */
double kerr_horizon_area(const KerrBlackHole *bh);

/**
 * Horizon angular velocity Omega_H = a / (r_plus^2 + a^2).
 * This is the angular velocity of the horizon as measured from infinity.
 * For Schwarzschild: Omega_H = 0.
 */
double kerr_horizon_angular_velocity(const KerrBlackHole *bh);

/**
 * Surface gravity at outer horizon:
 * kappa = (r_plus - M) / (r_plus^2 + a^2) = sqrt(M^2 - a^2) / (2 M r_plus)
 * For extremal Kerr (a=M): kappa = 0 (zero temperature).
 */
double kerr_surface_gravity(const KerrBlackHole *bh);

/**
 * Bekenstein-Hawking entropy: S_BH = k_B A / (4 l_Planck^2).
 * In geometric units (G=c=1): S_BH = A/4 = pi(r_plus^2 + a^2).
 * Returns entropy in dimensionless geometric units.
 */
double kerr_entropy_geometric(const KerrBlackHole *bh);

/**
 * Hawking temperature in geometric units:
 * T_H = kappa / (2 pi) = sqrt(M^2 - a^2) / (4 pi M r_plus).
 * For extremal Kerr: T_H = 0.
 */
double kerr_temperature_geometric(const KerrBlackHole *bh);

/**
 * Irreducible mass: M_irr^2 = A / (16 pi) = (r_plus^2 + a^2)/4.
 * M_irr is the minimum mass achievable by extracting rotational energy.
 * M^2 = M_irr^2 + J^2/(4 M_irr^2).
 */
double kerr_irreducible_mass(const KerrBlackHole *bh);

/**
 * Rotational energy available for extraction:
 * E_rot = M - M_irr = M - sqrt((r_plus^2 + a^2)/4).
 * For extremal Kerr: up to ~29.3% of M is extractable.
 */
double kerr_rotational_energy(const KerrBlackHole *bh);

/** Fill complete horizon data structure */
int kerr_horizon_data(const KerrBlackHole *bh, KerrHorizonData *data);

/* ==========================================================================
 * L3 -- Ergosphere
 * ========================================================================== */

/**
 * Ergosphere boundary: r_ergo(theta) = M + sqrt(M^2 - a^2 cos^2(theta)).
 * Between r_plus and r_ergo lies the ergosphere where no observer
 * can remain static (must rotate with the black hole).
 */
double kerr_ergosphere_radius(const KerrBlackHole *bh, double theta);

/**
 * Check if a point is inside the ergosphere.
 * Condition: r_plus < r < r_ergo(theta).
 */
int kerr_is_inside_ergosphere(const KerrBlackHole *bh,
                              const KerrBLPoint *pt);

/**
 * Ergosphere boundary data at a point.
 * Returns r_ergo, whether the point is inside, and g_tt sign.
 */
int kerr_ergosphere_data(const KerrBlackHole *bh, const KerrBLPoint *pt,
                         KerrErgosphereData *data);

/**
 * Maximum extent of the ergosphere (at equator, theta=pi/2):
 * r_ergo_max = M + sqrt(M^2) = 2M for a <= M.
 * This is independent of spin for a <= M (always 2M at equator).
 */
double kerr_ergosphere_equatorial_radius(const KerrBlackHole *bh);

/* ==========================================================================
 * L4 -- Ring Singularity
 * ========================================================================== */

/**
 * The Kerr singularity is a ring of coordinate radius a in the
 * equatorial plane (theta = pi/2) at r = 0 (Boyer-Lindquist).
 *
 * In Kerr-Schild Cartesian coordinates, r=0 corresponds to the disk
 * x^2 + y^2 <= a^2, z = 0. The ring x^2 + y^2 = a^2 is singular.
 */
int kerr_ring_singularity_data(const KerrBlackHole *bh,
                               KerrRingSingularity *ring);

/**
 * Compute the Kretschmann scalar near the ring singularity.
 * Returns infinity (DBL_MAX) if too close to the ring.
 * Safe distance threshold: 1e-10 in geometric units.
 */
double kerr_kretschmann_safe(const KerrBlackHole *bh, const KerrBLPoint *pt,
                             double safety_threshold);

/**
 * Check if a point is in the region of closed timelike curves (CTCs).
 * In Kerr, CTCs exist for r < 0 near the ring singularity where
 * g_phiphi can become negative.
 *
 * Returns: 1 if CTCs possible, 0 if chronology is safe.
 */
int kerr_has_ctc_region(const KerrBlackHole *bh, const KerrBLPoint *pt);

/* ==========================================================================
 * L5 -- BH Thermodynamics (Kerr-specific)
 * ========================================================================== */

/**
 * First law of BH thermodynamics for Kerr:
 * dM = (kappa/8pi) dA + Omega_H dJ
 *
 * Verify by numerical differentiation at current state.
 * Returns residual |dM - (kappa/8pi)dA - Omega_H dJ|.
 */
double kerr_smarr_formula_check(const KerrBlackHole *bh);

/**
 * Smarr formula for Kerr: M = 2 Omega_H J + (kappa/4pi) A.
 * Returns residual |M - 2 Omega_H J - (kappa/4pi) A|.
 */
double kerr_smarr_formula_residual(const KerrBlackHole *bh);

/**
 * Second law check: For any two Kerr BHs with M1,a1 and M2,a2 where
 * M_total >= M1+M2, check if area increases.
 * Returns area increase ratio A_total / (A1 + A2).
 */
double kerr_area_increase_ratio(double M1, double a1,
                                double M2, double a2,
                                double M_total);

/** Entropy bounds: S_Kerr(M,a) <= S_Schwarzschild(M) always holds */
double kerr_entropy_bound_ratio(const KerrBlackHole *bh);

/** Extremality bound: a/M <= 1. Returns a/M (should be <= 1). */
double kerr_extremality_parameter(const KerrBlackHole *bh);

#endif /* KERR_HORIZONS_H */