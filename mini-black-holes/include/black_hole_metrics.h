/**
 * @file black_hole_metrics.h
 * @brief Metric tensor definitions for Schwarzschild, Kerr, Reissner-Nordström, Kerr-Newman.
 *
 * Reference: Wald 1984, Ch. 6-7; MTW 1973, Ch. 33; Carroll 2004, Ch. 5
 *
 * This header defines the exact metric tensor components for all four
 * electrovac black hole solutions of the Einstein-Maxwell equations.
 *
 * L1 Definitions: g_μν for Schwarzschild / Kerr / RN / KN
 * L3 Math Structures: Christoffel symbols, Riemann/Ricci/Wein tensors, Kretschmann scalar
 * L6 Canonical Systems: All four BH solutions with horizon/isochrone/orbit computations
 */

#ifndef BLACK_HOLE_METRICS_H
#define BLACK_HOLE_METRICS_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 *  PHYSICAL CONSTANTS (SI)
 * ================================================================ */
#define BH_G                 6.67430e-11     /* Newton G [m^3 kg^-1 s^-2] */
#define BH_C                 2.99792458e8    /* Speed of light c [m/s] */
#define BH_HBAR              1.054571817e-34 /* Reduced Planck ħ [J·s] */
#define BH_KB                1.380649e-23    /* Boltzmann k_B [J/K] */
#define BH_SOLAR_MASS        1.98847e30      /* Solar mass [kg] */
#define BH_EPSILON0          8.8541878128e-12 /* Vacuum permittivity ε₀ [F/m] */

/* ================================================================
 *  4D METRIC TENSOR STRUCTURE
 * ================================================================ */

/**
 * @brief 4D Lorentzian metric tensor g_{μν}
 *
 * Index convention: 0=t (time), 1=r (radial), 2=θ (polar), 3=φ (azimuthal)
 * Signature: (-, +, +, +) — mostly-plus convention (spacelike)
 *
 * This is the fundamental geometric object in GR, encoding all
 * information about spacetime curvature.
 */
typedef struct {
    double g[4][4];        /**< Covariant metric components g_{μν} */
    double g_inv[4][4];    /**< Inverse contravariant metric g^{μν} */
    int signature;          /**< 0 = (-,+,+,+), 1 = (+,-,-,-) */
} Metric4D;

/* ================================================================
 *  BLACK HOLE PARAMETER STRUCTURES
 * ================================================================ */

/** @brief Schwarzschild (1916): spherical, uncharged, non-rotating */
typedef struct {
    double M;                /**< Mass in geometric units (meters) */
    double M_kg;             /**< Mass in SI (kg) */
    double r_s;              /**< Schwarzschild radius r_s = 2GM/c² */
    int geometric;           /**< Flag: 1 if using G=c=1 units */
} SchwarzschildBH;

/** @brief Kerr (1963): axisymmetric, uncharged, rotating */
typedef struct {
    double M;                /**< Mass (geometric meters) */
    double M_kg;             /**< Mass (kg) */
    double a;                /**< Spin parameter a = J/(Mc) (meters) */
    double J;                /**< Angular momentum (kg·m²/s) */
    double r_plus;           /**< Outer event horizon r_+ = M + sqrt(M²-a²) */
    double r_minus;          /**< Inner Cauchy horizon r_- = M - sqrt(M²-a²) */
    double r_ergo_eq;        /**< Equatorial ergosphere boundary = 2M */
    int geometric;
} KerrBH;

/** @brief Reissner-Nordström (1916,1918): spherical, charged, non-rotating */
typedef struct {
    double M;
    double M_kg;
    double Q;                /**< Charge in geometric units (meters) */
    double Q_coulomb;        /**< Charge in Coulombs */
    double r_plus;           /**< Outer horizon = M + sqrt(M² - Q²) */
    double r_minus;          /**< Inner horizon = M - sqrt(M² - Q²) */
    int geometric;
} ReissnerNordstromBH;

/** @brief Kerr-Newman (1965): most general electrovac BH */
typedef struct {
    double M;
    double M_kg;
    double a;                /**< Spin parameter (meters) */
    double J;                /**< Angular momentum (kg·m²/s) */
    double Q;                /**< Charge (geometric meters) */
    double Q_coulomb;
    double r_plus;           /**< Outer horizon = M + sqrt(M² - a² - Q²) */
    double r_minus;          /**< Inner Cauchy horizon */
    int geometric;
} KerrNewmanBH;

/** @brief Position in Boyer-Lindquist / Schwarzschild coordinates */
typedef struct {
    double t;                /**< Coordinate time */
    double r;                /**< Radial coordinate (areal radius) */
    double theta;            /**< Polar angle [0, π] */
    double phi;              /**< Azimuthal angle [0, 2π) */
} BLPosition;

/** @brief 4-velocity of an observer (timelike if u^μ u_μ = -1) */
typedef struct {
    double u[4];             /**< u^μ components (t, r, θ, φ) */
} FourVelocity;

/** @brief 4-momentum of a particle */
typedef struct {
    double p[4];             /**< p^μ components */
} FourMomentum;

/* ================================================================
 *  UNIT CONVERSIONS
 * ================================================================ */

/**
 * @brief Mass (kg) → geometric length (m): M_geom = GM/c².
 * For M_sun: GM_sun/c² ≈ 1477 m.
 */
double mass_to_geometric(double mass_kg);

/**
 * @brief Geometric mass (m) → SI mass (kg): M_kg = c²M_geom/G.
 */
double geometric_to_mass(double M_geom);

/**
 * @brief Angular momentum J (kg·m²/s) → spin parameter a (m): a = J/(Mc).
 */
double spin_to_geometric(double J_kg_m2_s, double M_kg);

/**
 * @brief Electric charge Q (C) → geometric charge (m).
 * Q_geom = sqrt(G/(4πε₀c⁴)) · Q. For 1C: ~ 8.61e-18 m.
 */
double charge_to_geometric(double Q_coulomb);

/**
 * @brief Geometric charge (m) → SI charge (C).
 */
double geometric_to_charge(double Q_geom);

/* ================================================================
 *  SCHWARZSCHILD METRIC FUNCTIONS
 * ================================================================ */

/**
 * @brief Initialize Schwarzschild black hole from SI mass.
 * Computes M_geom = GM/c², r_s = 2GM/c².
 */
void schwarzschild_init(SchwarzschildBH *bh, double mass_kg);

/**
 * @brief Schwarzschild radius: r_s = 2GM/c².
 * For a solar mass BH: r_s ≈ 2954 m ≈ 3 km.
 */
double schwarzschild_radius(double mass_kg);

/**
 * @brief Compute full metric tensor g_{μν} at position pos.
 *
 * g_tt = -(1 - r_s/r)
 * g_rr = (1 - r_s/r)^{-1}
 * g_θθ = r²
 * g_φφ = r² sin²θ
 * All off-diagonal = 0.
 *
 * Also computes the inverse metric g^{μν}.
 * Valid for r > r_s. The apparent singularity at r = r_s is
 * removable via coordinate transformation (Eddington-Finkelstein).
 *
 * Birkhoff (1923): Schwarzschild is the UNIQUE spherically symmetric
 * vacuum solution — any spherically symmetric gravitational field
 * in vacuum must be static and Schwarzschild.
 */
void schwarzschild_metric(const SchwarzschildBH *bh,
                          const BLPosition *pos,
                          Metric4D *metric);

/**
 * @brief Single metric component g_{μν} for Schwarzschild.
 * @param mu Index 0-3 (t,r,θ,φ)
 * @param nu Index 0-3
 */
double schwarzschild_metric_comp(const SchwarzschildBH *bh,
                                  double r, double theta,
                                  int mu, int nu);

/**
 * @brief Inverse metric component g^{μν} for Schwarzschild.
 */
double schwarzschild_inv_metric_comp(const SchwarzschildBH *bh,
                                      double r, double theta,
                                      int mu, int nu);

/* ================================================================
 *  KERR METRIC FUNCTIONS
 * ================================================================ */

/**
 * @brief Initialize Kerr BH from SI mass and angular momentum.
 * Checks cosmic censorship: |a| ≤ M, otherwise clamps a = M (extremal).
 */
void kerr_init(KerrBH *bh, double mass_kg, double J_kg_m2_s);

/**
 * @brief Compute Σ = r² + a² cos²θ (Kerr metric function).
 */
double kerr_sigma(double r, double theta, double a);

/**
 * @brief Compute Δ = r² - 2Mr + a² (Kerr horizon function).
 */
double kerr_delta(double r, double M, double a);

/**
 * @brief Outer event horizon: r_+ = M + sqrt(M² - a²).
 */
double kerr_horizon_radius(double M, double a);

/**
 * @brief Ergosphere boundary at polar angle θ.
 * r_E(θ) = M + sqrt(M² - a² cos²θ).
 * At equator (θ=π/2): r_E = 2M. At poles (θ=0): r_E = r_+.
 */
double kerr_ergosphere_radius(double M, double a, double theta);

/**
 * @brief Full Kerr metric in Boyer-Lindquist coordinates.
 *
 * Σ = r² + a²cos²θ,  Δ = r² - 2Mr + a²
 * g_tt = -(1 - 2Mr/Σ)
 * g_tφ = g_φt = -2aMr sin²θ / Σ
 * g_rr = Σ/Δ
 * g_θθ = Σ
 * g_φφ = (r² + a² + 2Ma²r sin²θ/Σ) sin²θ
 *
 * Kerr (1963): Only stationary, axisymmetric, asymptotically flat
 * vacuum solution — the rotating black hole.
 */
void kerr_metric(const KerrBH *bh,
                 const BLPosition *pos,
                 Metric4D *metric);

/**
 * @brief Single Kerr metric component.
 */
double kerr_metric_comp(const KerrBH *bh,
                         const BLPosition *pos,
                         int mu, int nu);

/**
 * @brief Frame-dragging angular velocity: ω = -g_{tφ} / g_{φφ}.
 *
 * This is the angular velocity of a Zero Angular Momentum Observer
 * (ZAMO), forced to rotate with the spacetime. At the horizon,
 * ω → ω_H = a/(2Mr_+), the horizon angular velocity.
 */
double kerr_frame_dragging_omega(const KerrBH *bh, const BLPosition *pos);

/**
 * @brief Horizon angular velocity: Ω_H = a/(2Mr_+) = ω|_{r=r_+}.
 */
double kerr_horizon_angular_velocity(const KerrBH *bh);

/**
 * @brief Surface gravity of Kerr BH: κ = sqrt(M² - a²)/(2M r_+).
 * κ → 0 as a → M (extremal limit: zero temperature).
 */
double kerr_surface_gravity(const KerrBH *bh);

/* ================================================================
 *  REISSNER-NORDSTRÖM METRIC FUNCTIONS
 * ================================================================ */

/**
 * @brief Initialize RN BH from SI mass and charge.
 */
void reissner_nordstrom_init(ReissnerNordstromBH *bh,
                              double mass_kg, double charge_coulomb);

/**
 * @brief Full RN metric.
 *
 * g_tt = -(1 - 2M/r + Q²/r²)
 * g_rr = (1 - 2M/r + Q²/r²)^{-1}
 * g_θθ = r²,  g_φφ = r² sin²θ
 *
 * Two horizons at r_± = M ± sqrt(M² - Q²).
 * |Q| ≤ M (cosmic censorship). For |Q| = M: extremal RN.
 * For |Q| > M: naked singularity (no horizon).
 */
void reissner_nordstrom_metric(const ReissnerNordstromBH *bh,
                                const BLPosition *pos,
                                Metric4D *metric);

double reissner_nordstrom_metric_comp(const ReissnerNordstromBH *bh,
                                       double r, double theta,
                                       int mu, int nu);

/**
 * @brief Electric potential at horizon: Φ_H = Q / r_+.
 */
double rn_horizon_electric_potential(const ReissnerNordstromBH *bh);

/* ================================================================
 *  KERR-NEWMAN METRIC FUNCTIONS
 * ================================================================ */

/**
 * @brief Initialize KN BH from SI mass, angular momentum, and charge.
 *
 * The most general stationary, axisymmetric, asymptotically flat
 * solution of Einstein-Maxwell equations (Newman et al., 1965).
 *
 * Characterized by only 3 parameters: M, J, Q (no-hair theorem).
 * r_± = M ± sqrt(M² - a² - Q²).
 * Requirement: a² + Q² ≤ M² (cosmic censorship).
 */
void kerr_newman_init(KerrNewmanBH *bh, double mass_kg,
                       double J_kg_m2_s, double charge_coulomb);

/**
 * @brief Full KN metric. Same form as Kerr but Δ = r² - 2Mr + a² + Q².
 */
void kerr_newman_metric(const KerrNewmanBH *bh,
                         const BLPosition *pos,
                         Metric4D *metric);

double kerr_newman_metric_comp(const KerrNewmanBH *bh,
                                const BLPosition *pos,
                                int mu, int nu);

/* ================================================================
 *  ORBITAL QUANTITIES
 * ================================================================ */

/**
 * @brief Photon sphere radius: r_ph = 3M (Schwarzschild).
 * For Kerr prograde: r_ph decreases from 3M (a=0) to M (a=M).
 */
double photon_sphere_radius(double M);

/**
 * @brief Photon orbit radius for Kerr (prograde equatorial).
 * r_ph = 2M [1 + cos(2/3 arccos(∓|a|/M))].
 * (-) for prograde, (+) for retrograde.
 */
double kerr_photon_orbit_radius(double M, double a, int prograde);

/**
 * @brief ISCO radius for Schwarzschild: r_ISCO = 6M.
 * @return 6M in geometric units.
 */
double schwarzschild_isco(double M);

/**
 * @brief ISCO radius for Kerr black hole.
 *
 * Prograde: r_ISCO decreases from 6M (a=0) to M (a=M, extremal).
 * Retrograde: r_ISCO increases from 6M to 9M.
 *
 * Using the Bardeen et al. (1972) formula:
 * Z1 = 1 + (1-a²/M²)^{1/3} [(1+a/M)^{1/3} + (1-a/M)^{1/3}]
 * Z2 = sqrt(3a²/M² + Z1²)
 * Prograde:  r_ISCO = M (3 + Z2 - sqrt((3-Z1)(3+Z1+2Z2)))
 * Retrograde: r_ISCO = M (3 + Z2 + sqrt((3-Z1)(3+Z1+2Z2)))
 */
double kerr_isco_radius(double M, double a, int prograde);

/**
 * @brief Innermost bound orbit for Kerr.
 * r_IBCO = 2M ∓ a + 2 sqrt(M² ∓ aM). (-) prograde, (+) retrograde.
 */
double kerr_ibco_radius(double M, double a, int prograde);

/**
 * @brief Keplerian orbital frequency at radius r (Schwarzschild).
 * Ω_K = sqrt(GM/r³). In geometric units: Ω_K = sqrt(M/r³).
 */
double schwarzschild_kepler_freq(double M, double r);

/* ================================================================
 *  CURVATURE QUANTITIES
 * ================================================================ */

/**
 * @brief Kretschmann scalar for Schwarzschild: K = R_{μνρσ} R^{μνρσ}.
 *
 * K = 48 M² / r⁶ (geometric units).
 * This is the simplest curvature invariant.
 * At r → 0: K → ∞ (genuine curvature singularity).
 * At r = r_s = 2M: K = 3/(4M⁴) — finite, horizon is regular.
 *
 * For rotating BHs: K = 48M² (r⁶ - 15a² r⁴ cos²θ + 15a⁴ r² cos⁴θ - a⁶ cos⁶θ) / (r²+a²cos²θ)⁶
 */
double schwarzschild_kretschmann(double M, double r);

/**
 * @brief Kretschmann scalar for Kerr metric at specified (r, θ).
 */
double kerr_kretschmann(double M, double a, double r, double theta);

/**
 * @brief Ricci scalar for vacuum solutions (Schwarzschild, Kerr).
 * Always R = 0 in vacuum. Non-zero only with matter/charge present.
 */
double vacuum_ricci_scalar(void);

/**
 * @brief Compute proper time interval dτ² = -ds² for a timelike observer.
 *
 * dτ² = -(g_tt dt² + 2g_tφ dt dφ + g_rr dr² + g_θθ dθ² + g_φφ dφ²)
 *
 * For a static observer at fixed (r,θ,φ) in Schwarzschild:
 * dτ = sqrt(1 - r_s/r) dt — gravitational time dilation.
 */
double proper_time_interval(const Metric4D *metric,
                             double dt, double dr,
                             double dtheta, double dphi);

/**
 * @brief Proper radial distance ∫ dr/√(1 - r_s/r) from r1 to r2.
 *
 * Closed form:
 * ℓ = [r√(1 - r_s/r) + (r_s/2) ln((√(r/r_s - 1) + √(r/r_s)) / (√(r/r_s - 1) - √(r/r_s)))]_{r1}^{r2}
 *   = [r√(1 - r_s/r) + r_s·arctanh(√(r_s/r))]_{r1}^{r2}
 */
double schwarzschild_proper_distance(double r_s, double r1, double r2);

/**
 * @brief Gravitational redshift: ν_obs/ν_em = √(g_tt(emitter)/g_tt(observer)).
 * For Schwarzschild with observer at infinity: z = 1/√(1 - r_s/r) - 1.
 */
double schwarzschild_redshift(double r_s, double r_emit);

/**
 * @brief Time dilation factor dτ/dt = √(-g_tt) for a static observer.
 */
double schwarzschild_time_dilation(double r_s, double r);

/**
 * @brief Print metric tensor to stdout in a formatted table.
 */
void metric_print(const Metric4D *metric, const char *label);

/**
 * @brief Validate that the metric satisfies g_{μν} g^{νρ} = δ_μ^ρ.
 * @return Max absolute deviation from identity (should be < 1e-10).
 */
double metric_inverse_check(const Metric4D *metric);

#endif /* BLACK_HOLE_METRICS_H */
