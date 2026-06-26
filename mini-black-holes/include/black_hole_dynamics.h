/**
 * @file black_hole_dynamics.h
 * @brief Geodesic motion, orbital dynamics, Penrose process, superradiance.
 *
 * Reference: Bardeen, Press, Teukolsky 1972; MTW Ch.33; Wald Ch.6-7
 *
 * This module implements particle and photon trajectories around
 * black holes, covering both timelike and null geodesics.
 *
 * L1 Definitions: Geodesic equation, Carter constant, ISCO, IBCO
 * L3 Math Structures: Hamilton-Jacobi separability, Killing vectors
 * L5 Computational Methods: RK4 geodesic integration, root-finding
 * L6 Canonical Systems: Schwarzschild orbits, Kerr equatorial orbits
 */

#ifndef BLACK_HOLE_DYNAMICS_H
#define BLACK_HOLE_DYNAMICS_H

#include "black_hole_metrics.h"

/* ================================================================
 *  GEODESIC EQUATION
 * ================================================================ */

/**
 * @brief Geodesic state: position and 4-velocity of a test particle.
 *
 * The geodesic equation:
 * d²x^μ/dλ² + Γ^μ_{νρ} (dx^ν/dλ)(dx^ρ/dλ) = 0
 *
 * where λ is an affine parameter: proper time τ for timelike,
 * or any affine parameter for null geodesics.
 */
typedef struct {
    double x[4];            /**< Position: t, r, θ, φ */
    double v[4];            /**< 4-velocity: dx^μ/dλ */
    double affine;          /**< Affine parameter λ */
} GeodesicState;

/**
 * @brief Constants of motion for Kerr geodesics.
 *
 * Kerr spacetime admits four integrals of motion:
 * - μ = -p_μ p^μ = m² (particle mass, 0 for photons)
 * - E = -p_t (energy, conserved via ∂_t Killing vector)
 * - L_z = p_φ (axial angular momentum, conserved via ∂_φ)
 * - Q = Carter constant (1968): separability of Hamilton-Jacobi equation
 *
 * Carter (1968): Kerr metric admits a Killing-Yano tensor K_{μν},
 * giving the Carter constant Q = K_{μν} p^μ p^ν.
 * This makes the geodesic equation completely integrable.
 */
typedef struct {
    double E;               /**< Energy at infinity */
    double L_z;             /**< Axial angular momentum */
    double Q;               /**< Carter constant */
    double mu;              /**< Rest mass (0 for photons) */
} GeodesicConstants;

/* ================================================================
 *  EFFECTIVE POTENTIALS
 * ================================================================ */

/**
 * @brief Effective potential for radial motion in Schwarzschild.
 *
 * For timelike geodesics (massive particles):
 * V_eff(r) = sqrt((1 - r_s/r)(1 + L²/r²))
 *
 * where L = L_z/E is the angular momentum per unit energy.
 *
 * Orbits exist where E ≥ V_eff(r). The extrema of V_eff give
 * circular orbits.
 *
 * @param r_s Schwarzschild radius
 * @param r Radial coordinate
 * @param L_sq L² = (L_z/E)²
 * @return V_eff
 */
double schwarzschild_effective_potential(double r_s, double r, double L_sq);

/**
 * @brief Derivative of Schwarzschild effective potential dV_eff/dr.
 * Used for finding circular orbits via root-finding.
 */
double schwarzschild_eff_potential_deriv(double r_s, double r, double L_sq);

/**
 * @brief Effective potential for equatorial (θ=π/2) Kerr geodesics.
 *
 * For equatorial orbits, Q=0 and the radial equation becomes:
 * (dr/dτ)² = (E² - μ²) + 2M(μ² r + (E - aL_z/r)²)/r - (L_z² - a²(E²-μ²))/r²
 *
 * @param M Mass
 * @param a Spin
 * @param r Radius
 * @param E Energy
 * @param L_z Angular momentum
 * @param mu Rest mass (1 for massive, 0 for photon)
 * @return V_eff²
 */
double kerr_effective_potential_sq(double M, double a, double r,
                                    double E, double L_z, double mu);

/* ================================================================
 *  CIRCULAR ORBITS
 * ================================================================ */

/**
 * @brief Specific angular momentum L/E for circular Schwarzschild orbit at r.
 *
 * L_circ = r sqrt(r_s / (2r - 3r_s))
 *
 * For r → 3r_s/2 (photon sphere): L → ∞ (no massive circular orbit possible).
 * For r → ∞: L ≈ sqrt(M r) (Newtonian limit).
 *
 * Stable orbits exist for r > 6M. Between 3M and 6M: unstable circular orbits.
 * At r=6M (ISCO): L = 2√3 M, the marginally stable orbit.
 */
double schwarzschild_circular_L(double r_s, double r);

/**
 * @brief Energy per unit mass for circular Schwarzschild orbit at r.
 *
 * E_circ = (r - r_s) / sqrt(r(r - 1.5 r_s))
 *
 * At ISCO (r=6M): E = 2√2/3 ≈ 0.943.
 * Binding energy: 1 - E ≈ 5.72% of rest mass released.
 */
double schwarzschild_circular_E(double r_s, double r);

/**
 * @brief Angular velocity Ω = dφ/dt for circular orbit in Schwarzschild.
 * Ω = sqrt(M/r³) = c sqrt(r_s/(2r³)) — same form as Kepler's third law!
 */
double schwarzschild_circular_omega(double r_s, double r);

/**
 * @brief Find the radius of a circular orbit given L in Schwarzschild.
 * Solves dV_eff/dr = 0 using Newton-Raphson.
 *
 * @param r_s Schwarzschild radius
 * @param L_sq L²
 * @param r_guess Initial guess for radius
 * @return Circular orbit radius (negative if no solution)
 */
double schwarzschild_find_circular_orbit(double r_s, double L_sq, double r_guess);

/**
 * @brief Radial epicyclic frequency ω_r for slightly perturbed circular orbit.
 *
 * ω_r² = Ω² (1 - 6M/r) = (M/r³)(1 - 6M/r)
 *
 * At ISCO (r=6M): ω_r = 0 → marginal stability.
 * For r > 6M: ω_r > 0 → stable (oscillates).
 * For r < 6M: ω_r² < 0 → unstable (exponential divergence).
 */
double schwarzschild_epicyclic_freq(double r_s, double r);

/**
 * @brief Vertical epicyclic frequency ω_θ for Kerr.
 * ω_θ² = Ω² (1 - 4aΩ + 3a²/r²)
 *
 * In Kerr, ω_r ≠ ω_θ ≠ Ω, giving three distinct frequencies
 * observable in QPOs (quasi-periodic oscillations) from BH X-ray binaries.
 */
double kerr_vertical_epicyclic_freq(double M, double a, double r);

/**
 * @brief Radial epicyclic frequency for Kerr.
 * ω_r² = Ω² (1 - 6M/r + 8a√(M/r³) - 3a²/r²)
 */
double kerr_radial_epicyclic_freq(double M, double a, double r);

/* ================================================================
 *  GEODESIC INTEGRATION (RK4)
 * ================================================================ */

/**
 * @brief Compute geodesic right-hand side: dv^μ/dλ = -Γ^μ_{νρ} v^ν v^ρ.
 *
 * For given metric, position, and velocity, computes the acceleration
 * due to spacetime curvature.
 *
 * @param metric_func Function to evaluate metric at a point
 * @param bh Black hole parameters
 * @param state Geodesic state (pos + vel)
 * @param acc Output: dv^μ/dλ
 * @param h Finite difference step for Christoffel symbols
 */
void geodesic_rhs(void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
                  const void *bh,
                  const GeodesicState *state,
                  double acc[4],
                  double h);

/**
 * @brief Single RK4 step for geodesic integration.
 *
 * Advances the geodesic state by one step dλ using the
 * classic 4th-order Runge-Kutta method.
 *
 * @return 0 on success, -1 if particle enters horizon (r ≤ r_s + ε)
 */
int geodesic_rk4_step(
    void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
    const void *bh,
    GeodesicState *state,
    double dlambda,
    double h);

/**
 * @brief Integrate a geodesic from start to end affine parameter.
 *
 * Stops if particle reaches horizon or escapes to large r.
 * Stores the trajectory in pre-allocated arrays.
 *
 * @param metric_func Metric function
 * @param bh BH parameters
 * @param state Initial state (modified in place)
 * @param lambda_max Maximum affine parameter
 * @param dlambda Step size
 * @param traj_r Output: radial coordinate history
 * @param traj_phi Output: azimuthal angle history
 * @param max_steps Maximum array size
 * @param n_steps Output: actual number of steps taken
 * @return 0 if completed, 1 if hit horizon, 2 if escaped
 */
int geodesic_integrate(
    void (*metric_func)(const void *bh, const BLPosition *p, Metric4D *m),
    const void *bh,
    GeodesicState *state,
    double lambda_max, double dlambda,
    double *traj_r, double *traj_phi,
    int max_steps, int *n_steps);

/* ================================================================
 *  PENROSE PROCESS
 * ================================================================ */

/**
 * @brief Maximum energy extraction efficiency of the Penrose process.
 *
 * Penrose (1969): A particle falling into a Kerr BH can split in the
 * ergosphere (r_+ < r < r_ergo), with one fragment falling into the
 * BH with NEGATIVE energy (as measured at infinity), while the other
 * escapes with MORE energy than the original particle.
 *
 * Maximum efficiency (Bardeen et al. 1972):
 * η_max = (E_out - E_in) / E_in = (1/√2 - 1) ≈ 20.7% for extremal Kerr.
 *
 * With subsequent interactions: up to 29% of BH mass is extractable
 * (equivalent to M - M_irr).
 *
 * @param a Spin parameter
 * @param M Mass
 * @return Maximum efficiency η ∈ [0, 0.207]
 */
double penrose_process_efficiency(double M, double a);

/**
 * @brief Check if a position is inside the ergosphere (energy can be negative).
 *
 * Ergosphere: r_+ < r < r_ergo(θ).
 * Within this region, the Killing vector ∂_t becomes spacelike,
 * meaning particle energy E = -p_μ ξ^μ can be negative.
 *
 * @return 1 if inside ergosphere, 0 otherwise
 */
int is_inside_ergosphere(const KerrBH *bh, const BLPosition *pos);

/**
 * @brief Energy at infinity for a particle at given position and 4-velocity.
 *
 * E = -p_μ ξ^μ = -g_{μν} ξ^μ p^ν = -(g_{tt} p^t + g_{tφ} p^φ)
 *
 * where ξ^μ = ∂_t = (1,0,0,0) is the timelike Killing vector.
 *
 * @param metric Metric at particle position
 * @param velocity 4-velocity u^μ
 * @param rest_mass Particle rest mass
 * @return E
 */
double particle_energy_at_infinity(const Metric4D *metric,
                                    const FourVelocity *velocity,
                                    double rest_mass);

/**
 * @brief Angular momentum of a particle: L_z = p_μ η^μ = g_{φν} p^ν.
 * where η^μ = ∂_φ = (0,0,0,1) is the axial Killing vector.
 */
double particle_angular_momentum(const Metric4D *metric,
                                  const FourVelocity *velocity,
                                  double rest_mass);

/* ================================================================
 *  SUPERADIANCE (Zeldovich 1971, Starobinsky 1973)
 * ================================================================ */

/**
 * @brief Superradiance condition for scalar waves incident on Kerr BH.
 *
 * An incident wave of frequency ω and azimuthal number m extracts
 * rotational energy from the BH if:
 * ω < m Ω_H
 *
 * where Ω_H = a/(r_+² + a²) is the horizon angular velocity.
 *
 * The reflected wave is AMPLIFIED — rotational energy of the BH
 * is converted into wave energy. This is the wave analog of the
 * Penrose process.
 *
 * @param omega Wave frequency
 * @param m Azimuthal quantum number
 * @param M BH mass
 * @param a Spin parameter
 * @return 1 if superradiant, 0 otherwise
 */
int superradiance_condition(double omega, int m, double M, double a);

/**
 * @brief Maximum superradiant amplification factor for scalar waves.
 *
 * For the l=m=1 mode, max amplification ~ 4.4% for extremal Kerr (a=M).
 * For gravitational waves (s=2): up to 138% amplification possible
 * (Teukolsky & Press 1974).
 */
double superradiance_max_amplification(double M, double a, int l, int m);

/* ================================================================
 *  MASS AND SPIN EVOLUTION
 * ================================================================ */

/**
 * @brief Bardeen's accretion model: spin evolution due to accretion.
 *
 * Thin disk accretion spins up the BH. Bardeen (1970) showed that
 * a BH accreting from a thin disk reaches a/M ≈ 0.998 (Thorne limit),
 * not a/M = 1, because photons emitted by the disk are captured
 * by the BH at high spin, exerting a spin-down torque.
 *
 * @param M_initial Initial mass (kg)
 * @param a_initial Initial spin parameter (m)
 * @param accreted_mass Mass accreted (kg)
 * @param a_final Output: final spin parameter (m)
 * @return Final mass (kg)
 */
double bardeen_spin_evolution(double M_initial, double a_initial,
                               double accreted_mass, double *a_final);

/* ================================================================
 *  PHOTON TRAJECTORIES AND SHADOW
 * ================================================================ */

/**
 * @brief Critical impact parameter for photon capture by Schwarzschild BH.
 *
 * b_crit = 3√3 M ≈ 5.196 M
 *
 * Photons with impact parameter b < b_crit are captured.
 * Photons with b = b_crit orbit at r = 3M (photon sphere: unstable).
 * Photons with b > b_crit escape to infinity after deflection.
 *
 * This defines the black hole shadow: a circle of radius b_crit
 * in the observer's sky.
 */
double schwarzschild_critical_impact(void);

/**
 * @brief Black hole shadow radius for Kerr (as seen by distant observer).
 *
 * The shadow is NOT circular for a ≠ 0 — it's displaced and D-shaped.
 * The shadow size is characterized by:
 * - R_shadow ≈ 2√3 M for Schwarzschild
 * - For Kerr, it's parametrized by celestial coordinates (α, β).
 *
 * @param M Mass
 * @param a Spin
 * @param theta_o Observer inclination angle
 * @param R_mean Output: mean shadow radius
 * @param asymmetry Output: shadow asymmetry parameter
 */
void kerr_shadow_radius(double M, double a, double theta_o,
                         double *R_mean, double *asymmetry);

/**
 * @brief Deflection angle of a photon passing near a Schwarzschild BH.
 *
 * In the weak-field limit (b >> M):
 * Δφ ≈ 4M/b (Einstein 1915 prediction for solar deflection).
 *
 * Strong field: need to integrate the null geodesic equation.
 * This function uses the exact integral:
 * Δφ = 2 ∫_{r_min}^{∞} dr / (r² sqrt(1/b² - (1-r_s/r)/r²)) - π
 *
 * @param b Impact parameter
 * @param r_s Schwarzschild radius
 * @return Deflection angle [rad]
 */
double schwarzschild_deflection_angle(double b, double r_s);

#endif /* BLACK_HOLE_DYNAMICS_H */
