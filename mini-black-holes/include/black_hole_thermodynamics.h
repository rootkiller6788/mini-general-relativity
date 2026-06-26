/**
 * @file black_hole_thermodynamics.h
 * @brief Black hole thermodynamics: Hawking radiation, Bekenstein-Hawking entropy,
 *        four laws of BH mechanics, evaporation.
 *
 * Reference: Hawking 1974,1975; Bekenstein 1972,1973; Wald 1984 Ch.7, Ch.14
 *
 * This is the core of modern black hole physics — the discovery that
 * black holes are thermal objects with temperature and entropy,
 * unifying GR, QFT, and thermodynamics.
 *
 * L1 Definitions: T_H, S_BH, surface gravity κ
 * L2 Core Concepts: BH thermodynamics, Hawking radiation, information paradox
 * L4 Fundamental Laws: Four laws of black hole mechanics ↔ thermodynamics
 * L6 Canonical Systems: Schwarzschild/Kerr BH thermodynamics
 */

#ifndef BLACK_HOLE_THERMODYNAMICS_H
#define BLACK_HOLE_THERMODYNAMICS_H

#include "black_hole_metrics.h"

/* ================================================================
 *  THERMODYNAMIC STATE OF A BLACK HOLE
 * ================================================================ */

/**
 * @brief Thermodynamic state of a Schwarzschild black hole.
 *
 * Bekenstein (1972): Black hole must have entropy proportional to area,
 * or else one could violate the second law by throwing entropy into it.
 *
 * Hawking (1974): Black holes radiate! T_H = ħκ/(2πk_B c).
 * For Schwarzschild: T_H = ħc³/(8πGMk_B).
 */
typedef struct {
    double mass_kg;          /**< Mass [kg] */
    double temperature;      /**< Hawking temperature T_H [K] */
    double entropy;          /**< Bekenstein-Hawking entropy S_BH [J/K] */
    double area;             /**< Horizon area A = 4πr_s² [m²] */
    double surface_gravity;  /**< κ = c⁴/(4GM) [m/s²] */
    double luminosity;       /**< Hawking radiation power P_H [W] */
    double evaporation_time; /**< Lifetime τ_evap [s] */
    double heat_capacity;    /**< C = dM/dT (negative for Schwarzschild!) [J/K] */
} SchwarzschildThermo;

/**
 * @brief Thermodynamic state of a Kerr-Newman black hole.
 *
 * General first law:
 * dM = T_H dS + Ω_H dJ + Φ_H dQ
 *
 * where S = A/(4G) (in Planck units), κ = surface gravity,
 * Ω_H = a/(r_+² + a²) = horizon angular velocity,
 * Φ_H = Q r_+/(r_+² + a²) = horizon electric potential.
 */
typedef struct {
    double mass_kg;
    double J;                /**< Angular momentum [kg·m²/s] */
    double Q;                /**< Charge [C] */
    double temperature;      /**< T_H [K] */
    double entropy;          /**< S_BH [J/K] */
    double area;             /**< Horizon area [m²] */
    double surface_gravity;  /**< κ [m/s²] */
    double omega_H;          /**< Horizon angular velocity Ω_H [rad/s] */
    double phi_H;            /**< Horizon electric potential Φ_H [V] */
    double irreducible_mass; /**< M_irr — cannot be decreased classically */
    double evaporation_time;
    double heat_capacity;
} KerrNewmanThermo;

/* ================================================================
 *  HAWKING TEMPERATURE
 * ================================================================ */

/**
 * @brief Hawking temperature for Schwarzschild BH.
 *
 * T_H = ħ c³ / (8π G M k_B)
 *
 * For a solar mass BH: T_H ≈ 6.17 × 10⁻⁸ K (colder than CMB!).
 * For a lunar mass BH (~10²² kg): T_H ≈ 1.7 K.
 * For a 10¹² kg primordial BH: T_H ≈ 1.2 × 10¹¹ K (very hot!).
 *
 * Temperature is INVERSELY proportional to mass — BH gets HOTTER
 * as it evaporates (negative heat capacity → runaway).
 *
 * Reference: Hawking, Nature 248, 30 (1974); Comm. Math. Phys. 43, 199 (1975)
 */
double hawking_temperature_schwarzschild(double mass_kg);

/**
 * @brief Hawking temperature for general (Kerr-Newman) BH.
 *
 * T_H = ħκ/(2πk_B c)
 * where κ = (r_+ - M)/(r_+² + a²) for KN, and
 * r_+ = M + sqrt(M² - a² - Q²)
 *
 * For extremal BH (M² = a² + Q²): T_H = 0 (degenerate).
 * Third law: cannot reach T_H = 0 in finite steps.
 */
double hawking_temperature_kn(double mass_kg, double a_m, double Q_m);

/* ================================================================
 *  BEKENSTEIN-HAWKING ENTROPY
 * ================================================================ */

/**
 * @brief Bekenstein-Hawking entropy.
 *
 * S_BH = k_B c³ A / (4 G ħ) = k_B A / (4 l_P²)
 *
 * where l_P = sqrt(ħG/c³) ≈ 1.616 × 10⁻³⁵ m is the Planck length,
 * and A = 4πr_s² = 16πG²M²/c⁴ is the horizon area.
 *
 * For a solar mass BH: S_BH ≈ 1.1 × 10⁷⁷ k_B — ENORMOUS.
 * A BH is the most entropic object for a given size ('t Hooft 1993).
 *
 * The entropy is proportional to AREA, not volume — this is the
 * origin of the holographic principle.
 *
 * Reference: Bekenstein, Phys. Rev. D 7, 2333 (1973);
 *            Hawking, Phys. Rev. D 13, 191 (1976)
 */
double bekenstein_hawking_entropy(double horizon_area);

/**
 * @brief Bekenstein-Hawking entropy from mass (Schwarzschild case).
 * S = 4π k_B G M² / (ħ c)
 */
double bh_entropy_from_mass(double mass_kg);

/**
 * @brief Horizon area of Schwarzschild BH: A = 4π r_s² = 16π G² M² / c⁴.
 */
double schwarzschild_horizon_area(double mass_kg);

/**
 * @brief Horizon area of Kerr BH:
 * A = 4π (r_+² + a²) = 8π M (M + sqrt(M² - a²))
 *
 * Note: area depends on M and a, but NOT on Q in KN (adds Q²):
 * A = 4π (r_+² + a²) with r_+ = M + sqrt(M² - a² - Q²).
 */
double kerr_horizon_area(double M, double a);

/**
 * @brief Horizon area of Kerr-Newman BH.
 */
double kn_horizon_area(double M, double a, double Q);

/* ================================================================
 *  BLACK HOLE EVAPORATION
 * ================================================================ */

/**
 * @brief Hawking radiation power (luminosity) for Schwarzschild BH.
 *
 * P_H ≈ ħ c⁶ / (15360 π G² M²)   (simplified, massless scalar+photon+graviton)
 *
 * The actual spectrum is a greybody — the BH is not a perfect blackbody.
 * Greybody factors suppress emission at low frequencies (λ > r_s).
 *
 * Reference: Page, Phys. Rev. D 13, 198 (1976)
 */
double hawking_luminosity(double mass_kg);

/**
 * @brief Mass loss rate: dM/dt = -P_H / c².
 */
double mass_loss_rate(double mass_kg);

/**
 * @brief Evaporation lifetime: τ_evap = ∫₀^{M₀} c² dM / P_H.
 *
 * For Schwarzschild: τ ≈ 5120 π G² M³ / (ħ c⁴)
 *                     ≈ 2.1 × 10⁶⁷ (M/M_sun)³ years.
 *
 * A stellar-mass BH lives 10⁶⁷ years (>> age of universe ~ 1.38×10¹⁰ yr).
 * But a 10¹² kg primordial BH evaporates in ~ age of universe.
 *
 * Reference: Page, Phys. Rev. D 13, 198 (1976);
 *            MacGibbon, Nature 329, 308 (1987)
 */
double evaporation_time(double mass_kg);

/**
 * @brief Temperature evolution of an evaporating BH.
 * T(t) = T₀ / (1 - t/τ_evap)^{1/3}  (diverges at t → τ_evap!)
 *
 * The "final burst": in the last ~0.1s, T skyrockets and
 * the BH releases ~10²² J in gamma rays.
 */
double evaporating_temperature(double mass_initial_kg, double time_s);

/**
 * @brief Mass evolution of an evaporating BH.
 * M(t) = M₀ (1 - t/τ_evap)^{1/3}
 */
double evaporating_mass(double mass_initial_kg, double time_s);

/* ================================================================
 *  FOUR LAWS OF BLACK HOLE MECHANICS
 * ================================================================ */

/**
 * @brief The Four Laws of Black Hole Mechanics ↔ Thermodynamics.
 *
 * | Law  | Mechanics (Bardeen-Carter-Hawking 1973) | Thermodynamics |
 * |------|---------------------------------------|----------------|
 * | 0th  | κ is constant on the horizon           | T uniform at equilibrium |
 * | 1st  | dM = (κ/8π)dA + Ω_H dJ + Φ_H dQ       | dE = T dS + μ dN + ... |
 * | 2nd  | δA ≥ 0 (area theorem)                 | ΔS ≥ 0 |
 * | 3rd  | κ → 0 cannot be reached in finite steps| T=0 unreachable |
 *
 * The identification T = ħκ/(2π) and S = A/(4G) makes the
 * analogy exact — BH thermodynamics is REAL thermodynamics,
 * not an analogy.
 */

/**
 * @brief Zeroth law: Check that surface gravity κ is uniform on horizon.
 * For stationary BH, this is an identity (κ = const on Killing horizon).
 * @return 1 if κ is uniform (always for stationary BH), 0 otherwise
 */
int zeroth_law_check(double kappa1, double kappa2);

/**
 * @brief First law of BH mechanics (Schwarzschild):
 * dM = (κ/8πG) dA = T_H dS
 *
 * Verify: dM - T_H * dS ≈ 0 for small mass changes.
 * @param M Mass
 * @param dM Small mass change
 * @return Residual |dM - T_H dS|
 */
double first_law_schwarzschild_check(double M, double dM);

/**
 * @brief First law of BH mechanics (general):
 * dM = T_H dS + Ω_H dJ + Φ_H dQ
 *
 * Checks the differential relation for a Kerr-Newman BH.
 * @return Residual magnitude
 */
double first_law_kn_check(double M, double a, double Q,
                           double dM, double dJ, double dQ);

/**
 * @brief Second law (area theorem): δA ≥ 0 for any classical process.
 *
 * Hawking's area theorem (1971): In any classical process,
 * the total area of all black hole horizons never decreases.
 *
 * This was Hawking's key insight that led Bekenstein to propose
 * BH entropy S ∝ A. The area theorem is the GR analog of ΔS ≥ 0.
 *
 * Proof relies on: (1) null energy condition, (2) cosmic censorship,
 * (3) Raychaudhuri equation → expansion θ of null geodesics cannot
 * become positive if initially zero/non-positive → trapped surfaces.
 *
 * @param A_initial Initial total horizon area
 * @param A_final Final total horizon area
 * @return A_final - A_initial (should be ≥ 0 classically)
 */
double area_increase_theorem(double A_initial, double A_final);

/**
 * @brief Third law: κ → 0 as extremality approached.
 * For a KN BH: κ = sqrt(M² - a² - Q²) / (2M² - Q² + 2M sqrt(M² - a² - Q²)).
 * As a²+Q² → M², κ → 0. Cannot reach in finite steps (Israel 1986).
 *
 * @param M Mass
 * @param a Spin
 * @param Q Charge
 * @return Surface gravity κ
 */
double surface_gravity_kn(double M, double a, double Q);

/* ================================================================
 *  BLACK HOLE HEAT CAPACITY AND STABILITY
 * ================================================================ */

/**
 * @brief Heat capacity of Schwarzschild BH: C = dM/dT.
 *
 * C = -8π² k_B G M² / (ħ c) = -(ħ c⁵) / (8π G k_B T²)
 *
 * NEGATIVE! This means the BH gets hotter as it loses energy —
 * thermodynamically unstable. A BH in a radiation bath of the
 * same temperature will either evaporate completely or grow
 * without bound (runaway).
 *
 * This is why astrophysical BHs (M > M_sun) are stable against
 * evaporation — CMB photons are hotter, so they gain mass.
 */
double schwarzschild_heat_capacity(double mass_kg);

/**
 * @brief Heat capacity of Kerr-Newman BH.
 * Can become positive for sufficiently large a or Q near extremality.
 */
double kn_heat_capacity(double M, double a, double Q);

/**
 * @brief Critical mass where T_BH = T_CMB (cosmic microwave background, ~2.725 K).
 * Below this mass, BH evaporates; above, it grows by accretion of CMB photons.
 *
 * M_crit ≈ 4.5 × 10²² kg ≈ 0.75% lunar mass.
 * Only primordial BHs with M < M_crit can have evaporated by now.
 */
double critical_mass_cmb(double T_cmb);

/**
 * @brief Compute the irreducible mass of a Kerr-Newman BH.
 *
 * M_irr² = (1/2) (M² - Q²/2 + M sqrt(M² - a² - Q²))
 *
 * M_irr CANNOT decrease in any classical process (Christodoulou 1970).
 * Fraction of extractable energy: (M - M_irr)/M ≤ 29% for extremal Kerr.
 *
 * This is the mass-energy bound to the BH that cannot be extracted
 * via Penrose process or superradiance.
 */
double irreducible_mass(double M, double a, double Q);

/**
 * @brief Maximum extractable energy from a Kerr BH via Penrose process.
 * E_extractable = M - M_irr. For extremal Kerr (a=M): up to 29% of M c².
 */
double extractable_energy_kerr(double M, double a);

/* ================================================================
 *  THERMODYNAMIC STATE FUNCTIONS
 * ================================================================ */

/**
 * @brief Compute full thermodynamic state of Schwarzschild BH.
 * Fills all fields in SchwarzschildThermo struct.
 */
void schwarzschild_thermo_compute(const SchwarzschildBH *bh,
                                   SchwarzschildThermo *thermo);

/**
 * @brief Compute thermodynamic state of Kerr-Newman BH.
 */
void kn_thermo_compute(const KerrNewmanBH *bh,
                        KerrNewmanThermo *thermo);

/**
 * @brief Print thermodynamic state to stdout.
 */
void thermo_print_schwarzschild(const SchwarzschildThermo *st);

/**
 * @brief Print Kerr-Newman thermodynamic state.
 */
void thermo_print_kn(const KerrNewmanThermo *kt);

/**
 * @brief Hawking particle emission spectrum (Planckian with greybody factor).
 *
 * d²N/(dt dω) = (1/(2π)) Γ(ω) / (exp(ħω/k_B T_H) - (-1)^{2s})
 *
 * where Γ(ω) is the greybody factor (transmission probability through
 * the BH potential barrier), s is the spin of the emitted particle.
 *
 * This function computes the spectral flux for massless scalar particles
 * (s=0, bosonic).
 *
 * @param omega Angular frequency [rad/s]
 * @param T_H Hawking temperature [K]
 * @param greybody_factor Γ(ω) ∈ [0,1]
 * @return Spectral flux d²N/(dt dω)
 */
double hawking_spectrum_scalar(double omega, double T_H,
                                double greybody_factor);

/**
 * @brief Page time: when BH has radiated half its initial entropy.
 *
 * t_Page ≈ M₀³ × (5120π G² / (ħ c⁴))
 *
 * Page (1993) showed that information begins to emerge from the BH
 * after the Page time, resolving the apparent paradox — the
 * entanglement entropy of the radiation follows the Page curve:
 * first increases, then decreases after t_Page.
 *
 * @param mass_initial_kg Initial mass
 * @return Page time [s]
 */
double page_time(double mass_initial_kg);

/**
 * @brief Page curve: entanglement entropy of Hawking radiation.
 *
 * Before Page time: S_rad(t) ≈ 4πk_B (M₀² - M(t)²) / (ħ c)
 * After Page time:  S_rad(t) ≈ 4πk_B M(t)² / (ħ c)
 *
 * At t = t_Page, S_rad = S_BH(initial)/2.
 *
 * @param mass_initial M₀ [kg]
 * @param time Current time [s]
 * @return S_rad [J/K]
 */
double page_curve_entropy(double mass_initial_kg, double time_s);

#endif /* BLACK_HOLE_THERMODYNAMICS_H */
