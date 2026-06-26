/**
 * @file black_hole_advanced.h
 * @brief Advanced black hole topics: area theorem, cosmic censorship,
 *        holographic principle, information paradox, ER=EPR, firewall.
 *
 * Reference: Hawking & Ellis 1973; 't Hooft 1993; Susskind 1995;
 *            Maldacena 1997; Almheiri, Marolf, Polchinski, Sully 2012;
 *            Maldacena & Susskind 2013
 *
 * L2 Core Concepts: Information paradox, holographic principle
 * L4 Fundamental Laws: Area theorem proof elements
 * L8 Advanced Topics: AdS/CFT, ER=EPR, firewalls, Page curve
 * L9 Research Frontiers: Quantum gravity implications for BHs
 */

#ifndef BLACK_HOLE_ADVANCED_H
#define BLACK_HOLE_ADVANCED_H

#include "black_hole_metrics.h"
#include "black_hole_thermodynamics.h"
#include <stdint.h>

/* ================================================================
 *  NO-HAIR THEOREM
 * ================================================================ */

/**
 * @brief No-hair theorem: a stationary BH is completely characterized
 * by only 3 externally observable parameters:
 *
 * - Mass M (measured via Kepler orbits at infinity)
 * - Angular momentum J (measured via frame-dragging/Lense-Thirring)
 * - Electric charge Q (measured via Gauss flux at infinity)
 *
 * Israel (1967, 1968): Static vacuum BH must be Schwarzschild.
 * Carter (1971): Stationary vacuum BH must be Kerr.
 * Robinson (1975): Stationary electrovac BH must be Kerr-Newman.
 *
 * This means ALL other information about the matter that formed
 * the BH is lost to the exterior — the origin of the "hair" metaphor
 * (Wheeler: "A black hole has no hair").
 */

/**
 * @brief Multipole moments of a Kerr BH.
 *
 * All multipole moments are determined by M and J alone:
 * M_l + i S_l = M (i a)^l
 *
 * Mass quadrupole: M_2 = -M a²
 * This "no-hair" relation is testable with GW observations!
 * Any deviation would signal new physics.
 */
typedef struct {
    double M;               /**< Mass monopole */
    double J;               /**< Angular momentum dipole S_1 = J */
    double M2;              /**< Mass quadrupole = -M a² */
    double S3;              /**< Current octupole = -J a² */
    double M4;              /**< Mass hexadecapole = M a⁴ */
} KerrMultipoles;

/**
 * @brief Compute Kerr multipole moments up to l=4.
 */
void kerr_multipoles(double M, double a, KerrMultipoles *mp);

/**
 * @brief Test the no-hair theorem: compare multipole ratio.
 * For Kerr: M_2 / (J²/M) = -1.
 * For a non-Kerr object (e.g., neutron star), this ratio differs.
 * LIGO-Virgo uses this to test GR (e.g., Isi et al. 2019 for GW150914).
 *
 * @param M2 Measured mass quadrupole
 * @param M Mass
 * @param J Angular momentum
 * @return M_2 / (J²/M) — should be -1 for Kerr
 */
double no_hair_test_quadrupole(double M2, double M, double J);

/**
 * @brief Check if given (M, J, Q) satisfies cosmic censorship.
 * @return 1 if a²+Q² ≤ M² (BH exists), 0 if naked singularity
 */
int cosmic_censorship_check(double M, double a, double Q);

/* ================================================================
 *  AREA THEOREM (Hawking 1971)
 * ================================================================ */

/**
 * @brief The area theorem: For any classical process involving black holes,
 * the total area of all BH horizons never decreases: δA_total ≥ 0.
 *
 * This is a theorem in classical GR (no quantum effects).
 *
 * Proof structure (Hawking & Ellis 1973, §9.2):
 * 1. Null energy condition: T_{μν} l^μ l^ν ≥ 0 for any null l^μ
 * 2. Null geodesic congruence expansion θ satisfies Raychaudhuri:
 *    dθ/dλ = -θ²/2 - σ_{μν}σ^{μν} + ω_{μν}ω^{μν} - R_{μν}l^μ l^ν ≤ 0
 * 3. θ ≤ 0 initially on trapped surfaces
 * 4. If θ < 0, geodesics converge → trapped surface → singularity (Penrose 1965)
 * 5. Horizon generators have θ ≥ 0 (otherwise horizon would enter trapped region)
 * 6. Therefore area A = ∫ θ dA dλ → dA/dλ ≥ 0
 *
 * Corollary: In BH mergers, A_final ≥ A₁ + A₂.
 * For GW150914: A_initial ≈ 2.4×10⁵ km², A_final ≈ 3.6×10⁵ km² → ΔA > 0 ✓
 */

/**
 * @brief Check area increase for a 2-BH merger: A_final ≥ A1 + A2.
 * @param M1, a1, Q1 Parameters of first BH
 * @param M2, a2, Q2 Parameters of second BH
 * @param M_final, a_final, Q_final Parameters of final BH
 * @return 1 if area theorem holds, 0 if violated
 */
int area_theorem_merger_check(double M1, double a1, double Q1,
                               double M2, double a2, double Q2,
                               double M_final, double a_final, double Q_final);

/**
 * @brief Maximum mass that can be extracted via BH merger (from area theorem).
 *
 * For two Schwarzschild BHs of masses m₁, m₂:
 * A_final ≥ A₁ + A₂ → M_final ≥ sqrt(m₁² + m₂²)
 *
 * Maximum radiated energy: E_max = (m₁+m₂) - sqrt(m₁²+m₂²)
 *
 * For equal masses m₁=m₂=m: E_max = 2m - m√2 ≈ 0.586 m c² ≈ 29% of total.
 * This is the Hawking bound on GW energy emission.
 */
double hawking_area_bound_energy(double m1, double m2);

/* ================================================================
 *  HOLOGRAPHIC PRINCIPLE
 * ================================================================ */

/**
 * @brief Bekenstein bound: maximum entropy in a region of size R
 * containing energy E: S ≤ 2π k_B R E / (ħ c).
 *
 * @param R Region radius [m]
 * @param E Energy contained [J]
 * @return Maximum entropy [J/K]
 */
double bekenstein_bound(double R, double E);

/**
 * @brief Holographic bound ('t Hooft 1993, Susskind 1995):
 * The maximum information content of a region of space is
 * proportional to its surface AREA (not volume), measured in Planck units:
 *
 * N_bits ≤ A / (4 l_P² ln 2)
 *
 * where l_P = sqrt(ħG/c³) ≈ 1.616 × 10⁻³⁵ m.
 *
 * This is the most radical implication of BH thermodynamics:
 * our 3D world may be a holographic projection of a 2D boundary theory.
 *
 * @param area_surface Boundary area [m²]
 * @return Maximum number of bits
 */
double holographic_bound_bits(double area_surface);

/**
 * @brief Planck area: l_P² = ħG/c³ ≈ 2.612 × 10⁻⁷⁰ m².
 * One bit of information requires 4 ln(2) l_P² ≈ 2.77 l_P² of area.
 *
 * For a 1 cm² surface: ~10⁶⁵ bits — a 1 TB hard drive needs only ~10⁻⁵⁴ cm².
 */
double planck_area(void);

/**
 * @brief Spherical bound: maximum information in a spherical region
 * of radius R is S ≤ A/(4 l_P²) = π R² / l_P² (in Planck units).
 *
 * Compare to naive volumetric bound for a quantum system:
 * V/l_P³ >> A/l_P² for R >> l_P → holography is much more restrictive!
 */
double spherical_holographic_bound(double radius_m);

/**
 * @brief Check if a physical system satisfies the holographic bound.
 *
 * Given a system with entropy S and bounding area A:
 * Returns A/(4 l_P² k_B) - S. Positive → bound satisfied.
 */
double holographic_bound_margin(double entropy, double area_m2);

/* ================================================================
 *  BLACK HOLE INFORMATION PARADOX
 * ================================================================ */

/**
 * @brief Information paradox states (Hawking 1976).
 *
 * The paradox: In semiclassical gravity, BH formation and evaporation
 * appears to evolve a pure quantum state into a mixed thermal state,
 * violating unitarity of quantum mechanics.
 *
 * Key tension:
 * - QM: information is conserved (unitarity) → pure → pure
 * - GR + QFT on curved spacetime: Hawking radiation is exactly thermal
 *   → pure → mixed = INFORMATION LOSS
 *
 * Resolution proposals:
 * 1. Information is lost (Hawking's original position, conceded 2004)
 * 2. Information is encoded in Hawking radiation correlations
 *    (Page 1993: information begins to emerge after Page time)
 * 3. Remnant scenario: evaporation stops at Planck mass, storing info
 * 4. AdS/CFT: Unitary CFT → information preserved (Maldacena 1997)
 * 5. Firewall: AMPS (2012) — horizon is not a smooth surface
 * 6. ER=EPR: Maldacena & Susskind (2013) — entanglement = wormhole
 */

/**
 * @brief Page information: the amount of information (in bits)
 * that has emerged in Hawking radiation by time t.
 *
 * I(t) ≈ S_initial - S_BH(t)  (before Page time)
 * I(t) ≈ S_BH(t)              (after Page time)
 *
 * Before t_Page: almost no information emerges (Hawking radiation
 * is nearly thermal). After t_Page: each emitted quantum carries
 * ~1 bit of information, purifying the earlier radiation.
 */
double page_information(double mass_initial_kg, double time_s);

/**
 * @brief Check if the Page time has been reached.
 * @return 1 if t ≥ t_Page, 0 otherwise
 */
int has_reached_page_time(double mass_initial_kg, double time_s);

/**
 * @brief Entanglement entropy between BH and its radiation.
 * Follows the Page curve (inverted V-shape):
 * S_ent(t) = min(S_BH(t), S_rad(t))
 * Peaks at t = t_Page with S_ent = S_BH(initial)/2.
 */
double page_entanglement_entropy(double mass_initial_kg, double time_s);

/* ================================================================
 *  ER = EPR CONJECTURE
 * ================================================================ */

/**
 * @brief ER=EPR (Maldacena & Susskind 2013):
 * Einstein-Rosen bridges (wormholes) are equivalent to
 * Einstein-Podolsky-Rosen entanglement.
 *
 * The conjecture: Any two entangled quantum systems are connected
 * by a non-traversable wormhole (ER bridge). The geometry of
 * spacetime emerges from quantum entanglement.
 *
 * Evidence from AdS/CFT: The Ryu-Takayanagi formula relates
 * entanglement entropy in the boundary CFT to minimal surface
 * area in the bulk AdS — geometry = entanglement!
 */

/**
 * @brief Compute the throat radius of an ER bridge between two
 * entangled BHs (or BH + radiation).
 *
 * In the thermofield double state (dual to eternal AdS BH),
 * the two sides are connected by an ER bridge.
 *
 * Throat size ~ l_P × sqrt(S_entanglement / k_B)
 *
 * @param entanglement_entropy S_ent [J/K]
 * @return Throat radius [m]
 */
double er_bridge_throat_radius(double entanglement_entropy);

/**
 * @brief Einstein-Rosen bridge length between two entangled subsystems.
 *
 * For the thermofield double: length grows linearly with time
 * (the "growing wormhole" of Hartman & Maldacena 2013):
 * L(t) ≈ (2π/β) t, where β = 1/T_H.
 *
 * This is the geometric dual of quantum circuit complexity growth.
 */
double er_bridge_length(double temperature, double time_s);

/* ================================================================
 *  FIREWALL PARADOX (AMPS 2012)
 * ================================================================ */

/**
 * @brief The firewall argument (Almheiri, Marolf, Polchinski, Sully 2012):
 *
 * For an old BH (past Page time), the following three statements
 * cannot all be true:
 * (1) Unitarity: information is preserved
 * (2) Effective Field Theory: no drama at the horizon
 * (3) Monogamy of entanglement: a mode cannot be maximally
 *     entangled with two different systems simultaneously
 *
 * AMPS argued that (1)+(3) → (2) must be violated → "firewall"
 * (high-energy quanta at the horizon, destroying infalling observers).
 *
 * This is arguably the most important theoretical tension in
 * theoretical physics today.
 */

/**
 * @brief Check if an "old" BH (past Page time) requires a firewall
 * given assumptions about unitarity and EFT validity.
 *
 * @param mass_initial Initial BH mass [kg]
 * @param time_s Current time [s]
 * @param assume_unitarity 1 if assuming unitarity
 * @param assume_eft 1 if assuming EFT valid at horizon
 * @return 0 (no firewall needed), 1 (firewall required), -1 (paradox unresolved)
 */
int firewall_required(double mass_initial_kg, double time_s,
                       int assume_unitarity, int assume_eft);

/* ================================================================
 *  ADS/CFT AND BLACK HOLES
 * ================================================================ */

/**
 * @brief AdS-Schwarzschild black hole temperature.
 *
 * In AdS_d+1 spacetime, large BHs have T ∝ r_+ (positive heat capacity!)
 * and are THERMODYNAMICALLY STABLE — unlike asymptotically flat BHs.
 *
 * This makes AdS BHs the natural arena for studying BH thermodynamics
 * and the holographic dual of thermal states in the boundary CFT.
 *
 * T = (d r_+² + (d-2)L²) / (4π L² r_+)
 * where L is the AdS curvature radius.
 *
 * For large BHs (r_+ >> L): T ≈ d r_+/(4π L²) (positive C)
 * For small BHs (r_+ << L): T ≈ (d-2)/(4π r_+) (negative C, like flat space)
 *
 * The Hawking-Page transition (1983): at T_HP, the dominant phase
 * switches from thermal AdS (no BH) to large AdS-BH.
 * This is the gravitational dual of the confinement/deconfinement
 * transition in the boundary gauge theory (Witten 1998).
 */
double ads_schwarzschild_temperature(double r_plus, double L, int d);

/**
 * @brief Hawking-Page transition temperature for AdS_d+1.
 * T_HP = (d-1)/(2π L) for d=4 (AdS₅/CFT₄, the N=4 SYM case).
 */
double hawking_page_temperature(double L, int d);

/**
 * @brief Bekenstein-Hawking entropy from AdS/CFT perspective.
 *
 * For N=4 SU(N) SYM at strong coupling:
 * S = (π²/2) N² V T³ (matches AdS₅ BH entropy via holography)
 *
 * The factor N² (number of degrees of freedom) emerges geometrically
 * from the AdS radius: L⁴/l_P⁴ ∝ N².
 */
double ads_cft_entropy(double N, double volume, double temperature);

/* ================================================================
 *  QUANTUM GRAVITY PHENOMENOLOGY
 * ================================================================ */

/**
 * @brief Planck mass: m_P = sqrt(ħc/G) ≈ 2.176 × 10⁻⁸ kg ≈ 1.22 × 10¹⁹ GeV/c².
 *
 * This is the mass scale where quantum gravity effects become O(1).
 * A BH of mass m_P has r_s ≈ l_P: quantum and gravitational regimes meet.
 */
double planck_mass(void);

/**
 * @brief Planck length: l_P = sqrt(ħG/c³) ≈ 1.616 × 10⁻³⁵ m.
 */
double planck_length(void);

/**
 * @brief Planck time: t_P = sqrt(ħG/c⁵) ≈ 5.391 × 10⁻⁴⁴ s.
 */
double planck_time(void);

/**
 * @brief Planck temperature: T_P = sqrt(ħc⁵/Gk_B²) ≈ 1.417 × 10³² K.
 */
double planck_temperature(void);

/**
 * @brief Remnant mass: the mass at which semiclassical evaporation
 * breaks down and quantum gravity effects dominate.
 *
 * M_remnant ≈ m_P (Planck mass) is a natural guess.
 * Some models predict M_remnant >> m_P (large extra dimensions).
 *
 * If stable remnants exist, they could be dark matter candidates
 * (primordial BH remnants).
 */
double remnant_mass(void);

/**
 * @brief Minimum BH mass for semiclassical approximation validity.
 * Below this mass, Hawking's calculation breaks down because
 * T_BH becomes comparable to the mass of emitted particles.
 *
 * M_min ≈ 10¹² kg (primordial BH that would evaporate today).
 */
double semiclassical_mass_limit(void);

#endif /* BLACK_HOLE_ADVANCED_H */
