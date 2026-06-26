/**
 * @file black_hole_advanced.c
 * @brief Advanced black hole topics: no-hair, area theorem, holography,
 *        information paradox, ER=EPR, firewall, AdS/CFT.
 *
 * This file covers L8 Advanced Topics and L9 Research Frontiers.
 * Each function implements an independent theoretical concept at the
 * frontier of modern black hole physics.
 *
 * L4 Fundamental Laws: Area theorem validation
 * L8 Advanced Topics: Holographic principle, Information paradox, ER=EPR
 * L9 Research Frontiers: Firewall, Quantum gravity phenomenology
 */

#include "black_hole_advanced.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 *  NO-HAIR THEOREM
 * ================================================================ */

void kerr_multipoles(double M, double a, KerrMultipoles *mp) {
    /* Kerr multipole moments (Hansen 1974):
     * M_l + i S_l = M (i a)^l
     *
     * Mass moments (even l):
     * M₀ = M,  M₂ = -M a²,  M₄ = M a⁴, ...
     *
     * Current moments (odd l):
     * S₁ = M a,  S₃ = -M a³, ...
     *
     * All determined by just two numbers: M and a → No-hair theorem. */
    assert(mp != NULL);
    mp->M = M;
    mp->J = M * a;   /* S₁ = M a */
    mp->M2 = -M * a * a;
    mp->S3 = -M * a * a * a;
    mp->M4 = M * a * a * a * a;
}

double no_hair_test_quadrupole(double M2, double M, double J) {
    /* For a Kerr BH: M₂ = -J²/M, so M₂·M/J² = -1.
     *
     * This is a TESTABLE prediction of the no-hair theorem.
     * GW observations can measure M₂ (mass quadrupole) independently
     * of M and J from the inspiral-merger-ringdown waveform.
     *
     * Deviations from -1 would indicate:
     * - Non-Kerr compact object (e.g., boson star, gravastar)
     * - Modified gravity (e.g., scalar-tensor, dynamical Chern-Simons)
     * - Quantum gravity corrections to BH structure
     *
     * LIGO-Virgo test (Isi et al. 2019, PRL 123, 111102):
     * For GW150914: M₂/(J²/M) = -1.0 ± 0.3 (consistent with Kerr). */

    if (M <= 0 || J == 0) return 0.0;
    return M2 * M / (J * J);
}

int cosmic_censorship_check(double M, double a, double Q) {
    /* Cosmic censorship (Penrose 1969): Singularities formed in
     * gravitational collapse are always hidden behind an event horizon.
     *
     * For KN BH: a² + Q² ≤ M² → horizon exists.
     * For a² + Q² > M² → naked singularity (no horizon).
     *
     * Status: Cosmic censorship remains a CONJECTURE, not proven
     * in full generality. Counterexamples exist in specialized
     * circumstances (e.g., Choptuik critical collapse, certain
     * super-extremal charged configurations).
     *
     * Weak cosmic censorship (Penrose): "Nature abhors a naked singularity."
     * Strong cosmic censorship (Penrose 1979): Singularities are always
     * spacelike (inextendable), so no observer can "see" them. */

    return (M * M >= a * a + Q * Q) ? 1 : 0;
}

/* ================================================================
 *  AREA THEOREM
 * ================================================================ */

int area_theorem_merger_check(double M1, double a1, double Q1,
                               double M2, double a2, double Q2,
                               double M_final, double a_final, double Q_final) {
    /* Check A_final ≥ A₁ + A₂ for binary BH merger.
     *
     * For GW150914:
     * M₁ ≈ 36 M_sun, a₁ ≈ 0.3 (estimated)
     * M₂ ≈ 29 M_sun, a₂ ≈ 0.3
     * M_final ≈ 62 M_sun, a_final ≈ 0.68
     *
     * A₁ + A₂ ≈ 2.4×10⁵ km²
     * A_final ≈ 3.6×10⁵ km² → ΔA > 0 ✓
     *
     * LIGO-Virgo directly tested the area theorem with GW150914 data
     * (Isi et al., PRL 127, 011103, 2021): confirmed with ~97% probability. */

    double A1 = kn_horizon_area(M1, a1, Q1);
    double A2 = kn_horizon_area(M2, a2, Q2);
    double Af = kn_horizon_area(M_final, a_final, Q_final);

    double delta_A = Af - (A1 + A2);
    return (delta_A >= -1e-10) ? 1 : 0;
}

double hawking_area_bound_energy(double m1, double m2) {
    /* From area theorem: M_final² ≥ m₁² + m₂² (for Schwarzschild)
     *
     * E_rad_max = (m₁ + m₂) - sqrt(m₁² + m₂²)
     *
     * This is the HAWKING BOUND on GW energy emission.
     *
     * For equal masses m₁=m₂=m:
     * E_max = 2m - m√2 = (2 - √2)m ≈ 0.586m
     * As a fraction of total: (2-√2)/2 ≈ 0.293 = 29.3%
     *
     * This is an UPPER BOUND from classical GR. Quantum effects
     * (Hawking radiation) can add extra energy loss but it's
     * negligible for astrophysical BHs. */

    if (m1 <= 0 || m2 <= 0) return 0.0;
    return (m1 + m2) - sqrt(m1 * m1 + m2 * m2);
}

/* ================================================================
 *  HOLOGRAPHIC PRINCIPLE
 * ================================================================ */

double bekenstein_bound(double R, double E) {
    /* Bekenstein (1981): For any weakly gravitating system:
     * S ≤ 2π k_B R E / (ħ c)
     *
     * This is the UNIVERSAL entropy bound. It was derived from
     * the GSL (Generalized Second Law): ΔS_total = ΔS_matter + ΔS_BH ≥ 0.
     *
     * For a human (R≈1m, E≈10¹⁰ J): S_max ≈ 2.6×10⁴² k_B,
     * actual entropy ~10²⁸ k_B → far below the bound.
     *
     * The bound is saturated by BHs: S_BH = 2π k_B (2M) E / (ħ c)
     * with R = 2M, E = M c² → S = 4π k_B M² = A/(4 l_P²) ✓ */
    if (R <= 0 || E <= 0) return 0.0;
    return 2.0 * M_PI * BH_KB * R * E / (BH_HBAR * BH_C);
}

double holographic_bound_bits(double area_surface) {
    /* N_bits_max = A / (4 l_P² ln 2)
     *
     * One bit of information requires AT LEAST 4 ln(2) l_P² of area.
     * Since l_P² ≈ 2.612×10⁻⁷⁰ m², this is TINY.
     *
     * For a 1 cm² surface: N_bits ≈ 1.0×10⁶⁵ — enormous!
     * The holographic bound is NEVER saturated in ordinary matter.
     * Only BHs approach it (by definition: S_BH = A/(4 l_P²)). */

    if (area_surface <= 0) return 0.0;
    double lP = planck_length();
    return area_surface / (4.0 * lP * lP * log(2.0));
}

double planck_area(void) {
    /* l_P² = ħ G / c³ ≈ 2.612×10⁻⁷⁰ m²
     *
     * The Planck area is the quantum of area in many QG approaches.
     * Loop Quantum Gravity predicts discrete area spectrum:
     * A = 8π γ l_P² Σ √[j(j+1)], with γ ≈ 0.2375 (Immirzi parameter). */
    return BH_HBAR * BH_G / (BH_C * BH_C * BH_C);
}

double spherical_holographic_bound(double radius_m) {
    /* For a sphere of radius R: S_max ≤ A/(4 l_P²) = π R² / l_P²
     *
     * The maximum entropy scales as R² (area), not R³ (volume).
     * This is RADICAL — it means our 3D world is holographically
     * encoded on a 2D surface. */
    if (radius_m <= 0) return 0.0;
    return M_PI * radius_m * radius_m / planck_area();
}

double holographic_bound_margin(double entropy, double area_m2) {
    /* Tolerated entropy: S_max = k_B A / (4 l_P²)
     * Margin = S_max - S. Positive → bound holds. */
    double S_max = BH_KB * area_m2 / (4.0 * planck_area());
    return S_max - entropy;
}

/* ================================================================
 *  BLACK HOLE INFORMATION PARADOX
 * ================================================================ */

double page_information(double mass_initial_kg, double time_s) {
    /* I(t) = information emerged in Hawking radiation by time t.
     *
     * Before Page time: I(t) ≈ 0 (radiation is nearly thermal,
     *   carries minimal information)
     *
     * After Page time: I(t) ≈ S_initial - S_BH(t) ≈ S_rad(t)
     *   (information emerges at rate ~1 bit per emitted quantum)
     *
     * In bits: I_bits = I(t) / (k_B ln 2) */

    double S_initial = bh_entropy_from_mass(mass_initial_kg);
    double M_t = evaporating_mass(mass_initial_kg, time_s);
    double S_bh_t = (M_t > 0) ? bh_entropy_from_mass(M_t) : 0.0;
    double t_page = page_time(mass_initial_kg);

    if (time_s >= evaporation_time(mass_initial_kg)) {
        /* BH fully evaporated: all information returned */
        return S_initial / (BH_KB * log(2.0)); /* In bits */
    }

    if (time_s < t_page) {
        /* Before Page time: only ~O(1) bits emerged */
        double frac = time_s / t_page;
        return S_initial * frac * 0.01 / (BH_KB * log(2.0));
    } else {
        /* After Page time: S_initial - S_bh_t information emerged */
        return (S_initial - S_bh_t) / (BH_KB * log(2.0));
    }
}

int has_reached_page_time(double mass_initial_kg, double time_s) {
    return (time_s >= page_time(mass_initial_kg)) ? 1 : 0;
}

double page_entanglement_entropy(double mass_initial_kg, double time_s) {
    /* S_entanglement(t) = min(S_BH(t), S_radiation(t))
     *
     * The PAGE CURVE:
     * - Starts at 0 (pure BH + no radiation)
     * - Rises to S_initial/2 at t = t_Page
     * - Falls to 0 as M → 0 (pure radiation, information preserved)
     *
     * This inverted-V shape is the key prediction that resolves
     * the apparent information paradox within semiclassical gravity
     * (Page 1993, PRL 71, 3743). */

    double S_initial = bh_entropy_from_mass(mass_initial_kg);
    double M_t = evaporating_mass(mass_initial_kg, time_s);
    if (M_t <= 0) return 0.0; /* Fully evaporated: pure state */

    double S_bh_t = bh_entropy_from_mass(M_t);
    double S_rad = S_initial - S_bh_t;

    return (S_bh_t < S_rad) ? S_bh_t : S_rad;
}

/* ================================================================
 *  ER = EPR CONJECTURE
 * ================================================================ */

double er_bridge_throat_radius(double entanglement_entropy) {
    /* ER = EPR (Maldacena & Susskind 2013):
     * Einstein-Rosen bridges (wormholes) are geometrically equivalent
     * to Einstein-Podolsky-Rosen quantum entanglement.
     *
     * A maximally entangled pair of BHs are connected by an ER bridge.
     * The throat radius r_throat ~ l_P × √(S_entanglement / k_B).
     *
     * For two maximally entangled BHs of entropy S each:
     * r_throat ~ √(4G S / π c³) → the same as the BH horizon radius!
     *
     * This is the most concrete realization of "it from qubit"
     * (Wheeler's vision): spacetime geometry emerges from
     * quantum information. */

    if (entanglement_entropy <= 0) return 0.0;
    double lP = planck_length();
    return lP * sqrt(entanglement_entropy / BH_KB);
}

double er_bridge_length(double temperature, double time_s) {
    /* ER bridge length grows linearly with time in the
     * thermofield double state (dual to eternal AdS BH):
     *
     * L(t) ≈ (2π/β) t = 2π k_B T t / ħ
     *
     * This geometric growth is the holographic dual of
     * quantum circuit complexity growth ("Complexity = Volume"
     * or "Complexity = Action" conjectures).
     *
     * The wormhole grows FOREVER even though the dual QFT
     * has thermalized — this is the "dual" of the fact that
     * quantum complexity continues to increase long after
     * thermalization (Susskind 2014). */

    if (temperature <= 0) return 0.0;
    double beta = BH_HBAR / (BH_KB * temperature); /* ħ/(k_B T) */
    if (beta <= 0) return 0.0;
    return (2.0 * M_PI / beta) * time_s * BH_C; /* Convert to length */
}

/* ================================================================
 *  FIREWALL PARADOX (AMPS 2012)
 * ================================================================ */

int firewall_required(double mass_initial_kg, double time_s,
                       int assume_unitarity, int assume_eft) {
    /* The AMPS (Almheiri, Marolf, Polchinski, Sully 2012) argument:
     *
     * For an OLD black hole (past Page time), consider three
     * mutually inconsistent assumptions:
     *
     * (A) Unitarity: Information is preserved (pure→pure evolution)
     * (B) EFT at horizon: No drama, vacuum state for infalling observer
     * (C) Monogamy of entanglement: Late radiation mode B cannot be
     *     simultaneously maximally entangled with both:
     *     - Early radiation mode A (required by unitarity + Page curve)
     *     - Interior mode C (required by EFT at horizon)
     *
     * AMPS: (A)+(C) → ¬(B) → FIREWALL (high-energy quanta at horizon)
     * Others: (B)+(C) → ¬(A) → Information loss (Hawking's original view)
     * Maldacena+Susskind: ER=EPR resolves the paradox — B and C are
     *   the SAME degree of freedom (ER bridge = entanglement).
     *
     * Current status: UNRESOLVED. This is THE central puzzle of
     * theoretical physics in the 21st century. */

    int is_old = has_reached_page_time(mass_initial_kg, time_s);

    if (!is_old) return 0; /* Young BH: no paradox */

    if (assume_unitarity && assume_eft) {
        /* AMPS: Both cannot hold simultaneously for old BH
         * → firewall required (or new physics) */
        return 1;
    }
    if (!assume_unitarity) {
        /* Information is lost (Hawking's concession withdrawn in 2004,
         * but some still explore this) */
        return 0;
    }
    if (!assume_eft) {
        /* EFT breaks down (e.g., fuzzball, firewall, nonlocal effects) */
        return 1;
    }
    return -1; /* Unresolved */
}

/* ================================================================
 *  ADS/CFT AND BLACK HOLES
 * ================================================================ */

double ads_schwarzschild_temperature(double r_plus, double L, int d) {
    /* AdS_d+1 Schwarzschild BH temperature:
     *
     * T = (d r_+² + (d-2)L²) / (4π L² r_+)
     *
     * For large BHs (r_+ >> L): T ≈ d r_+/(4π L²) → C_V > 0 (stable!)
     * For small BHs (r_+ << L): T ≈ (d-2)/(4π r_+) → C_V < 0 (like flat space)
     *
     * This is CRUCIAL: large AdS BHs are thermodynamically STABLE.
     * They can exist in equilibrium with their Hawking radiation
     * (AdS acts as a confining box). This enables the holographic
     * duality with thermal CFT states. */

    if (r_plus <= 0 || L <= 0 || d < 2) return 0.0;
    return (d * r_plus * r_plus + (d - 2.0) * L * L)
           / (4.0 * M_PI * L * L * r_plus);
}

double hawking_page_temperature(double L, int d) {
    /* Hawking-Page transition temperature:
     *
     * T_HP = (d-1)/(2π L) for d=4 (AdS₅/CFT₄)
     *
     * Below T_HP: thermal AdS dominates (no BH, confined phase in CFT)
     * Above T_HP: large AdS BH dominates (deconfined phase)
     *
     * This is the gravitational dual of the confinement/deconfinement
     * phase transition in large-N gauge theories (Witten 1998).
     *
     * For QCD (N=3, not large-N): the deconfinement transition
     * is seen at T ≈ 150-170 MeV in heavy-ion collisions at RHIC/LHC. */
    if (L <= 0 || d < 2) return 0.0;
    return (d - 1.0) / (2.0 * M_PI * L);
}

double ads_cft_entropy(double N, double volume, double temperature) {
    /* N=4 SU(N) SYM at strong coupling ('t Hooft coupling λ >> 1):
     * S = (π²/2) N² V T³
     *
     * This matches the Bekenstein-Hawking entropy of an AdS₅ BH
     * at Hawking temperature T, computed via the Bekenstein-Hawking
     * formula S = A/(4G₅), confirming the AdS/CFT correspondence
     * at the level of thermodynamics.
     *
     * The N² scaling is the hallmark of large-N gauge theories:
     * N² "gluon" degrees of freedom per unit volume. */

    if (N <= 0 || volume <= 0 || temperature <= 0) return 0.0;
    double pi2 = M_PI * M_PI;
    return (pi2 / 2.0) * N * N * volume * temperature * temperature * temperature;
}

/* ================================================================
 *  QUANTUM GRAVITY PHENOMENOLOGY
 * ================================================================ */

double planck_mass(void) {
    /* m_P = √(ħ c / G) ≈ 2.176×10⁻⁸ kg ≈ 1.22×10¹⁹ GeV/c²
     *
     * About 22 micrograms — comparable to a flea egg.
     *
     * The Planck mass is the ONLY mass scale that can be formed
     * from ħ, c, G alone (dimensional analysis).
     *
     * Hierarchy problem: Why is m_Planck >> m_Higgs?
     * m_P / m_weak ≈ 10¹⁶ — this enormous ratio is unexplained
     * in the Standard Model. */
    return sqrt(BH_HBAR * BH_C / BH_G);
}

double planck_length(void) {
    /* l_P = √(ħ G / c³) ≈ 1.616×10⁻³⁵ m
     *
     * Planck length is the scale where quantum fluctuations of
     * spacetime become O(1): δg_{μν} ~ l_P / r ~ 1.
     *
     * If extra dimensions exist with radius R (ADD model, Arkani-Hamed et al. 1998),
     * the true Planck scale can be much lower (potentially TeV scale),
     * explaining the hierarchy problem geometrically. */
    return sqrt(BH_HBAR * BH_G / (BH_C * BH_C * BH_C));
}

double planck_time(void) {
    /* t_P = √(ħ G / c⁵) ≈ 5.391×10⁻⁴⁴ s
     *
     * The Planck time is the time for light to cross the Planck length.
     * It's the natural "clock tick" of quantum gravity.
     *
     * Our universe is ~8×10⁶⁰ t_P old — that's alot of ticks. */
    return sqrt(BH_HBAR * BH_G / (BH_C * BH_C * BH_C * BH_C * BH_C));
}

double planck_temperature(void) {
    /* T_P = √(ħ c⁵ / G k_B²) ≈ 1.417×10³² K
     *
     * This is the temperature at which thermal energy k_B T equals
     * the Planck energy. At T_P, quantum gravity effects are unsuppressed.
     *
     * The universe was at T_P at t ~ t_P (the Planck epoch),
     * when all four fundamental forces may have been unified. */
    double num = BH_HBAR * BH_C * BH_C * BH_C * BH_C * BH_C / BH_G;
    return sqrt(num) / BH_KB;
}

double remnant_mass(void) {
    /* M_remnant ≈ m_P (Planck mass)
     *
     * When an evaporating BH reaches M ~ m_P, the semiclassical
     * approximation (Hawking 1974-1975) breaks down.
     *
     * Possibilities for the end state:
     * 1. Complete evaporation (no remnant): M → 0
     * 2. Stable Planck-mass remnant: evaporation stops
     *    → Dark matter candidate (MacGibbon 1987)
     * 3. Naked singularity: horizon disappears at finite mass
     *
     * If Planck-mass remnants are stable and were produced in the
     * early universe, they could constitute dark matter. However,
     * strong constraints exist from gamma-ray observations. */

    return planck_mass();
}

double semiclassical_mass_limit(void) {
    /* M_min for Hawking's semiclassical calculation to be valid.
     *
     * Condition: T_H << m_π c²/k_B (pion rest mass)
     * → M >> ~10¹² kg
     *
     * A 10¹² kg primordial BH would be evaporating TODAY.
     * Observed gamma-ray background constrains the abundance
     * of such PBHs (Carr et al. 2010).
     *
     * Constraints from extragalactic gamma-ray background (EGRB)
     * and 511 keV Galactic positron line restrict PBH DM fraction
     * to f_PBH < 10⁻⁸ for M_PBH ~ 10¹⁵ g. */

    /* M where T_H = m_e c²/k_B (electron rest mass),
     * above which electron-positron pair emission is suppressed */
    double m_e = 9.10938356e-31; /* electron mass [kg] */
    double T_e = m_e * BH_C * BH_C / BH_KB; /* ~5.93×10⁹ K */
    return BH_HBAR * BH_C * BH_C * BH_C
           / (8.0 * M_PI * BH_G * BH_KB * T_e);
}

/* Computed once at the beginning to provide information.
 * We are using cbrt() for cube root. */
#ifndef cbrt
#define cbrt(x) pow((x), 1.0/3.0)
#endif
