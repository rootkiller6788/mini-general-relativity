/**
 * @file bh_thermo.c
 * @brief Black hole thermodynamics: Bekenstein-Hawking entropy, temperature.
 *
 * Reference: Wald (1984), Ch.7; Carroll (2004), Ch.8 appendix
 *            Bekenstein "Black Holes and Entropy" (1973)
 *            Hawking "Particle Creation by Black Holes" (1975)
 *            Bardeen, Carter, Hawking "The Four Laws of BH Mechanics" (1973)
 *
 * Black hole thermodynamics establishes a deep connection between GR,
 * quantum mechanics, and statistical mechanics.
 */

#include <math.h>
#include "coordinate.h"
#include "einstein.h"

/* Physical constants in SI */
#define HBAR 1.054571817e-34  /* Planck's reduced constant (J·s) */
#define KB   1.380649e-23     /* Boltzmann constant (J/K) */

/* Planck quantities */
#define PLANCK_MASS    2.176434e-8    /* kg */
#define PLANCK_LENGTH  1.616255e-35   /* m */
#define PLANCK_TIME    5.391247e-44   /* s */
#define PLANCK_TEMP    1.416784e32    /* K */

/* ========================================================================
 * L8: Schwarzschild black hole thermodynamics
 * ========================================================================*/

/**
 * Bekenstein-Hawking entropy (Schwarzschild):
 *   S_BH = k_B * A / (4 ℓ_P²)
 *
 * where A = 16π M² is the horizon area (geometric units) and
 * ℓ_P = √(ℏ G / c³) is the Planck length.
 *
 * In SI units:
 *   S_BH = k_B * c³ * A / (4 G ℏ)
 *        = (k_B * π * c³ / (G ℏ)) * r_s²
 *        = k_B * A / (4 ℓ_P²)
 *
 * Knowledge point: Black hole entropy (Bekenstein 1973, Hawking 1974)
 * is proportional to the horizon AREA (not volume), suggesting that
 * the information content of a black hole is encoded on its surface.
 * This is the origin of the holographic principle.
 *
 * For a solar-mass BH (M~2e30 kg): S_BH ~ 10^77 k_B (enormous!)
 * Compare: Sun's thermal entropy ~ 10^57 k_B
 */
double bekenstein_hawking_entropy_schwarzschild(double M_kg)
{
    /* r_s = 2 G M / c² */
    double r_s = 2.0 * G_NEWTON * M_kg / (C_LIGHT * C_LIGHT);

    /* A = 4π r_s² */
    double A = 4.0 * M_PI * r_s * r_s;

    /* S = k_B * c³ * A / (4 G ℏ) */
    double S = KB * C_LIGHT * C_LIGHT * C_LIGHT * A
               / (4.0 * G_NEWTON * HBAR);

    return S; /* J/K */
}

/**
 * Hawking temperature (Schwarzschild):
 *   T_H = ℏ c³ / (8π G M k_B)
 *
 * In terms of surface gravity κ = c⁴/(4GM):
 *   T_H = ℏ κ / (2π k_B c)
 *
 * For a solar-mass BH: T_H ~ 6 × 10^{-8} K (colder than CMB!)
 * For a micro BH with M = 10^{12} kg: T_H ~ 10^{11} K (hot!)
 *
 * Knowledge point: Hawking radiation (1974) arises from quantum field
 * theory in curved spacetime. The temperature is inversely proportional
 * to mass — small black holes are hotter and evaporate faster.
 *
 * This is NOT a purely classical GR effect — it requires quantum fields
 * on a classical curved background. Full theory of quantum gravity is
 * needed to describe the final stages of evaporation.
 */
double hawking_temperature_schwarzschild(double M_kg)
{
    return HBAR * C_LIGHT * C_LIGHT * C_LIGHT
           / (8.0 * M_PI * G_NEWTON * M_kg * KB);
}

/**
 * Surface gravity (Schwarzschild):
 *   κ = c⁴ / (4GM)
 *
 * The surface gravity is the force per unit mass that must be exerted
 * at infinity to hold a test particle stationary at the horizon.
 *
 * Knowledge point: The surface gravity κ is constant over the horizon
 * (zeroth law of BH mechanics). For Schwarzschild: κ = 1/(4M) in
 * geometric units. It plays the role of temperature in the classical
 * laws of BH mechanics.
 */
double surface_gravity_schwarzschild(double M_kg)
{
    double r_s = 2.0 * G_NEWTON * M_kg / (C_LIGHT * C_LIGHT);
    return C_LIGHT * C_LIGHT * C_LIGHT * C_LIGHT
           / (4.0 * G_NEWTON * M_kg);
    /* Equivalent: c⁴ / (2GM) * (1/2) → c⁴/(4GM) */
    /* Actually standard formula: κ = c⁴/(4GM) */
    /* In geometric units: κ = 1/(4M) */
    /* Check: at r_s = 2GM/c², κ = c²/(2r_s) = c⁴/(4GM) ✓ */
    (void)r_s; /* r_s computed but formula uses M directly */
}

/**
 * Black hole evaporation time (Schwarzschild):
 *
 * From Stefan-Boltzmann law for a black body:
 *   dM/dt = -ℏ c⁴ / (15360 π G² M²)
 *
 * Integrating from initial mass M to 0:
 *   t_evap = (5120 π G² M³) / (ℏ c⁴)
 *
 * For M = M_sun: t_evap ~ 10^67 years (far longer than universe age ~14 Gyr)
 * For M = 10^{12} kg (primordial BH): t_evap ~ age of universe (~14 Gyr)
 *
 * Knowledge point: BH evaporation provides a natural cutoff for
 * primordial black hole dark matter candidates.
 * PBHs with M < 10^{12} kg would have evaporated by now.
 * PBHs with M ~ 10^{15} kg could be dark matter (evaporation time > age).
 */
double black_hole_evaporation_time_s(double M_kg)
{
    double G2 = G_NEWTON * G_NEWTON;
    double M3 = M_kg * M_kg * M_kg;
    return 5120.0 * M_PI * G2 * M3
           / (HBAR * pow(C_LIGHT, 4.0));
}

/**
 * Black hole luminosity (Hawking radiation power):
 *   P = ℏ c⁶ / (15360 π G² M²)
 *
 * For Schwarzschild, treating the BH as a black body with temperature T_H
 * and area A: P = σ A T_H⁴ where σ = π² k_B⁴ / (60 ℏ³ c²).
 *
 * Knowledge point: The power radiated scales as 1/M², so as the BH
 * evaporates and M decreases, the power increases dramatically
 * (runaway evaporation), culminating in a final burst.
 */
double hawking_luminosity_schwarzschild(double M_kg)
{
    return HBAR * pow(C_LIGHT, 6.0)
           / (15360.0 * M_PI * G_NEWTON * G_NEWTON * M_kg * M_kg);
}

/**
 * Entropy for Kerr black hole:
 *   S = (k_B c³ / (4Gℏ)) * A
 *   A = 8π M (M + √(M² - a²))  (Bekenstein-Hawking area)
 *
 * Note: For extremal Kerr (a=M): A = 8πM² (half of Schwarzschild!)
 *       S_extremal = (1/2) S_schwarzschild
 *
 * Knowledge point: Kerr BHs have lower entropy than Schwarzschild of
 * same mass because rotation "stores" energy that reduces the needed
 * horizon area. The area is a measure of the BH's hidden information.
 */
double kerr_horizon_area(double M, double a)
{
    double r_plus = kerr_horizon_radius(M, a);
    /* Standard formula: A = 8π M r_+  (equivalent to A = 4π(r_+² + a²)) */
    double A = 8.0 * M_PI * M * r_plus;
    return A;
}

double bekenstein_hawking_entropy_kerr(double M, double a)
{
    double A = kerr_horizon_area(M, a);
    return KB * pow(C_LIGHT, 3.0) * A / (4.0 * G_NEWTON * HBAR);
}

/**
 * Hawking temperature for Kerr black hole:
 *   T_H = (ℏ / (2π k_B)) * κ
 *   κ = √(M² - a²) / (2M (M + √(M² - a²)))
 *
 * In terms of r_±:
 *   κ = (r_+ - r_-) / (2(r_+² + a²))
 *
 * For extremal Kerr (a→M): κ → 0, T_H → 0
 * → Extremal BHs are stable against Hawking radiation.
 *
 * Knowledge point: The third law of BH thermodynamics (analogous to
 * Nernst's theorem) states that κ cannot be reduced to zero in a
 * finite number of physical processes (though extremal BHs have κ=0,
 * they require infinite proper time to form).
 */
double hawking_temperature_kerr(double M, double a)
{
    if (a >= M) return 0.0; /* extremal: T=0 */

    double r_plus = M + sqrt(M*M - a*a);
    double r_minus = M - sqrt(M*M - a*a);

    /* Surface gravity κ = (r_+ - r_-) / (2(r_+² + a²)) */
    double kappa = (r_plus - r_minus) / (2.0 * (r_plus*r_plus + a*a));

    /* T = ℏ κ / (2π k_B c) */
    return HBAR * kappa / (2.0 * M_PI * KB * C_LIGHT);
}

/* ========================================================================
 * L8: Four laws of black hole mechanics
 * ========================================================================*/

/**
 * Zeroth law: Surface gravity κ is constant over the event horizon
 * of a stationary black hole.
 *
 * This is analogous to the zeroth law of thermodynamics (temperature
 * is constant throughout a system in thermal equilibrium).
 *
 * Verification for Kerr: κ depends only on M and a, which are constant
 * over the horizon. ✓
 *
 * @return 1 if law holds, 0 if violated (always 1 for Kerr).
 */
int zeroth_law_bh_mechanics(double M, double a)
{
    /* Surface gravity κ for Kerr is uniform over the horizon.
     * It depends only on M and a, which are global parameters.
     * This is always true for the Kerr solution. */
    (void)M; (void)a;
    return 1;
}

/**
 * First law (differential mass formula):
 *   dM = (κ/8π) dA + Ω_H dJ + Φ_H dQ
 *
 * where:
 *   Ω_H = a/(r_+² + a²) = angular velocity of the horizon
 *   Φ_H = Q r_+/(r_+² + a²) = electrostatic potential of horizon
 *
 * For uncharged Kerr (Q=0):
 *   dM = T_H dS + Ω_H dJ
 *   with T_H dS = (κ/2π) * (k_B c³ A / (4Gℏ)) * ... → κ dA/(8πG)
 *
 * Knowledge point: The first law of BH mechanics is formally identical
 * to the first law of thermodynamics dE = T dS + Ω dJ + Φ dQ when:
 *   T_H = ℏκ/(2πk_B c), S_BH = k_B A/(4ℓ_P²).
 *
 * @return energy change dM from given dA and dJ (geometric units G=c=1).
 */
double first_law_bh_mechanics(double M, double a, double dA, double dJ)
{
    double r_plus = M + sqrt(M*M - a*a);
    double kappa = (r_plus - (M - sqrt(M*M - a*a))) / (2.0*(r_plus*r_plus + a*a));
    double Omega_H = a / (r_plus * r_plus + a * a);

    return (kappa / (8.0 * M_PI)) * dA + Omega_H * dJ;
}

/**
 * Second law (area theorem):
 *   The area A of the event horizon never decreases in any classical process:
 *   dA/dt ≥ 0
 *
 * (Hawking 1971)
 *
 * This is the analog of the second law of thermodynamics.
 * Quantum effects (Hawking radiation) can decrease the area,
 * leading to the generalized second law:
 *   d(S_BH + S_matter)/dt ≥ 0
 *
 * Knowledge point: The area theorem is a purely classical GR result.
 * It follows from the Raychaudhuri equation + null energy condition.
 * The fact that Hawking radiation violates it was the original clue
 * that BHs have temperature: if dA can decrease, and dA ~ dS, then
 * BHs must radiate to satisfy the second law.
 *
 * @return 1 if area never decreases in the process, 0 if it does.
 */
int second_law_bh_check(double A_initial, double A_final)
{
    return (A_final >= A_initial - 1e-15);
}

/**
 * Third law:
 *   It is impossible by any physical process to reduce κ to zero
 *   in a finite sequence of operations.
 *
 * (Israel 1986)
 *
 * Analogy: Nernst's theorem — absolute zero is unreachable in a
 * finite number of thermodynamic processes.
 *
 * Extremal Kerr (κ=0) requires a=M, which would require infinite
 * proper time to achieve via finite physical processes.
 *
 * @return 1 if process is consistent with third law, 0 if it would
 *         violate (i.e., drive κ to exactly zero).
 */
int third_law_bh_check(double kappa_initial, double kappa_final)
{
    /* κ > 0 for non-extremal, cannot reach exactly 0 in finite steps */
    if (kappa_final <= 0.0 && kappa_initial > 0.0) return 0;
    return 1;
}

/**
 * Information paradox summary function.
 *
 * The BH information paradox arises from the apparent conflict between:
 *   1. Unitarity of quantum mechanics (information is preserved)
 *   2. Hawking radiation is purely thermal (information appears lost)
 *
 * Proposed resolutions:
 *   - Information is encoded in subtle correlations of Hawking radiation
 *   - Information is stored in a Planck-scale remnant
 *   - Information escapes through the firewall at the horizon
 *   - AdS/CFT: information is preserved in the dual CFT
 *
 * This is a L9 research frontier topic.
 *
 * @return 0 (unresolved — open problem in quantum gravity)
 */
int information_paradox_status(void)
{
    /* The paradox remains an active area of research.
     * Recent consensus (Maldacena, Susskind, 't Hooft, Page, etc.)
     * leans toward information preservation via complementarity
     * or holography, but no complete resolution exists. */
    return 0; /* UNRESOLVED */
}
