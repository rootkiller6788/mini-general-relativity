/**
 * @file black_hole_thermodynamics.c
 * @brief Black hole thermodynamics implementations.
 *
 * Covers Hawking temperature, Bekenstein-Hawking entropy,
 * evaporation, four laws, heat capacity, and Page curve.
 *
 * L4 Fundamental Laws: BH thermodynamics laws
 * L6 Canonical Systems: Schwarzschild/Kerr thermodynamic states
 * L8 Advanced Topics: Page curve, information paradox
 */

#include "black_hole_thermodynamics.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 *  HAWKING TEMPERATURE
 * ================================================================ */

double hawking_temperature_schwarzschild(double mass_kg) {
    /* T_H = ħ c³ / (8π G M k_B)
     *
     * Derivation (Hawking 1974):
     * Bogoliubov transformation between in/out vacuum states →
     * thermal spectrum with T = κ/(2π) where κ = 1/(4M) in geometric units.
     * Converting to SI: T = ħc³/(8πGMk_B).
     *
     * For M = M_sun: T ≈ 6.17e-8 K
     * For M = 10^12 kg: T ≈ 1.23e11 K */
    if (mass_kg <= 0) return 0.0;
    return BH_HBAR * BH_C * BH_C * BH_C
           / (8.0 * M_PI * BH_G * mass_kg * BH_KB);
}

double hawking_temperature_kn(double mass_kg, double a_m, double Q_m) {
    /* General KN Hawking temperature:
     * T = ħ κ / (2π k_B c)
     * κ = (r_+ - M) / (r_+² + a²)
     *   = sqrt(M² - a² - Q²) / (2M² - Q² + 2M sqrt(M²-a²-Q²))
     *
     * For extremal (M² = a²+Q²): κ = 0 → T = 0 */
    if (mass_kg <= 0) return 0.0;

    double M = mass_to_geometric(mass_kg);
    double disc = M * M - a_m * a_m - Q_m * Q_m;
    if (disc <= 0) return 0.0; /* Extremal or naked singularity */

    double kappa = sqrt(disc) / (2.0 * M * M - Q_m * Q_m
                    + 2.0 * M * sqrt(disc));
    return BH_HBAR * kappa / (2.0 * M_PI * BH_KB * BH_C);
}

/* ================================================================
 *  BEKENSTEIN-HAWKING ENTROPY
 * ================================================================ */

double bekenstein_hawking_entropy(double horizon_area) {
    /* S_BH = k_B c³ A / (4 G ħ)
     *      = k_B A / (4 l_P²)
     *
     * This is ENORMOUS for astrophysical BHs:
     * M=10 M_sun → A ≈ 1.1×10⁶ m² → S ≈ 2.9×10⁵⁴ J/K
     * Compare: observable universe entropy ≈ 10¹⁰⁴ k_B → BH dominates. */
    if (horizon_area <= 0) return 0.0;
    return BH_KB * BH_C * BH_C * BH_C * horizon_area
           / (4.0 * BH_G * BH_HBAR);
}

double bh_entropy_from_mass(double mass_kg) {
    /* S = 4π k_B G M² / (ħ c) */
    if (mass_kg <= 0) return 0.0;
    return 4.0 * M_PI * BH_KB * BH_G * mass_kg * mass_kg
           / (BH_HBAR * BH_C);
}

double schwarzschild_horizon_area(double mass_kg) {
    /* A = 4π r_s² = 16π G² M² / c⁴ */
    if (mass_kg <= 0) return 0.0;
    double rs = schwarzschild_radius(mass_kg);
    return 4.0 * M_PI * rs * rs;
}

double kerr_horizon_area(double M, double a) {
    /* A = 4π (r_+² + a²)
     *   = 8π M (M + sqrt(M² - a²))
     *
     * For extremal Kerr (a=M): A = 8π M² (half the Schwarzschild area!) */
    if (M <= 0) return 0.0;
    double disc = M * M - a * a;
    if (disc < 0) return 0.0;
    double rp = M + sqrt(disc);
    return 4.0 * M_PI * (rp * rp + a * a);
}

double kn_horizon_area(double M, double a, double Q) {
    /* A = 4π (r_+² + a²) with r_+ = M + sqrt(M² - a² - Q²) */
    if (M <= 0) return 0.0;
    double disc = M * M - a * a - Q * Q;
    if (disc < 0) return 0.0;
    double rp = M + sqrt(disc);
    return 4.0 * M_PI * (rp * rp + a * a);
}

/* ================================================================
 *  BLACK HOLE EVAPORATION
 * ================================================================ */

double hawking_luminosity(double mass_kg) {
    /* P_H = ħ c⁶ / (15360 π G² M²)
     *
     * This is the total power radiated by a Schwarzschild BH,
     * integrating over the thermal spectrum with greybody factors.
     *
     * The greybody factor averaged over all species is ~ 1/(15360π).
     * Including all Standard Model degrees of freedom would increase
     * this by factor ~100 at high T (when T >> all particle masses). */
    if (mass_kg <= 0) return 0.0;
    return BH_HBAR * BH_C * BH_C * BH_C * BH_C * BH_C * BH_C
           / (15360.0 * M_PI * BH_G * BH_G * mass_kg * mass_kg);
}

double mass_loss_rate(double mass_kg) {
    /* dM/dt = -P_H / c² */
    if (mass_kg <= 0) return 0.0;
    return -hawking_luminosity(mass_kg) / (BH_C * BH_C);
}

double evaporation_time(double mass_kg) {
    /* τ = ∫₀^{M} c² dM̃ / P_H(M̃) = 5120 π G² M³ / (ħ c⁴)
     *
     * τ ≈ 2.1e67 × (M/M_sun)³ years
     *
     * For M = 1 kg: τ ≈ 8.4e-17 s (instantaneous explosion!)
     * For M = 10^12 kg: τ ≈ 2.5e10 yr (evaporating now)
     * For M = M_sun: τ ≈ 2.1e67 yr (stable against evaporation) */
    if (mass_kg <= 0) return 0.0;
    return 5120.0 * M_PI * BH_G * BH_G * mass_kg * mass_kg * mass_kg
           / (BH_HBAR * BH_C * BH_C * BH_C * BH_C);
}

double evaporating_temperature(double mass_initial_kg, double time_s) {
    /* T(t) = T₀ / (1 - t/τ)^{1/3}
     * where T₀ = T_H(M₀), τ = evaporation_time(M₀).
     *
     * As t → τ: T → ∞ (classical prediction — quantum gravity
     * effects become important in the final Planck-scale moments). */
    if (mass_initial_kg <= 0 || time_s < 0) return 0.0;
    double tau = evaporation_time(mass_initial_kg);
    if (time_s >= tau) return INFINITY;
    double T0 = hawking_temperature_schwarzschild(mass_initial_kg);
    double frac = 1.0 - time_s / tau;
    return T0 / cbrt(frac * frac * frac); /* = T0 / frac^{1/3} */
}

double evaporating_mass(double mass_initial_kg, double time_s) {
    /* M(t) = M₀ (1 - t/τ)^{1/3} */
    if (mass_initial_kg <= 0 || time_s < 0) return 0.0;
    double tau = evaporation_time(mass_initial_kg);
    if (time_s >= tau) return 0.0;
    return mass_initial_kg * cbrt(1.0 - time_s / tau);
}

/* ================================================================
 *  FOUR LAWS OF BLACK HOLE MECHANICS
 * ================================================================ */

int zeroth_law_check(double kappa1, double kappa2) {
    /* The 0th law: Surface gravity κ is constant over the horizon
     * of a stationary black hole. This is an identity for stationary BHs
     * (Bardeen, Carter, Hawking 1973).
     *
     * For non-stationary (dynamical) horizons, κ may vary spatially,
     * but approaches constant as the BH settles to stationarity. */
    double tol = 1e-10;
    return (fabs(kappa1 - kappa2) < tol * (1.0 + fabs(kappa1))) ? 1 : 0;
}

double first_law_schwarzschild_check(double M_geom, double dM) {
    /* First law (Schwarzschild): dM = (κ/8π) dA = T_H dS
     *
     * Verify: dM - T_H dS ≈ 0
     *
     * For Schwarzschild:
     * κ = 1/(4M), A = 16πM², S = A/(4) = 4πM²
     * → T_H dS = (κ/8π) · 8πM dM = (1/(4M))/8π · 8πM · dM = dM ✓ */
    if (M_geom <= 0) return 0.0;

    /* T_H = κ/(2π) = 1/(8πM) in geometric units */
    double TH = 1.0 / (8.0 * M_PI * M_geom);
    /* dS = 8πM dM (in geometric units, G=ħ=c=k_B=1) */
    double dS = 8.0 * M_PI * M_geom * dM;
    return fabs(dM - TH * dS);
}

double first_law_kn_check(double M, double a, double Q,
                           double dM, double dJ, double dQ) {
    /* First law (KN): dM = T_H dS + Ω_H dJ + Φ_H dQ
     *
     * We compute each term numerically and return the residual:
     * |dM - (T_H dS + Ω_H dJ + Φ_H dQ)| / |dM| */

    double disc = M * M - a * a - Q * Q;
    if (disc <= 0) return -1.0;

    double rp = M + sqrt(disc);

    /* T_H = sqrt(M² - a² - Q²) / [4π M (M + sqrt(M²-a²-Q²))] */
    double TH = sqrt(disc) / (4.0 * M_PI * M * rp);

    /* dS/dM: S = π (r_+² + a²) → dS ≈ 2π r_+ dM (approximate) */
    double dS = 2.0 * M_PI * rp * dM;

    /* Ω_H = a / (r_+² + a²) */
    double OH = a / (rp * rp + a * a);

    /* Φ_H = Q r_+ / (r_+² + a²) */
    double PH = Q * rp / (rp * rp + a * a);

    double residual = dM - (TH * dS + OH * dJ + PH * dQ);
    return (dM != 0) ? fabs(residual / dM) : fabs(residual);
}

double area_increase_theorem(double A_initial, double A_final) {
    /* δA ≥ 0 (Hawking 1971)
     *
     * For BH mergers that produce GW150914-like events:
     * A_initial = A₁ + A₂ (two progenitor BHs)
     * A_final = area of remnant BH
     *
     * Classical GR requires A_final - A_initial ≥ 0.
     *
     * This is the ONLY classical "arrow of time" that is not
     * statistical — it's a geometrical theorem! */
    return A_final - A_initial;
}

double surface_gravity_kn(double M, double a, double Q) {
    /* κ = sqrt(M² - a² - Q²) / (2M² - Q² + 2M sqrt(M²-a²-Q²))
     *
     * As a²+Q² → M² (extremal limit): κ → 0
     * Third law: cannot reach κ=0 in finite physical steps
     * (analog of Nernst's unattainability principle). */
    double disc = M * M - a * a - Q * Q;
    if (disc <= 0) return 0.0;
    double sqrt_disc = sqrt(disc);
    return sqrt_disc / (2.0 * M * M - Q * Q + 2.0 * M * sqrt_disc);
}

/* ================================================================
 *  HEAT CAPACITY AND STABILITY
 * ================================================================ */

double schwarzschild_heat_capacity(double mass_kg) {
    /* C = dE/dT where E = M c²
     *
     * T_H ∝ 1/M → dT_H/dM ∝ -1/M² → C = dM/dT_H ∝ -M²
     *
     * C = -8π² k_B G M² / (ħ c) = -ħ c⁵ / (8π G k_B T²) < 0
     *
     * NEGATIVE heat capacity is the hallmark of self-gravitating
     * systems. It means: lose energy → get HOTTER → radiate more →
     * lose energy faster → thermal runaway. */
    if (mass_kg <= 0) return 0.0;
    return -8.0 * M_PI * M_PI * BH_KB * BH_G * mass_kg * mass_kg
           / (BH_HBAR * BH_C);
}

double kn_heat_capacity(double M, double a, double Q) {
    /* For KN BHs, heat capacity can be positive or negative depending
     * on (M,a,Q). Near extremality, it often becomes positive.
     *
     * Simplified formula for Schwarzschild-like behavior: */
    double TH = hawking_temperature_kn(geometric_to_mass(M), a, Q);
    if (TH <= 0) return 0.0; /* Extremal: C can be finite or diverge */

    /* C = dE/dT, approximate via numerical derivative */
    double dM = M * 1e-6;
    double TH2 = hawking_temperature_kn(geometric_to_mass(M + dM), a, Q);
    double dT = TH2 - TH;
    if (fabs(dT) < 1e-30) return 0.0;
    return dM * BH_C * BH_C / dT; /* Convert geometric mass to energy */
}

double critical_mass_cmb(double T_cmb) {
    /* M_crit where T_H = T_CMB:
     * M_crit = ħ c³ / (8π G k_B T_cmb)
     *
     * For T_cmb = 2.725 K: M_crit ≈ 4.5e22 kg ≈ 0.75% lunar mass.
     *
     * Only primordial BHs lighter than this could have evaporated
     * by now. Heavier BHs gain mass from CMB accretion.
     *
     * This constrains PBH dark matter: M_pbh > ~10^15 g to survive
     * to the present epoch. */
    if (T_cmb <= 0) return INFINITY;
    return BH_HBAR * BH_C * BH_C * BH_C
           / (8.0 * M_PI * BH_G * BH_KB * T_cmb);
}

double irreducible_mass(double M, double a, double Q) {
    /* M_irr² = (1/2) [M² - Q²/2 + M sqrt(M² - a² - Q²)]
     *
     * Christodoulou (1970): M_irr cannot decrease via any classical process.
     * For Schwarzschild: M_irr = M (no extractable energy).
     * For extremal Kerr: M_irr = M/√2 → 29.3% of mass is extractable. */
    double disc = M * M - a * a - Q * Q;
    if (disc < 0) return 0.0;
    double inner = M * M - Q * Q / 2.0 + M * sqrt(disc);
    if (inner < 0) return 0.0;
    return sqrt(inner / 2.0);
}

double extractable_energy_kerr(double M, double a) {
    /* E_extractable = M - M_irr
     * For extremal Kerr (a=M): E/M ≈ 0.293 */
    double Mirr = irreducible_mass(M, a, 0.0);
    if (Mirr <= 0) return 0.0;
    return (M - Mirr) * BH_C * BH_C; /* Convert to joules */
}

/* ================================================================
 *  THERMODYNAMIC STATE COMPUTATIONS
 * ================================================================ */

void schwarzschild_thermo_compute(const SchwarzschildBH *bh,
                                   SchwarzschildThermo *thermo) {
    assert(bh != NULL && thermo != NULL);

    thermo->mass_kg = bh->M_kg;
    thermo->temperature = hawking_temperature_schwarzschild(bh->M_kg);
    thermo->area = schwarzschild_horizon_area(bh->M_kg);
    thermo->entropy = bekenstein_hawking_entropy(thermo->area);
    thermo->surface_gravity = BH_C * BH_C * BH_C * BH_C
                              / (4.0 * BH_G * bh->M_kg);

    /* For Schwarzschild: κ = c⁴/(4GM) */
    thermo->luminosity = hawking_luminosity(bh->M_kg);
    thermo->evaporation_time = evaporation_time(bh->M_kg);
    thermo->heat_capacity = schwarzschild_heat_capacity(bh->M_kg);
}

void kn_thermo_compute(const KerrNewmanBH *bh,
                        KerrNewmanThermo *thermo) {
    assert(bh != NULL && thermo != NULL);

    double M = bh->M;
    double a = bh->a;
    double Q = bh->Q;

    thermo->mass_kg = bh->M_kg;
    thermo->J = bh->J;
    thermo->Q = bh->Q_coulomb;
    thermo->area = kn_horizon_area(M, a, Q);
    thermo->entropy = bekenstein_hawking_entropy(
        thermo->area * BH_C * BH_C * BH_C * BH_C / (BH_G * BH_G));
    thermo->surface_gravity = surface_gravity_kn(M, a, Q) * BH_C * BH_C;
    thermo->omega_H = bh->a / (bh->r_plus * bh->r_plus + bh->a * bh->a)
                      * BH_C; /* Convert geometric angular velocity to rad/s */
    thermo->phi_H = rn_horizon_electric_potential(
        (const ReissnerNordstromBH *)bh); /* Same form */
    thermo->irreducible_mass = irreducible_mass(M, a, Q) * BH_C * BH_C / BH_G;
    thermo->evaporation_time = evaporation_time(bh->M_kg);
    thermo->heat_capacity = kn_heat_capacity(M, a, Q);
    thermo->temperature = hawking_temperature_kn(bh->M_kg, a, Q);
}

void thermo_print_schwarzschild(const SchwarzschildThermo *st) {
    if (st == NULL) return;
    printf("\n========================================\n");
    printf("  SCHWARZSCHILD BLACK HOLE THERMODYNAMICS\n");
    printf("========================================\n");
    printf("Mass:              %12.4e kg\n", st->mass_kg);
    printf("Horizon area:      %12.4e m²\n", st->area);
    printf("Entropy S_BH:      %12.4e J/K\n", st->entropy);
    printf("Temperature T_H:   %12.4e K\n", st->temperature);
    printf("Surface gravity κ: %12.4e m/s²\n", st->surface_gravity);
    printf("Luminosity P_H:    %12.4e W\n", st->luminosity);
    printf("Lifetime τ_evap:   %12.4e s  (%12.4e yr)\n",
           st->evaporation_time, st->evaporation_time / (365.25*86400));
    printf("Heat capacity C:   %12.4e J/K  %s\n",
           st->heat_capacity, st->heat_capacity < 0 ? "(NEGATIVE!)" : "");
}

void thermo_print_kn(const KerrNewmanThermo *kt) {
    if (kt == NULL) return;
    printf("\n========================================\n");
    printf("  KERR-NEWMAN BLACK HOLE THERMODYNAMICS\n");
    printf("========================================\n");
    printf("Mass:              %12.4e kg\n", kt->mass_kg);
    printf("Angular momentum:  %12.4e kg·m²/s\n", kt->J);
    printf("Charge:            %12.4e C\n", kt->Q);
    printf("Horizon area:      %12.4e m²\n", kt->area);
    printf("Entropy S_BH:      %12.4e J/K\n", kt->entropy);
    printf("Temperature T_H:   %12.4e K\n", kt->temperature);
    printf("Ω_H:               %12.4e rad/s\n", kt->omega_H);
    printf("Φ_H:               %12.4e V\n", kt->phi_H);
    printf("Irreducible mass:  %12.4e kg\n", kt->irreducible_mass);
    printf("Lifetime τ_evap:   %12.4e s\n", kt->evaporation_time);
}

/* ================================================================
 *  HAWKING SPECTRUM AND PAGE CURVE
 * ================================================================ */

double hawking_spectrum_scalar(double omega, double T_H,
                                double greybody_factor) {
    /* d²N/(dt dω) = (1/(2π)) · Γ(ω) / (exp(ħω/k_B T_H) - 1)
     *
     * This is the Planck spectrum multiplied by the greybody factor Γ.
     *
     * At low frequencies (ħω << k_B T_H):
     * d²N/(dt dω) ≈ (1/(2π)) Γ(ω) k_B T_H / (ħω)
     *
     * At high frequencies: exponential cutoff exp(-ħω/k_B T_H). */
    if (omega <= 0 || T_H <= 0) return 0.0;
    double x = BH_HBAR * omega / (BH_KB * T_H);
    if (x > 700) return 0.0; /* Avoid overflow */
    double planck = 1.0 / (exp(x) - 1.0);
    return 0.5 / M_PI * greybody_factor * planck;
}

double page_time(double mass_initial_kg) {
    /* t_Page ≈ M₀³ × 5120π G²/(ħ c⁴)
     *
     * This is the time when the BH has radiated HALF its initial entropy.
     * Page (1993): before t_Page, Hawking radiation carries almost no
     * information; after t_Page, the radiation purifies the earlier emission.
     *
     * The "Page curve" solves the apparent paradox: entropy of radiation
     * first rises (thermal phase), then falls (information-emergence phase). */
    if (mass_initial_kg <= 0) return 0.0;
    /* At t_Page, S_BH = S₀/2 → M(t) = M₀/√2 → (1-t/τ)^{1/3} = 1/√2
     * → (1-t/τ) = (1/2)^{3/2} → t/τ = 1 - (1/2)^{3/2} ≈ 0.646 */
    return evaporation_time(mass_initial_kg)
           * (1.0 - pow(0.5, 1.5));
}

double page_curve_entropy(double mass_initial_kg, double time_s) {
    /* S_rad(t):
     * Before t_Page: ≈ S_BH(initial) - S_BH(t) (thermal growth)
     * After t_Page:  ≈ S_BH(t) (information emergence)
     *
     * The Page curve has the famous "inverted V" shape:
     * rises linearly from 0 to S_BH(initial)/2 at t_Page,
     * then falls back to 0 as M(t) → 0. */

    double S_initial = bh_entropy_from_mass(mass_initial_kg);
    double M_t = evaporating_mass(mass_initial_kg, time_s);
    if (M_t <= 0) return S_initial; /* All radiated */

    double S_bh_t = bh_entropy_from_mass(M_t);
    double t_page = page_time(mass_initial_kg);

    if (time_s < t_page) {
        /* Before Page time: radiation entropy = lost BH entropy */
        return S_initial - S_bh_t;
    } else {
        /* After Page time: radiation entropy = remaining BH entropy */
        return S_bh_t;
    }
}
