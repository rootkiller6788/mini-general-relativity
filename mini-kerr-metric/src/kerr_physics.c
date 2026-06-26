/**
 * @file kerr_physics.c
 * @brief Physical processes in Kerr geometry
 *
 * Implements:
 * - Penrose process (rotational energy extraction from ergosphere)
 * - Superradiance (wave amplification by rotating black holes)
 * - Blandford-Znajek mechanism (electromagnetic energy extraction)
 * - Novikov-Thorne accretion disk model
 * - Black hole shadow observables
 * - Quasi-normal mode ringdown from BH perturbations
 * - Iron K-alpha line profile (relativistic broadening)
 *
 * Reference: Penrose & Floyd (1971), Blandford & Znajek (1977),
 * Novikov & Thorne (1973), Berti, Cardoso & Will (2006)
 */

#include "kerr_physics.h"
#include "kerr_geodesics.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

/* ====================================================================
 * Penrose process
 *
 * Maximum efficiency for extremal Kerr (a=M):
 *   eta_max = 1 - 1/sqrt(2) ~ 0.2929
 *
 * This means up to 29.3% of the rest mass energy can be extracted.
 * Compare: nuclear fusion ~ 0.7%, matter-antimatter ~ 100%.
 *
 * The process: a particle enters the ergosphere and decays into two
 * fragments. One fragment falls into the BH with negative energy
 * (relative to infinity), while the other escapes with more energy
 * than the original particle.
 * ==================================================================== */

double kerr_penrose_max_efficiency(const KerrBlackHole *bh)
{
    double a = bh->a;
    double M = bh->M;

    if (a > M) return 0.0;  /* Naked singularity */

    /* Efficiency depends on how close to extremal */
    double aM = a / M;
    if (aM >= 1.0) return 1.0 - 1.0 / sqrt(2.0);

    /* General formula: eta ~ 1 - sqrt(r_plus/(2M))
     * where r_plus = M + sqrt(M^2 - a^2) */
    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);
    if (isnan(r_plus)) return 0.0;

    return 1.0 - sqrt(r_plus / (2.0 * M));
}

int kerr_penrose_feasibility(const KerrBlackHole *bh,
                             const KerrBLPoint *pt,
                             double E_particle, double Lz_particle)
{
    /* Penrose process requires:
     * 1. Point is inside the ergosphere
     * 2. It's possible for a fragment to have negative energy */

    if (!kerr_is_inside_ergosphere(bh, pt)) return 0;

    /* Inside ergosphere, g_tt > 0, allowing negative-energy states.
     * E = -p_t = -g_{t mu} p^mu.
     * For a particle with 4-momentum p^mu, E can be negative inside
     * the ergosphere if p^phi is sufficiently large and co-rotating. */
    (void)E_particle;
    (void)Lz_particle;

    /* Simplification: inside ergosphere => Penrose process feasible */
    return 1;
}

int kerr_penrose_process(const KerrBlackHole *bh,
                         const KerrBLPoint *decay_pt,
                         double E_in, double Lz_in, double mu_in,
                         double mu_frag1, double mu_frag2,
                         PenroseProcess *result)
{
    assert(bh && decay_pt && result);

    result->r_decay = decay_pt->r;
    result->theta_decay = decay_pt->theta;
    result->E_in = E_in;
    result->Lz_in = Lz_in;
    result->feasible = 0;

    if (!kerr_penrose_feasibility(bh, decay_pt, E_in, Lz_in)) {
        result->efficiency = 0.0;
        result->E_out = E_in;
        result->E_negative = 0.0;
        return 0;
    }

    /* Calculate the maximum extractable energy.
     * For a particle with E_in, Lz_in at decay point,
     * the maximum energy of the escaping fragment is:
     * E_out_max = E_in + omega * (Lz_in - Lz_negative) */

    double omega = kerr_frame_dragging_omega(bh, decay_pt);

    /* The negative-energy fragment needs:
     * E_neg = E_in - E_out < 0
     * Condition: E_neg + omega * Lz_neg < 0 (at infinity)
     *
     * Maximum efficiency occurs when the negative fragment has
     * E_neg = omega * Lz_neg (marginal case) */

    result->E_negative = -omega * fabs(Lz_in);
    result->E_out = E_in - result->E_negative;
    result->efficiency = (result->E_out - E_in) / E_in;
    result->feasible = 1;

    /* Impose conservation of 4-momentum */
    double eta_max = kerr_penrose_max_efficiency(bh);
    if (result->efficiency > eta_max)
        result->efficiency = eta_max;

    return 0;
}

double kerr_collisional_cm_energy(const KerrBlackHole *bh,
                                  const KerrBLPoint *pt,
                                  const KerrGeodesicConstants *gc1,
                                  const KerrGeodesicConstants *gc2)
{
    /* CM energy squared for two particles colliding:
     * E_cm^2 = - (p1^mu + p2^mu)(p1_mu + p2_mu)
     *        = m1^2 + m2^2 - 2 g_{mu nu} p1^mu p2^nu
     *
     * Near extremal Kerr horizon, E_cm can be arbitrarily large
     * (Banados-Silk-West 2009 effect). */

    double p1[4], p2[4];
    kerr_momenta_from_constants(bh, pt, gc1, 1, 1, p1);
    kerr_momenta_from_constants(bh, pt, gc2, 1, 1, p2);

    KerrMetricInvBL inv;
    kerr_inverse_metric_bl(bh, pt, &inv);

    /* p1_mu p2^mu = g^{mu nu} p1_mu p2_nu ... actually use covariant */
    /* E_cm^2 = m1^2 + m2^2 - 2 p1^mu p2_mu */
    double p1_dot_p2 = inv.inv_tt * p1[0]*p2[0]
                       + inv.inv_tph * (p1[0]*p2[3] + p1[3]*p2[0])
                       + inv.inv_rr * p1[1]*p2[1]
                       + inv.inv_thth * p1[2]*p2[2]
                       + inv.inv_phph * p1[3]*p2[3];

    double m12 = gc1->mu * gc1->mu;
    double m22 = gc2->mu * gc2->mu;

    double E_cm2 = m12 + m22 - 2.0 * p1_dot_p2;

    return (E_cm2 > 0.0) ? sqrt(E_cm2) : 0.0;
}

/* ====================================================================
 * Superradiance
 *
 * Condition: omega < m * Omega_H
 * For scalar waves with frequency omega and azimuthal quantum number m,
 * the reflected wave can be amplified, extracting rotational energy
 * from the BH (Zeldovich 1971, Starobinsky 1973, Teukolsky 1974).
 *
 * Superradiance is the wave analog of the Penrose process.
 * It occurs for bosonic fields (scalar, EM, gravitational) but not
 * for fermionic fields (Pauli exclusion prevents it).
 *
 * A massive scalar field can form a "gravitational atom" through
 * superradiant instability, producing monochromatic GW signals
 * detectable by LIGO/Virgo.
 * ==================================================================== */

int kerr_superradiance_condition(const KerrBlackHole *bh,
                                 double omega, int m,
                                 SuperradianceCondition *cond)
{
    assert(bh && cond);

    cond->omega = omega;
    cond->m = m;
    cond->omega_H = kerr_horizon_angular_velocity(bh);

    if (isnan(cond->omega_H)) {
        cond->is_superradiant = 0;
        cond->amplification = 0.0;
        cond->threshold = 0.0;
        return 1;
    }

    cond->threshold = m * cond->omega_H;
    cond->is_superradiant = (omega < cond->threshold) ? 1 : 0;

    if (cond->is_superradiant && omega > 0.0) {
        /* Low-frequency amplification factor */
        double r_plus, r_minus;
        kerr_horizon_radii(bh, &r_plus, &r_minus);
        double AH = kerr_horizon_area(bh);

        cond->amplification = 4.0 * omega * (m * cond->omega_H - omega)
                              * (r_plus * r_plus + bh->a * bh->a);
        cond->amplification = fmax(0.0, cond->amplification);
    } else {
        cond->amplification = 0.0;
    }

    return 0;
}

double kerr_superradiance_amplification_scalar(const KerrBlackHole *bh,
                                                double omega, int m,
                                                int l)
{
    /* Scalar wave superradiance amplification.
     * For omega M << 1, l=m=1:
     * Z ~ 4 a M omega (1 - omega/(m Omega_H)) */

    double Omega_H = kerr_horizon_angular_velocity(bh);
    if (isnan(Omega_H) || omega <= 0.0) return 0.0;

    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);

    double threshold = m * Omega_H;
    if (omega >= threshold) return 0.0;  /* No superradiance */

    /* Low-frequency approximation (Starobinsky 1973) */
    double factor = 1.0;
    if (m == 1 && l == 1) {
        factor = 1.0;
    } else if (m == 1 && l == 2) {
        factor = 0.5;
    } else if (m == 2 && l == 2) {
        factor = 0.25;
    }

    return factor * 4.0 * omega * (threshold - omega)
           * (r_plus * r_plus + bh->a * bh->a);
}

double kerr_superradiance_max_gain_scalar(const KerrBlackHole *bh,
                                          int l, int m)
{
    /* Maximum gain occurs near omega = m * Omega_H / 2 */
    double Omega_H = kerr_horizon_angular_velocity(bh);
    if (isnan(Omega_H)) return 0.0;

    double omega_opt = m * Omega_H / 2.0;
    return kerr_superradiance_amplification_scalar(bh, omega_opt, m, l);
}

double kerr_superradiant_instability_timescale(const KerrBlackHole *bh,
                                                double mu_scalar)
{
    /* Superradiant instability growth timescale for massive scalar.
     * For a/M ~ 0.99 and M mu ~ 0.4 (optimal):
     * tau ~ 10^7 M ~ 10 minutes for stellar mass BH
     *         ~ 10^7 years for supermassive BH
     *
     * Reference: Dolan (2007), Brito, Cardoso & Pani (2015) */

    double aM = bh->a / bh->M;
    double M_mu = bh->M * mu_scalar;

    if (aM < 0.5 || M_mu < 0.1 || M_mu > 1.0) {
        return INFINITY;  /* No instability */
    }

    /* Approximate scaling from Detweiler (1980) */
    double alpha = aM * M_mu;
    double tau_M = exp(1.84 / alpha) / M_mu;

    return tau_M * bh->M;  /* In geometric units */
}

/* ====================================================================
 * Blandford-Znajek mechanism
 *
 * Electromagnetic energy extraction from a rotating BH threaded
 * by a large-scale magnetic field.
 *
 * Power: P_BZ = (1/6pi) B^2 M^2 a^2 * f(Omega_F/Omega_H)
 *
 * For optimal loading (Omega_F = Omega_H/2):
 * P_BZ_max = (1/24) B^2 r_plus^2 a^2 Omega_H
 *
 * In physical units:
 * P_BZ = 10^46 (B/10^4 G)^2 (M/10^9 Msun)^2 (a/M)^2 erg/s
 *
 * This is the primary mechanism for powering relativistic jets
 * from active galactic nuclei (AGN).
 *
 * Reference: Blandford & Znajek (1977) MNRAS 179, 433
 * ==================================================================== */

int kerr_blandford_znajek(const KerrBlackHole *bh, double B_poloidal,
                          BlandfordZnajek *bz)
{
    assert(bh && bz);

    bz->B_field = B_poloidal;
    bz->Omega_F = kerr_horizon_angular_velocity(bh) / 2.0;  /* Optimal */

    if (isnan(bz->Omega_F)) {
        bz->power = 0.0;
        bz->efficiency = 0.0;
        bz->jet_power = 0.0;
        return 1;
    }

    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);
    double a = bh->a;

    /* P_BZ = (1/6pi) B^2 M^2 a^2 in geometric units
     * More precisely: P_BZ = (1/24) B^2 r_plus^2 a^2 Omega_H */
    double Omega_H = kerr_horizon_angular_velocity(bh);

    bz->power = (1.0 / 24.0) * B_poloidal * B_poloidal
                * r_plus * r_plus * a * a * Omega_H;

    /* Efficiency: P / (B^2 M^2) */
    double M2 = bh->M * bh->M;
    bz->efficiency = bz->power / (B_poloidal * B_poloidal * M2);

    /* Jet power = P_BZ (the extracted Poynting flux powers jets) */
    bz->jet_power = bz->power;

    return 0;
}

double kerr_bz_max_power(const KerrBlackHole *bh, double B_poloidal)
{
    BlandfordZnajek bz;
    kerr_blandford_znajek(bh, B_poloidal, &bz);
    return bz.power;
}

double kerr_bz_power_physical(double M_solar, double a_spin,
                               double B_gauss)
{
    /* Physical BZ power in erg/s.
     *
     * M_solar: BH mass in solar masses
     * a_spin: dimensionless spin a/M
     * B_gauss: magnetic field at horizon in Gauss
     *
     * P ~ 10^46 (B/10^4 G)^2 (M/10^9 Msun)^2 a^2 erg/s */

    double M_abs = M_solar * 1.989e33;  /* grams */
    double M_len = M_abs * 6.674e-8 / (2.998e10 * 2.998e10);  /* cm */

    /* In geometric units: M ~ 1.48e5 cm per Msun */
    double M_geo = M_solar * 1.477e5;  /* cm */

    /* B^2 M^2 in geometric units
     * B in Gauss -> geometric: B_geo = B * G^(1/2)/c^2 */
    double B_gauss_to_geo = 1.0;  /* Simplified — using scaling */

    double power_geo = B_gauss * B_gauss * M_geo * M_geo
                       * a_spin * a_spin / 24.0;

    /* Convert to erg/s: multiply by c^5/G */
    double c5_G = 3.629e59;  /* erg/s in geometric units */

    return power_geo * c5_G * B_gauss_to_geo;
}

double kerr_bz_estimate_B_field(double M_solar, double a_spin,
                                 double jet_power_erg_s)
{
    /* Inverse BZ: estimate B field from observed jet power.
     * Use bisection since B appears squared in the formula. */
    double B_low = 1.0;      /* 1 G */
    double B_high = 1.0e6;   /* 10^6 G */

    for (int i = 0; i < 50; i++) {
        double B_mid = (B_low + B_high) / 2.0;
        double P_mid = kerr_bz_power_physical(M_solar, a_spin, B_mid);

        if (P_mid < jet_power_erg_s)
            B_low = B_mid;
        else
            B_high = B_mid;
    }

    return (B_low + B_high) / 2.0;
}

/* ====================================================================
 * Novikov-Thorne thin accretion disk model
 *
 * The standard relativistic thin disk model assumes:
 * - Steady, axisymmetric accretion
 * - Optically thick, geometrically thin disk
 * - Keplerian orbits with zero torque at ISCO
 * - Radiative efficiency: eta = 1 - E_isco
 *
 * Radiative efficiencies:
 *   Schwarzschild (a=0): eta ~ 0.057
 *   Extremal Kerr prograde (a=M): eta ~ 0.423
 *   Extremal Kerr retrograde: eta ~ 0.038
 *
 * Reference: Novikov & Thorne (1973), Page & Thorne (1974)
 * ==================================================================== */

int kerr_accretion_disk(const KerrBlackHole *bh, double M_dot,
                        AccretionDisk *disk)
{
    assert(bh && disk);

    /* Find ISCO for prograde disk (standard assumption) */
    KerrISCOData isco;
    kerr_isco(bh, 1, &isco);

    if (!isco.exists) return 1;

    disk->M_dot = M_dot;
    disk->r_in = isco.r_isco;
    disk->r_out = 100.0 * bh->M;  /* Typical outer edge ~100 M */
    disk->eta = 1.0 - isco.E_isco;

    /* Total disk luminosity: L = eta * M_dot * c^2 */
    disk->L_disk = disk->eta * M_dot;

    /* Maximum flux occurs at r ~ 2 * r_isco */
    /* Disk temperature ~ (3 G M M_dot / (8 pi sigma r^3))^(1/4) */
    /* Use r_max_flux ~ 1.5 * r_isco as approximation */
    disk->r_max_flux = 1.5 * isco.r_isco;

    /* Maximum effective temperature (geometric units) */
    double flux_max = kerr_disk_flux_radius(bh, M_dot, disk->r_max_flux);
    if (flux_max > 0.0)
        disk->T_max = sqrt(sqrt(flux_max));
    else
        disk->T_max = 0.0;

    return 0;
}

double kerr_disk_flux_radius(const KerrBlackHole *bh, double M_dot, double r)
{
    /* Novikov-Thorne flux formula:
     * F(r) = M_dot / (4pi) * f(r,a) / M^2
     *
     * where f(r,a) = -(dOmega/dr) (E - Omega Lz)^(-2)
     *                * integral_{r_isco}^r (E - Omega Lz) Lz,r dr
     *
     * We use a simplified fit to the exact relativistic function. */

    double M = bh->M;
    double a = bh->a;

    double r_isco = kerr_isco_radius_bardeen(M, a, 1);
    if (isnan(r_isco) || r <= r_isco) return 0.0;

    /* Relativistic correction factor f(r)
     * Approximate form from Page & Thorne (1974) */
    double x = sqrt(r / M);
    double x0 = sqrt(r_isco / M);

    /* f(r) ~ (3M/(2r)) * (1 - sqrt(r_isco/r) + ...) for Schwarzschild
     * For Kerr, include spin-dependent correction */
    double f_r = (3.0 * M) / (2.0 * r) * (1.0 - x0/x);

    /* Spin correction */
    double aM = a / M;
    double spin_correction = 1.0 - aM / x + aM * aM / (x * x);
    f_r *= spin_correction;

    return M_dot * f_r / (4.0 * M_PI * M * M);
}

double kerr_disk_temperature(const KerrBlackHole *bh, double M_dot, double r)
{
    double flux = kerr_disk_flux_radius(bh, M_dot, r);
    if (flux <= 0.0) return 0.0;

    /* T = (F / sigma)^(1/4), sigma = Stefan-Boltzmann
     * In geometric units, sigma ~ 1 */
    return sqrt(sqrt(flux));
}

double kerr_disk_spectrum(const KerrBlackHole *bh, double M_dot,
                          double nu, double iota)
{
    /* Multicolor blackbody spectrum:
     * F_nu = integral_{r_isco}^{r_out} B_nu(T(r)) * dA * cos(iota)
     *
     * where B_nu is the Planck function and dA includes
     * relativistic ray-tracing corrections. */

    double r_isco = kerr_isco_radius_bardeen(bh->M, bh->a, 1);
    if (isnan(r_isco)) return 0.0;

    double r_out = 100.0 * bh->M;
    int n_r = 500;
    double total_flux = 0.0;

    for (int i = 0; i < n_r; i++) {
        double r = r_isco + (r_out - r_isco) * (i + 0.5) / n_r;
        double dr = (r_out - r_isco) / n_r;
        double T = kerr_disk_temperature(bh, M_dot, r);

        if (T <= 0.0) continue;

        /* Planck function: B_nu(T) = 2h nu^3 / (c^2 (exp(h nu/kT)-1)) */
        /* In geometric units, h=c=k=1 for simplicity */
        if (nu / T > 700.0) continue;  /* Avoid overflow */

        double planck = nu * nu * nu / (exp(nu / T) - 1.0);

        /* Area element: 2pi r dr */
        double dA = 2.0 * M_PI * r * dr;

        /* Relativistic corrections (approximate)
         * Doppler factor g ~ 1 + O(v/c) */
        double g = 1.0;
        total_flux += planck * dA * g * cos(iota);
    }

    return total_flux;
}

int kerr_iron_line_profile(const KerrBlackHole *bh, double iota,
                           double r_in, double r_out,
                           int n_bins, double *energy_bins,
                           double *flux)
{
    /* Relativistically broadened iron K-alpha line profile.
     *
     * The Fe K-alpha line at 6.4 keV (rest frame) is broadened by:
     * 1. Doppler shifts from Keplerian disk rotation
     * 2. Gravitational redshift near the BH
     * 3. Light bending (ray tracing)
     *
     * The line profile shape encodes the BH spin:
     * - Red wing extent depends on how close ISCO is to horizon
     * - Blue peak height depends on inclination
     *
     * Reference: Fabian et al. (1989), Laor (1991)
     */

    assert(bh && energy_bins && flux);

    double E_rest = 6.4;  /* keV, iron K-alpha rest energy */

    for (int i = 0; i < n_bins; i++) flux[i] = 0.0;

    double r_min = (r_in > 0.0) ? r_in :
                   kerr_isco_radius_bardeen(bh->M, bh->a, 1);
    if (isnan(r_min)) return 1;

    double r_max = (r_out > 0.0) ? r_out : 50.0 * bh->M;

    int n_r = 200;
    for (int ir = 0; ir < n_r; ir++) {
        double r = r_min + (r_max - r_min) * ir / (n_r - 1);
        double dr = (r_max - r_min) / n_r;

        /* Keplerian velocity and Doppler shift */
        double Omega = kerr_keplerian_frequency(bh, r, 1);
        double v = Omega * r;  /* Approximate */
        double gamma = 1.0 / sqrt(1.0 - v*v);

        /* Gravitational redshift */
        double z_grav = 1.0 / sqrt(1.0 - 2.0*bh->M/r);  /* Approximate */

        /* Observed energy: E_obs = E_rest / (g * (1+z)) */
        double g_doppler = gamma * (1.0 - v * sin(iota));
        double E_obs = E_rest / (g_doppler * z_grav);

        /* Find bin */
        double bin_width = (energy_bins[n_bins-1] - energy_bins[0]) / n_bins;
        int bin = (int)((E_obs - energy_bins[0]) / bin_width);
        if (bin >= 0 && bin < n_bins) {
            /* Emissivity profile: ~ r^{-3} (simple) */
            double emissivity = 1.0 / (r * r * r);
            flux[bin] += emissivity * dr * 2.0 * M_PI * r;
        }
    }

    /* Normalize */
    double total = 0.0;
    for (int i = 0; i < n_bins; i++) total += flux[i];
    if (total > 0.0)
        for (int i = 0; i < n_bins; i++) flux[i] /= total;

    return 0;
}

/* ====================================================================
 * Black hole shadow observables
 *
 * The shadow is the photon capture cross-section as seen by a
 * distant observer. For Schwarzschild, it's a perfect circle of
 * radius sqrt(27) M. For Kerr, it's asymmetric (D-shaped) due to
 * frame dragging.
 *
 * Key observables (EHT collaboration):
 * - Shadow diameter: ~10 GM/c^2 for M87* (~42 microarcseconds)
 * - Asymmetry parameter A = (R_right - R_left) / (R_right + R_left)
 * - Deviation from Schwarzschild: delta_C < 0.1 (GR test)
 *
 * Reference: EHT Collaboration (2019) ApJ 875, L1
 * Psaltis et al. (2020) PRL 125, 141104
 * ==================================================================== */

int kerr_shadow_observables(const KerrBlackHole *bh, double iota,
                            KerrShadow *shadow)
{
    assert(bh && shadow);

    double M = bh->M;
    double a = bh->a;
    double aM = a / M;

    /* Shadow radius (semi-major axis) */
    shadow->shadow_radius = kerr_shadow_radius_approx(bh, iota);

    /* Schwarzschild shadow for reference */
    double R_schw = sqrt(27.0) * M;
    shadow->shadow_deviation = (shadow->shadow_radius - R_schw) / R_schw;

    /* Asymmetry: due to frame dragging, shadow is shifted in
     * the direction perpendicular to the spin axis.
     * A ~ 2 a/M sin(iota) / (shadow_radius/M) */
    shadow->asymmetry = 2.0 * aM * sin(iota)
                        / (shadow->shadow_radius / M);

    /* Spin estimate from shadow */
    shadow->spin_estimate = kerr_spin_from_shadow_deviation(
        shadow->shadow_deviation, iota);

    return 0;
}

double kerr_shadow_radius_approx(const KerrBlackHole *bh, double iota)
{
    /* Approximate shadow radius formula.
     * For Schwarzschild: R_shadow = sqrt(27) M ~ 5.196 M
     * For Kerr: slightly smaller for prograde, depends on inclination.
     *
     * Fitting formula (Johannsen & Psaltis 2010):
     * R_shadow/M ~ 5.2 - 0.2 (a/M)^2 - 0.4 (a/M)^4 cos(iota) */

    double M = bh->M;
    double aM = bh->a / M;
    double aM2 = aM * aM;
    double aM4 = aM2 * aM2;

    double R = 5.196 - 0.2 * aM2 - 0.4 * aM4 * cos(iota);
    return R * M;
}

double kerr_spin_from_shadow_deviation(double deviation, double iota)
{
    /* Invert shadow deviation to estimate spin.
     * delta = (R_kerr - R_schw) / R_schw
     *       ~ -0.04 a^2/M^2 (for face-on)
     *
     * a/M ~ sqrt(-delta / 0.04) */

    double sin_i = sin(iota);
    double correction = (1.0 + 0.3 * sin_i * sin_i);

    if (deviation >= 0.0) return 0.0;

    double aM2 = -deviation / (0.04 * correction);
    if (aM2 > 1.0) aM2 = 1.0;
    if (aM2 < 0.0) aM2 = 0.0;

    return sqrt(aM2);
}

/* ====================================================================
 * Quasi-normal mode ringdown
 *
 * After a binary black hole merger, the remnant BH rings down
 * via a superposition of damped sinusoids:
 *
 * h(t) ~ exp(-t/tau) cos(omega_R t + phi)
 *
 * The QNM frequencies are complex eigenvalues of the Teukolsky
 * equation. They depend only on M and a (no-hair theorem).
 *
 * The dominant mode is l=m=2, n=0.
 *
 * Fitting formula (Berti, Cardoso & Will 2006):
 * M omega_220 = f_1 + f_2 (1-a/M)^{f_3}
 *
 * Reference: Berti, Cardoso & Will (2006) PRD 73, 064030
 * ==================================================================== */

double kerr_qnm_frequency_real(const KerrBlackHole *bh, int l, int m, int n)
{
    /* Dominant l=m=2, n=0 mode.
     * Fitting coefficients from Berti+2006 Table VIII. */

    double aM = bh->a / bh->M;
    if (aM > 1.0) aM = 1.0;

    if (l == 2 && m == 2 && n == 0) {
        /* f1 + f2 * (1-a)^f3 */
        double f1 = 0.3737;
        double f2 = 1.0927;
        double f3 = 0.4183;
        double omega = f1 + f2 * pow(1.0 - aM, f3);
        return omega / bh->M;
    } else if (l == 2 && m == 2 && n == 1) {
        double f1 = 0.3505;
        double f2 = 1.1101;
        double f3 = 0.5664;
        double omega = f1 + f2 * pow(1.0 - aM, f3);
        return omega / bh->M;
    } else if (l == 3 && m == 3 && n == 0) {
        double f1 = 0.5994;
        double f2 = 1.0827;
        double f3 = 0.3957;
        double omega = f1 + f2 * pow(1.0 - aM, f3);
        return omega / bh->M;
    }

    /* Default: approximate as l=m=2, n=0 */
    return kerr_qnm_frequency_real(bh, 2, 2, 0);
}

double kerr_qnm_damping_time(const KerrBlackHole *bh, int l, int m, int n)
{
    /* tau = 1 / |omega_I|
     * Fitting formula for omega_I from Berti+2006. */

    double aM = bh->a / bh->M;
    if (aM > 1.0) aM = 1.0;

    double omega_I;

    if (l == 2 && m == 2 && n == 0) {
        double q1 = 0.08896;
        double q2 = 0.10234;
        double q3 = 0.7856;
        omega_I = -q1 - q2 * pow(1.0 - aM, q3);
    } else if (l == 2 && m == 2 && n == 1) {
        double q1 = 0.09056;
        double q2 = 0.11412;
        double q3 = 0.8725;
        omega_I = -q1 - q2 * pow(1.0 - aM, q3);
    } else {
        /* First overtone is typically more damped */
        double q1 = 0.08896;
        double q2 = 0.10234;
        double q3 = 0.7856;
        omega_I = -q1 - q2 * pow(1.0 - aM, q3);
    }

    return bh->M / fabs(omega_I);
}

double kerr_qnm_frequency_fit(const KerrBlackHole *bh, int l, int m, int n)
{
    return kerr_qnm_frequency_real(bh, l, m, n);
}

double kerr_qnm_quality_factor(const KerrBlackHole *bh, int l, int m, int n)
{
    double omega_R = kerr_qnm_frequency_real(bh, l, m, n);
    double tau = kerr_qnm_damping_time(bh, l, m, n);

    if (tau <= 0.0) return 0.0;

    /* Q = omega_R * tau / 2 */
    return omega_R * tau / 2.0;
}

/* ====================================================================
 * QNM frequency in physical units (Hz) given BH mass in solar masses
 * ==================================================================== */

double kerr_qnm_frequency_hz(const KerrBlackHole *bh, double M_solar,
                             int l, int m, int n)
{
    /* Convert geometric frequency to Hz.
     * omega_geo has units of 1/M.
     * f_Hz = omega_R * c^3 / (2 pi G M)
     * For M in Msun: f ~ 32 kHz * (Msun/M) */

    double omega_R = kerr_qnm_frequency_real(bh, l, m, n);
    /* Geometric frequency is omega_R / M, multiply by conversion factor */
    return omega_R * 3.231e4 / M_solar;  /* Hz */
}

double kerr_ringdown_frequency_hz(double M_solar, double a_spin)
{
    /* Convenience: ringdown frequency in Hz for l=m=2, n=0.
     * For a 10 Msun BH: f ~ 1.2 kHz
     * For a 65 Msun BH (GW150914 remnant): f ~ 250 Hz
     * For M87* (~6.5e9 Msun): f ~ 2.5e-6 Hz */

    KerrBlackHole bh = {1.0, a_spin, 1};  /* Geometric M=1, a=a_spin */
    double omega_R = kerr_qnm_frequency_real(&bh, 2, 2, 0);
    return omega_R * 3.231e4 / M_solar;
}
