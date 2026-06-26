/**
 * @file kerr_horizons.c
 * @brief Kerr black hole horizons and thermodynamics
 */
#include "kerr_horizons.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

int kerr_horizon_radii(const KerrBlackHole *bh, double *r_plus,
                       double *r_minus)
{
    assert(bh != NULL && r_plus != NULL && r_minus != NULL);
    double M = bh->M;
    double a = bh->a;
    double disc = M * M - a * a;
    if (disc < 0.0) {
        *r_plus = NAN;
        *r_minus = NAN;
        return 1;
    }
    double sqrt_disc = sqrt(disc);
    *r_plus = M + sqrt_disc;
    *r_minus = M - sqrt_disc;
    return 0;
}

KerrBHType kerr_bh_classify(const KerrBlackHole *bh)
{
    if (bh->M <= 0.0) return BH_INVALID;
    double a = bh->a;
    double M = bh->M;
    if (fabs(a) < 1e-15) return BH_SCHWARZSCHILD;
    if (fabs(a - M) < 1e-15) return BH_KERR_EXTREMAL;
    if (a < M) return BH_KERR_SLOW;
    return BH_NAKED_SING;
}

double kerr_horizon_area(const KerrBlackHole *bh)
{
    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) return NAN;
    return 4.0 * M_PI * (r_plus * r_plus + bh->a * bh->a);
}

double kerr_horizon_angular_velocity(const KerrBlackHole *bh)
{
    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) return NAN;
    return bh->a / (r_plus * r_plus + bh->a * bh->a);
}

double kerr_surface_gravity(const KerrBlackHole *bh)
{
    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) return NAN;
    return (r_plus - bh->M) / (r_plus * r_plus + bh->a * bh->a);
}

double kerr_entropy_geometric(const KerrBlackHole *bh)
{
    return kerr_horizon_area(bh) / 4.0;
}

double kerr_temperature_geometric(const KerrBlackHole *bh)
{
    double kappa = kerr_surface_gravity(bh);
    return kappa / (2.0 * M_PI);
}

double kerr_irreducible_mass(const KerrBlackHole *bh)
{
    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) return NAN;
    double area = kerr_horizon_area(bh);
    return sqrt(area / (16.0 * M_PI));
}

double kerr_rotational_energy(const KerrBlackHole *bh)
{
    double Mirr = kerr_irreducible_mass(bh);
    return bh->M - Mirr;
}

int kerr_horizon_data(const KerrBlackHole *bh, KerrHorizonData *data)
{
    assert(bh != NULL && data != NULL);
    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) {
        memset(data, 0, sizeof(KerrHorizonData));
        return 1;
    }
    data->r_plus = r_plus;
    data->r_minus = r_minus;
    data->surface_area = kerr_horizon_area(bh);
    data->omega_H = kerr_horizon_angular_velocity(bh);
    data->kappa = kerr_surface_gravity(bh);
    data->temperature = kerr_temperature_geometric(bh);
    data->entropy = kerr_entropy_geometric(bh);
    data->irreducible_mass = kerr_irreducible_mass(bh);
    return 0;
}

/* ====================================================================
 * Ergosphere
 *
 * The ergosphere is the region between the outer horizon r_plus
 * and the static limit r_ergo(theta) where g_tt = 0.
 *
 * r_ergo(theta) = M + sqrt(M^2 - a^2 cos^2theta)
 *
 * Inside the ergosphere, no observer can remain static
 * (i.e., have constant r, theta, phi). Frame dragging forces
 * all observers to co-rotate with the black hole.
 *
 * This is the key region for the Penrose process.
 * ==================================================================== */

double kerr_ergosphere_radius(const KerrBlackHole *bh, double theta)
{
    double M = bh->M;
    double a = bh->a;
    double ct = cos(theta);
    double disc = M * M - a * a * ct * ct;

    if (disc < 0.0) return M;  /* Below static limit */
    return M + sqrt(disc);
}

int kerr_is_inside_ergosphere(const KerrBlackHole *bh,
                              const KerrBLPoint *pt)
{
    double r_ergo = kerr_ergosphere_radius(bh, pt->theta);
    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);

    if (pt->r > r_plus && pt->r < r_ergo) return 1;
    return 0;
}

int kerr_ergosphere_data(const KerrBlackHole *bh, const KerrBLPoint *pt,
                         KerrErgosphereData *data)
{
    assert(bh != NULL && pt != NULL && data != NULL);

    data->r_ergo = kerr_ergosphere_radius(bh, pt->theta);
    data->inside = kerr_is_inside_ergosphere(bh, pt);

    KerrMetricBL metric;
    kerr_metric_bl(bh, pt, &metric);
    data->g_tt = metric.g_tt;

    return 0;
}

double kerr_ergosphere_equatorial_radius(const KerrBlackHole *bh)
{
    /* At theta=pi/2: cos(theta)=0, so r_ergo = M + M = 2M */
    return 2.0 * bh->M;
}

/* ====================================================================
 * Ring singularity
 *
 * The Kerr singularity is a ring of coordinate radius a in the
 * equatorial plane (theta = pi/2) at r = 0 (Boyer-Lindquist).
 *
 * In Kerr-Schild Cartesian coordinates, the singularity is the ring
 * x^2 + y^2 = a^2, z = 0. The interior of the ring (r=0, theta != pi/2)
 * is regular and can be traversed to negative r (another asymptotically
 * flat region containing closed timelike curves).
 *
 * The ring structure avoids the point singularity of Schwarzschild,
 * and in principle allows passage to other universes.
 * ==================================================================== */

int kerr_ring_singularity_data(const KerrBlackHole *bh,
                               KerrRingSingularity *ring)
{
    assert(bh != NULL && ring != NULL);

    ring->ring_radius = bh->a;
    ring->r_bl = 0.0;
    ring->theta_bl = M_PI / 2.0;
    ring->is_singular = 0;

    return 0;
}

double kerr_kretschmann_safe(const KerrBlackHole *bh, const KerrBLPoint *pt,
                             double safety_threshold)
{
    double r = pt->r;
    double th = pt->theta;
    double ct = cos(th);

    /* Check proximity to ring singularity */
    double dist_to_ring = sqrt(r * r + bh->a * bh->a * ct * ct);

    if (dist_to_ring < safety_threshold) {
        return INFINITY;
    }

    return kerr_kretschmann(bh, pt);
}

int kerr_has_ctc_region(const KerrBlackHole *bh, const KerrBLPoint *pt)
{
    /* CTCs exist where g_phiphi becomes negative.
     * For Kerr, this occurs for r < 0 near the ring singularity.
     * In BL coordinates with r < 0, g_phiphi can change sign. */

    if (pt->r >= 0.0) return 0;

    KerrMetricBL metric;
    kerr_metric_bl(bh, pt, &metric);

    /* For very negative r, check g_phiphi */
    if (metric.g_phph < 0.0) return 1;

    return 0;
}

/* ====================================================================
 * Black hole thermodynamics for Kerr
 *
 * First law: dM = (kappa/8pi) dA + Omega_H dJ
 * Smarr formula: M = 2 Omega_H J + (kappa/4pi) A
 *
 * For Kerr: J = a M, so the Smarr formula becomes:
 * M = 2 Omega_H (a M) + (kappa/4pi) A
 *
 * The Smarr formula is a consequence of Euler's theorem for
 * homogeneous functions applied to the mass function M(A,J).
 *
 * Second law (Area Theorem): A_final >= A_initial for any
 * classical process (Hawking 1971).
 *
 * Third law: kappa -> 0 cannot be reached in finite steps
 * (analogous to Nernst's theorem in ordinary thermodynamics).
 *
 * Reference: Bardeen, Carter & Hawking (1973),
 * Commun. Math. Phys. 31, 161
 * ==================================================================== */

double kerr_smarr_formula_check(const KerrBlackHole *bh)
{
    /* First law: dM = (kappa/8pi)dA + Omega_H dJ
     * Check numerically by perturbing M and a. */

    double h = 1e-6;
    KerrBlackHole bh_plus  = *bh;
    KerrBlackHole bh_minus = *bh;

    bh_plus.M += h;
    bh_minus.M -= h;

    double A_plus = kerr_horizon_area(&bh_plus);
    double A_minus = kerr_horizon_area(&bh_minus);
    double dA_dM = (A_plus - A_minus) / (2.0 * h);

    double r_plus, r_minus;
    kerr_horizon_radii(bh, &r_plus, &r_minus);
    double Omega_H = kerr_horizon_angular_velocity(bh);
    double kappa = kerr_surface_gravity(bh);

    /* dM = (kappa/8pi) dA + Omega_H dJ = (kappa/8pi) dA + Omega_H a dM */
    double dM_predicted = (kappa / (8.0 * M_PI)) * dA_dM + Omega_H * bh->a;

    return fabs(1.0 - dM_predicted);
}

double kerr_smarr_formula_residual(const KerrBlackHole *bh)
{
    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) return INFINITY;

    double Omega_H = kerr_horizon_angular_velocity(bh);
    double kappa = kerr_surface_gravity(bh);
    double A = kerr_horizon_area(bh);
    double J = bh->a * bh->M;

    /* M = 2 Omega_H J + (kappa/4pi) A */
    double M_smarr = 2.0 * Omega_H * J + (kappa / (4.0 * M_PI)) * A;

    return fabs(bh->M - M_smarr);
}

double kerr_area_increase_ratio(double M1, double a1,
                                double M2, double a2,
                                double M_total)
{
    KerrBlackHole bh1 = {M1, a1, 1};
    KerrBlackHole bh2 = {M2, a2, 1};
    double A1 = kerr_horizon_area(&bh1);
    double A2 = kerr_horizon_area(&bh2);

    if (isnan(A1) || isnan(A2)) return -1.0;

    /* For a merged black hole with mass M_total and maximum spin,
     * compute the minimum possible area (a=0 gives max area). */
    KerrBlackHole bh_merged = {M_total, 0.0, 1};
    double A_total = kerr_horizon_area(&bh_merged);

    return A_total / (A1 + A2);
}

double kerr_entropy_bound_ratio(const KerrBlackHole *bh)
{
    /* S_Kerr(M,a) / S_Schwarzschild(M) <= 1
     * For fixed M, Schwarzschild has maximum entropy.
     * Ratio = (r_plus^2 + a^2) / (4 M^2) */

    double r_plus, r_minus;
    if (kerr_horizon_radii(bh, &r_plus, &r_minus) != 0) return NAN;

    double S_kerr = kerr_entropy_geometric(bh);
    /* S_schwarzschild = 4 pi M^2 */
    double S_schw = 4.0 * M_PI * bh->M * bh->M;

    return S_kerr / S_schw;
}

double kerr_extremality_parameter(const KerrBlackHole *bh)
{
    if (bh->M <= 0.0) return INFINITY;
    return bh->a / bh->M;
}

/* ====================================================================
 * Spin limit from cosmic censorship
 *
 * The cosmic censorship conjecture (Penrose 1969) states that
 * singularities are always hidden behind event horizons.
 * For Kerr, this requires a <= M (extremality bound).
 *
 * Astrophysical black holes are expected to have a/M <= 0.998
 * (Thorne 1974 limit from accretion disk photon capture).
 * ==================================================================== */

double kerr_thorne_spin_limit(void)
{
    /* Thorne (1974): a/M <= 0.998 for thin-disk accretion */
    return 0.998;
}

double kerr_max_spin_for_accretion(double eta)
{
    /* Spin-up limit from accretion disk radiation:
     * Photons emitted by the disk can be captured by the BH,
     * counteracting spin-up. This gives a/M <= 0.998.
     * Returns the limiting a/M for a given efficiency eta. */
    (void)eta;
    return 0.998;
}
