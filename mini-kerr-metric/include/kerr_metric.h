/**
 * @file kerr_metric.h
 * @brief Kerr metric components in Boyer-Lindquist and Kerr-Schild coordinates
 *
 * The Kerr metric (1963) is the unique stationary, axisymmetric, asymptotically
 * flat vacuum solution to Einstein equations describing a rotating black hole.
 *
 * Reference: R. P. Kerr, Phys. Rev. Lett. 11, 237 (1963).
 * MIT 8.962 (General Relativity), Wald Ch.7, Chandrasekhar Ch.6
 */

#ifndef KERR_METRIC_H
#define KERR_METRIC_H

#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1 -- Core Definitions
 * ========================================================================== */

#define KERR_G  (6.67430e-11)
#define KERR_C  (299792458.0)

/** Kerr black hole parameters */
typedef struct {
    double M;        /* Mass parameter */
    double a;        /* Spin parameter a = J/M */
    int geo_units;   /* 1 = geometric units (G=c=1), 0 = SI */
} KerrBlackHole;

/** Boyer-Lindquist coordinate point */
typedef struct {
    double t;        /* Coordinate time */
    double r;        /* Radial coordinate */
    double theta;    /* Polar angle in [0, pi] */
    double phi;      /* Azimuthal angle in [0, 2pi) */
} KerrBLPoint;

/** Kerr-Schild (Cartesian-like) coordinate point */
typedef struct {
    double t;
    double x, y, z;
} KerrKSPoint;

/** 4x4 metric tensor g_munu in Boyer-Lindquist coordinates */
typedef struct {
    double g_tt, g_tr, g_tth, g_tph;
    double g_rt, g_rr, g_rth, g_rph;
    double g_tht, g_thr, g_thth, g_thph;
    double g_pht, g_phr, g_phth, g_phph;
} KerrMetricBL;

/** Inverse metric g^munu in Boyer-Lindquist coordinates */
typedef struct {
    double inv_tt, inv_tr, inv_tth, inv_tph;
    double inv_rt, inv_rr, inv_rth, inv_rph;
    double inv_tht, inv_thr, inv_thth, inv_thph;
    double inv_pht, inv_phr, inv_phth, inv_phph;
} KerrMetricInvBL;

/** Petrov classification of the Weyl tensor */
typedef enum {
    PETROV_I = 0, PETROV_II = 1, PETROV_D = 2,
    PETROV_III = 3, PETROV_N = 4, PETROV_O = 5
} PetrovType;

/** Weyl scalars in Newman-Penrose formalism */
typedef struct {
    double Psi0, Psi1, Psi2, Psi3, Psi4;
} WeylScalars;

/** Newman-Penrose spin coefficients */
typedef struct {
    double kappa, sigma, rho, epsilon;
    double pi, alpha, beta, gamma;
    double tau, mu, nu, lambda;
} KerrSpinCoefficients;

/** Tetrad vectors for ZAMO (Zero Angular Momentum Observer) */
typedef struct {
    double e_t[4];  /* Timelike tetrad leg */
    double e_r[4];  /* Radial leg */
    double e_th[4]; /* Theta leg */
    double e_ph[4]; /* Phi leg */
} KerrZAMOTetrad;

/* ==========================================================================
 * L2 -- Core Metric Functions
 * ========================================================================== */

/**
 * Sigma = r^2 + a^2 cos^2(theta).
 * Fundamental metric function controlling oblateness of Kerr geometry.
 */
double kerr_sigma(double r, double theta, double a);

/**
 * Delta = r^2 - 2Mr + a^2.
 * Horizon function: Delta(r_plus) = Delta(r_minus) = 0.
 */
double kerr_delta(double r, double M, double a);

/**
 * A = (r^2 + a^2)^2 - a^2 Delta sin^2(theta).
 * Appears in g_phiphi component and frame-dragging formulas.
 */
double kerr_A_function(double r, double theta, double M, double a);

/**
 * omega_tilde = 2 M a r / Sigma.
 * Shorthand for the gtphi metric component.
 */
double kerr_omega_tilde(double r, double theta, double M, double a);

/** Compute the full Boyer-Lindquist metric tensor at a point.
 *  Line element:
 *  ds^2 = -(1-2Mr/Sigma)dt^2 - (4Mar sin^2theta/Sigma)dt dphi
 *         + (Sigma/Delta)dr^2 + Sigma dtheta^2
 *         + (r^2+a^2+2Ma^2r sin^2theta/Sigma) sin^2theta dphi^2 */
int kerr_metric_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                   KerrMetricBL *metric);

/** Kerr metric in Kerr-Schild Cartesian coordinates.
 *  g_munu = eta_munu + f k_mu k_nu  where k_mu is null.
 *  Regular at horizons. */
int kerr_metric_ks(const KerrBlackHole *bh, const KerrKSPoint *pt,
                   double g[4][4]);

/** Inverse Boyer-Lindquist metric g^munu */
int kerr_inverse_metric_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                           KerrMetricInvBL *inv);

/** Metric determinant in BL coordinates: det(g) = -Sigma^2 sin^2(theta) */
double kerr_metric_determinant_bl(const KerrBlackHole *bh,
                                  const KerrBLPoint *pt);

/** 3+1 ADM decomposition: lapse function N for Kerr */
double kerr_lapse_bl(const KerrBlackHole *bh, const KerrBLPoint *pt);

/** 3+1 ADM decomposition: shift vector beta^i for Kerr */
int kerr_shift_vector_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                         double beta[3]);

/** 3+1 ADM: spatial 3-metric gamma_ij induced on t=const hypersurface */
int kerr_three_metric_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                         double gamma[3][3]);

/** 3+1 ADM: extrinsic curvature K_ij for Kerr */
int kerr_extrinsic_curvature_bl(const KerrBlackHole *bh,
                                const KerrBLPoint *pt, double K[3][3]);

/* ==========================================================================
 * L3 -- Tensor Calculus
 * ========================================================================== */

/**
 * Compute Christoffel symbols Gamma^mu_nurho for Kerr in BL coords.
 * Stored as christoffel[16*mu + 4*nu + rho] = Gamma^mu_nurho.
 * 20 non-zero independent components (Chandrasekhar Ch.6).
 */
int kerr_christoffel_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                        double christoffel[64]);

/** Riemann curvature tensor R^mu_nurhosigma, flattened to 256 doubles */
int kerr_riemann_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                    double riemann[256]);

/** Ricci tensor R_munu (should be identically zero for vacuum Kerr) */
int kerr_ricci_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                  double ricci[16]);

/** Einstein tensor G_munu = R_munu - (1/2)R g_munu (zero for vacuum) */
int kerr_einstein_bl(const KerrBlackHole *bh, const KerrBLPoint *pt,
                     double einstein[16]);

/**
 * Kretschmann scalar K = R_{munurhosigma} R^{munurhosigma}.
 * For Kerr: K = 48M^2(r^6 - 15a^2 r^4 cos^2theta + 15a^4 r^2 cos^4theta
 *             - a^6 cos^6theta) / (r^2 + a^2 cos^2theta)^6
 */
double kerr_kretschmann(const KerrBlackHole *bh, const KerrBLPoint *pt);

/** Weyl scalars in Kinnersley tetrad (only Psi2 non-zero for Type D) */
int kerr_weyl_scalars(const KerrBlackHole *bh, const KerrBLPoint *pt,
                      WeylScalars *psi);

/** Christoffel symbols via numerical differentiation (central diff, O(h^2)) */
int kerr_christoffel_numerical(const KerrBlackHole *bh,
                               const KerrBLPoint *pt, double h,
                               double christoffel[64]);

/** Compute null tetrad legs in Kinnersley tetrad */
int kerr_kinnersley_tetrad(const KerrBlackHole *bh, const KerrBLPoint *pt,
                           double l[4], double n[4], double m[4],
                           double mbar[4]);

/** Compute spin coefficients in Kinnersley tetrad (Chandrasekhar Ch.7) */
int kerr_spin_coefficients(const KerrBlackHole *bh, const KerrBLPoint *pt,
                           KerrSpinCoefficients *sc);

/* ==========================================================================
 * Coordinate transformations
 * ========================================================================== */

/** BL to KS: x+iy = (r+ia) sin(theta) e^{iphi}, z = r cos(theta) */
int kerr_bl_to_ks(const KerrBLPoint *bl, double a, KerrKSPoint *ks);

/** KS to BL: solve quartic for r */
int kerr_ks_to_bl(const KerrKSPoint *ks, double a, KerrBLPoint *bl);

/** Ingoing Eddington-Finkelstein v coordinate for Kerr */
double kerr_ingress_v(const KerrBLPoint *pt, const KerrBlackHole *bh);

/** Check if point is at a coordinate singularity (Delta=0) */
int kerr_is_coord_singularity(const KerrBLPoint *pt, const KerrBlackHole *bh);

/** Compute ZAMO (Zero Angular Momentum Observer) tetrad */
int kerr_zamo_tetrad(const KerrBlackHole *bh, const KerrBLPoint *pt,
                     KerrZAMOTetrad *tetrad);

/** Compute proper time along a ZAMO worldline between two t coordinates */
double kerr_zamo_proper_time(const KerrBlackHole *bh, const KerrBLPoint *pt,
                             double dt);

/** Compute the locally non-rotating frame angular velocity Omega_ZAMO */
double kerr_zamo_angular_velocity(const KerrBlackHole *bh,
                                  const KerrBLPoint *pt);

#endif /* KERR_METRIC_H */