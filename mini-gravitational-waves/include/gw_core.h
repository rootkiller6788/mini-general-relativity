/**
 * @file    gw_core.h
 * @brief   Gravitational Waves — Core Definitions & Physical Constants
 *
 * Reference: Wald (1984) General Relativity, Ch.4
 *            Maggiore (2008) Gravitational Waves, Vol.1
 *            MIT 8.962 — General Relativity (Prof. Scott Hughes)
 *
 * L1 — Definitions: Physical constants (G,c,Msun,Mpc),
 *      wave tensor h_{mu nu}, TT gauge, quadrupole moment I_{ij},
 *      reduced quadrupole Q_{ij}, strain amplitude h,
 *      polarization basis e^+/e^×
 *
 * L2 — Core Concepts: Linearized gravity on Minkowski background,
 *      gauge freedom and TT projection, GW luminosity / energy flux
 *
 * L3 — Mathematical Structures: 3x3 symmetric trace-free tensor
 *      operations, TT projector Lambda_{ij,kl}
 */

#ifndef GW_CORE_H
#define GW_CORE_H

#include <stddef.h>

/* ================================================================
 * L1 — Physical Constants (CODATA 2018 / IAU 2015)
 * ================================================================ */

/** Newton's gravitational constant [m^3 kg^{-1} s^{-2}] */
#define GW_G      6.67430e-11

/** Speed of light in vacuum [m/s] */
#define GW_C      299792458.0

/** Solar mass [kg] (IAU 2015 nominal) */
#define GW_MSUN   1.98847e30

/** Megaparsec [m] */
#define GW_MPC    3.08567758149137e22

/** G/c^4 [s^2/(kg*m)] — appears in quadrupole formula h ~ (G/c^4) */
#define GW_G_OVER_C4  (GW_G / (GW_C * GW_C * GW_C * GW_C))

/** G/(5 c^5) — appears in quadrupole luminosity formula */
#define GW_G_OVER_5C5 (GW_G / (5.0 * GW_C * GW_C * GW_C * GW_C * GW_C))

/** pi */
#define GW_PI      3.14159265358979323846

/** pi^(2/3) — recurring in GW formulas */
#define GW_PI_2_3  2.1450293971110257

/** pi^(8/3) — recurring in PN evolution */
#define GW_PI_8_3  21.2445034929768

/* ================================================================
 * L1 — 3x3 Symmetric Tensor
 * ================================================================ */

/** 3x3 symmetric matrix for spatial tensors (h_{ij}, I_{ij}, Q_{ij}) */
typedef struct {
    double xx, xy, xz;
    double yx, yy, yz;
    double zx, zy, zz;
} GwTensor3;

/* ================================================================
 * L1 — Polarization Modes
 * ================================================================ */

/** GW polarization modes.
 *  GR predicts only + and x; additional modes test metric theories. */
typedef enum {
    GW_POL_PLUS,        /**< Plus (+) — GR */
    GW_POL_CROSS,       /**< Cross (x) — GR */
    GW_POL_SCALAR_B,    /**< Scalar breathing (Brans-Dicke etc.) */
    GW_POL_SCALAR_L,    /**< Scalar longitudinal */
    GW_POL_VECTOR_X,    /**< Vector-x mode */
    GW_POL_VECTOR_Y     /**< Vector-y mode */
} GwPolarization;

/* ================================================================
 * L1 — Binary Source Parameters
 * ================================================================ */

/** Compact binary coalescence (CBC) source parameters */
typedef struct {
    double m1;           /**< Primary mass [kg] */
    double m2;           /**< Secondary mass [kg] */
    double M_total;      /**< Total mass M = m1 + m2 [kg] */
    double chirp_mass;   /**< Chirp mass Mc [kg] */
    double eta;          /**< Symmetric mass ratio eta = m1*m2/(m1+m2)^2 */
    double D_L;          /**< Luminosity distance [m] */
    double iota;         /**< Inclination angle [rad] */
    double phi_ref;      /**< Reference phase [rad] */
    double f_ref;        /**< Reference frequency [Hz] */
    double t_c;          /**< Coalescence time [s] */
    double phi_c;        /**< Coalescence phase [rad] */
} GwBinaryParams;

/** Single-detector response */
typedef struct {
    double F_plus;       /**< Antenna pattern for + */
    double F_cross;      /**< Antenna pattern for x */
    double ra;           /**< Right ascension [rad] */
    double dec;          /**< Declination [rad] */
    double psi;          /**< Polarization angle [rad] */
    double gmst;         /**< Greenwich mean sidereal time [rad] */
} GwDetectorResponse;

/* ================================================================
 * L1 — GW Strain Data Structures
 * ================================================================ */

/** Time-domain GW strain */
typedef struct {
    size_t   n;          /**< Number of samples */
    double  *t;          /**< Time array [s] */
    double  *hp;         /**< h_plus(t) strain */
    double  *hx;         /**< h_cross(t) strain */
} GwStrainSeries;

/** Frequency-domain GW data */
typedef struct {
    size_t   n;          /**< Number of frequency bins */
    double  *f;          /**< Frequency array [Hz] */
    double  *re;         /**< Real part of strain */
    double  *im;         /**< Imaginary part of strain */
    double  *psd;        /**< One-sided PSD [Hz^{-1}] */
} GwFreqSeries;

/** Stochastic background parameters */
typedef struct {
    double Omega_ref;    /**< Omega_GW at reference frequency */
    double f_ref;        /**< Reference frequency [Hz] */
    double alpha;        /**< Spectral index: Omega(f) ~ f^alpha */
    double f_min;        /**< Minimum frequency [Hz] */
    double f_max;        /**< Maximum frequency [Hz] */
} GwStochasticParams;

/* ================================================================
 * L2 — Tensor Operations
 * ================================================================ */

void gw_tensor_zero(GwTensor3 *T);
double gw_tensor_trace(const GwTensor3 *T);
void gw_tensor_add(GwTensor3 *C, const GwTensor3 *A, const GwTensor3 *B);
void gw_tensor_sub(GwTensor3 *C, const GwTensor3 *A, const GwTensor3 *B);
void gw_tensor_scale(GwTensor3 *T, double alpha);
double gw_tensor_contract(const GwTensor3 *A, const GwTensor3 *B);

/**
 * TT gauge projection: h_{ij}^{TT} = Lambda_{ij,kl} S_{kl}
 *
 * Lambda_{ij,kl} = P_{ik} P_{jl} - (1/2) P_{ij} P_{kl}
 * P_{ij} = delta_{ij} - n_i n_j  (transverse projector)
 *
 * @param h_tt Output TT tensor
 * @param S    Input symmetric tensor
 * @param nx,ny,nz  Unit propagation direction
 */
void gw_tensor_tt_project(GwTensor3 *h_tt, const GwTensor3 *S,
                          double nx, double ny, double nz);

/**
 * Compute eigenvalues and eigenvectors of a 3x3 symmetric tensor.
 * Uses analytical cubic formula for 3x3 symmetric matrix.
 *
 * @param T          Input symmetric tensor
 * @param eigenvals  Output eigenvalues [3], sorted descending
 * @param eigenvecs  Output eigenvectors as 9 doubles (column-major 3x3)
 * @return 0 on success, -1 if algorithm fails
 */
int gw_tensor_eigensystem(const GwTensor3 *T,
                          double eigenvals[3], double eigenvecs[9]);

/* ================================================================
 * L2 — GW Strain Operations
 * ================================================================ */

int  gw_strain_alloc(GwStrainSeries *s, size_t n);
void gw_strain_free(GwStrainSeries *s);

/** h_det(t) = F_plus * h_plus(t) + F_cross * h_cross(t) */
double gw_detector_strain(double F_plus, double F_cross,
                          double h_plus, double h_cross);

/**
 * Ground-based interferometer antenna pattern (90 deg arms).
 *
 * F_plus  = 1/2 (1+cos^2 theta) cos(2 phi) cos(2 psi)
 *           - cos theta sin(2 phi) sin(2 psi)
 * F_cross = 1/2 (1+cos^2 theta) cos(2 phi) sin(2 psi)
 *           + cos theta sin(2 phi) cos(2 psi)
 */
void gw_antenna_pattern(double *F_plus, double *F_cross,
                        double theta, double phi, double psi);

/**
 * Space-based detector (LISA-like) antenna pattern.
 * LISA has 60 deg opening angle instead of 90 deg.
 */
void gw_antenna_pattern_lisa(double *F_plus, double *F_cross,
                             double theta, double phi, double psi);

/* ================================================================
 * L2 — Polarization Tensors
 * ================================================================ */

void gw_pol_tensor_plus(GwTensor3 *e_plus);
void gw_pol_tensor_cross(GwTensor3 *e_cross);
void gw_pol_rotate(GwTensor3 *ep_out, GwTensor3 *ex_out, double psi);
void gw_strain_decompose(const GwTensor3 *h, double *h_plus, double *h_cross);

/* ================================================================
 * L4 — Fundamental Scales
 * ================================================================ */

/** PN expansion parameter x = (pi G M f / c^3)^{2/3} */
double gw_pn_parameter(double M_total, double f_gw);

/** Orbital separation: a = (G M / omega^2)^{1/3}, omega = pi f */
double gw_orbital_separation(double M_total, double f_gw);

/** ISCO frequency: f_ISCO ~ c^3 / (6^{3/2} pi G M) */
double gw_isco_frequency(double M_total);

/** Schwarzschild radius Rs = 2 G M / c^2 */
double gw_schwarzschild_radius(double M);

/** Light-ring frequency: f_LR ~ c^3 / (3^{3/2} 2 pi G M) */
double gw_lightring_frequency(double M_total);

#endif /* GW_CORE_H */
