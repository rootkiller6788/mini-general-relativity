/**
 * @file dg_gravwave.c
 * @brief Gravitational wave strain computation — L7 Application
 *
 * Knowledge: L7 — Gravitational wave polarizations in TT gauge
 *            Uses curvature tensor and geodesic deviation for detector response.
 *
 * Reference: Wald Sec.4.4, Carroll Sec.7, Schutz Sec.9
 * Application: LIGO, Virgo, KAGRA gravitational wave detectors
 *              (Nobel Prize 2017: Barish, Thorne, Weiss)
 *
 * Gravitational waves in the transverse-traceless (TT) gauge:
 *   h_{mu nu} = h_plus * e^+_{mu nu} + h_cross * e^x_{mu nu}
 *
 * where e^+ and e^x are the plus and cross polarization tensors.
 * The geodesic deviation equation gives the detector strain:
 *   d^2 xi^i / dt^2 = -R^i_{0j0} xi^j  (in proper detector frame)
 *
 * For a GW propagating in the z-direction:
 *   h_plus(t,z) = A_plus * cos(omega*(t - z/c))
 *   h_cross(t,z) = A_cross * cos(omega*(t - z/c) + delta_phi)
 */

#include "dg_curvature.h"
#include "dg_geodesic.h"
#include "dg_metric.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Compute the gravitational wave metric perturbation in TT gauge.
 *
 * The spacetime metric is g_{mu nu} = eta_{mu nu} + h_{mu nu}
 * where |h_{mu nu}| << 1 is the perturbation.
 *
 * In TT gauge for a wave traveling in the +z direction:
 *   h_{mu nu} = [[0, 0, 0, 0],
 *                [0,  h_plus, h_cross, 0],
 *                [0, h_cross, -h_plus, 0],
 *                [0, 0, 0, 0]]
 *
 * @param h_plus   Plus polarization amplitude
 * @param h_cross  Cross polarization amplitude
 * @param g_out    Output: 4x4 metric components (flat array, dim*dim)
 * @param dim      Dimension (4)
 *
 * Course Ref: Carroll Eq.(7.147) — TT Gauge Metric
 */
void gravwave_metric_tt(double h_plus, double h_cross, double *g_out, int dim) {
    if (!g_out || dim != 4) return;

    /* Start with Minkowski metric */
    for (int i = 0; i < 16; i++) g_out[i] = 0.0;
    g_out[0*4+0] = -1.0;
    g_out[1*4+1] =  1.0 + h_plus;      /* g_xx = 1 + h_+ */
    g_out[2*4+2] =  1.0 - h_plus;      /* g_yy = 1 - h_+ */
    g_out[3*4+3] =  1.0;
    g_out[1*4+2] =  h_cross;           /* g_xy = h_x */
    g_out[2*4+1] =  h_cross;           /* g_yx = h_x (symmetric) */
}

/**
 * @brief Compute the detector response to a gravitational wave.
 *
 * The strain measured by an interferometric detector (LIGO-like):
 *   h(t) = F_plus(theta, phi, psi) * h_plus(t)
 *        + F_cross(theta, phi, psi) * h_cross(t)
 *
 * where F_plus, F_cross are the antenna pattern functions.
 *
 * For a detector with arms along the x and y axes (simplified):
 *   F_plus = (1/2)(1 + cos^2(theta)) cos(2*phi) cos(2*psi)
 *          - cos(theta) sin(2*phi) sin(2*psi)
 *   F_cross = (1/2)(1 + cos^2(theta)) cos(2*phi) sin(2*psi)
 *           + cos(theta) sin(2*phi) cos(2*psi)
 *
 * @param theta   Polar angle of GW source (radians)
 * @param phi     Azimuthal angle of GW source (radians)
 * @param psi     Polarization angle (radians)
 * @param F_plus  Output: plus polarization antenna factor
 * @param F_cross Output: cross polarization antenna factor
 *
 * Real-world reference: LIGO Scientific Collaboration,
 *   "Observation of Gravitational Waves from a Binary Black Hole Merger"
 *   PRL 116, 061102 (2016) — GW150914
 */
void gravwave_antenna_pattern(double theta, double phi, double psi,
                               double *F_plus, double *F_cross) {
    double ct = cos(theta);
    double c2p = cos(2.0 * phi);
    double s2p = sin(2.0 * phi);
    double c2psi = cos(2.0 * psi);
    double s2psi = sin(2.0 * psi);

    *F_plus  = 0.5 * (1.0 + ct*ct) * c2p * c2psi - ct * s2p * s2psi;
    *F_cross = 0.5 * (1.0 + ct*ct) * c2p * s2psi + ct * s2p * c2psi;
}

/**
 * @brief Compute the Riemann tensor components for a GW in TT gauge.
 *
 * In linearized theory, the Riemann tensor for a plane GW in TT gauge is:
 *   R_{i0j0} = -(1/2) d^2 h_{ij}^{TT} / dt^2
 *
 * For a monochromatic wave h_{ij} = A_{ij} cos(omega*t):
 *   R_{i0j0} = (omega^2/2) * A_{ij} * cos(omega*t)
 *
 * This curvature drives the geodesic deviation that gravitational wave
 * detectors measure.
 *
 * @param h_plus   Plus polarization amplitude
 * @param h_cross  Cross polarization amplitude
 * @param omega    Angular frequency
 * @param t        Time
 * @param R_comp   Output: R_{i0j0} components for i,j=1,2,3 (9 values)
 *
 * Course Ref: Wald Sec.4.4a — Linearized Einstein Equations
 */
void gravwave_riemann_tt(double h_plus, double h_cross,
                          double omega, double t, double *R_comp) {
    if (!R_comp) return;
    /* d^2 h_{ij}/dt^2 = -omega^2 * A_{ij} * cos(omega*t)
     * R_{i0j0} = +(omega^2/2) * A_{ij} * cos(omega*t) */
    double amp = 0.5 * omega * omega * cos(omega * t);

    /* Initialize to zero */
    for (int i = 0; i < 9; i++) R_comp[i] = 0.0;

    /* R_{x0x0} = (omega^2/2) * h_plus * cos(omega*t) */
    R_comp[0*3+0] =  amp * h_plus;
    /* R_{x0y0} = (omega^2/2) * h_cross * cos(omega*t) */
    R_comp[0*3+1] =  amp * h_cross;
    /* R_{y0x0} = same (symmetric) */
    R_comp[1*3+0] =  amp * h_cross;
    /* R_{y0y0} = -(omega^2/2) * h_plus * cos(omega*t) */
    R_comp[1*3+1] = -amp * h_plus;
    /* R_{z0z0} = 0 (TT gauge, transverse) */
}

/**
 * @brief Compute the detector strain from geodesic deviation.
 *
 * For a detector with arm length L, the strain is:
 *   h(t) = Delta L / L = (1/2) (h_{ij} n^i n^j)
 * where n^i is the unit vector along the detector arm.
 *
 * Using the geodesic deviation equation:
 *   d^2 xi^i / dt^2 = -R^i_{0j0} xi^j
 *
 * For small perturbations: xi^i(t) = xi_0^i + (1/2) h^i_j xi_0^j
 *
 * @param R_comp  R_{i0j0} components (9 values, flat)
 * @param arm_dir Detector arm direction unit vector (3 components)
 * @param arm_len Detector arm length
 * @return        Strain h = Delta L / L
 *
 * LIGO specifications: 4km arms (LIGO Hanford, LIGO Livingston)
 * Sensitivity: h ~ 10^{-21} at ~100 Hz
 * GW150914: peak strain ~ 10^{-21}, frequency ~ 150 Hz
 */
double gravwave_detector_strain(const double *R_comp,
                                 const double *arm_dir, double arm_len) {
    if (!R_comp || !arm_dir) return 0.0;

    /* h = (1/2) R_{i0j0} * arm_len^2 * n^i * n^j / arm_len
     * Actually, for geodesic deviation with initial rest:
     *   xi^i(t) approx (1/2) R^i_{0j0} * xi_0^j * t^2
     * Strain = Delta L / L = xi / arm_len
     *
     * Simplified: h = (1/2) * R_{i0j0} * n^i * n^j * arm_len
     */
    double h = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            h += R_comp[i*3+j] * arm_dir[i] * arm_dir[j];
        }
    }
    h *= 0.5 * arm_len;
    return h;
}

/**
 * @brief Compute the characteristic strain of a binary black hole merger.
 *
 * For a circular binary with masses m1, m2, the GW frequency evolves as:
 *   f(t) = (1/pi) * sqrt(G*M / r(t)^3)
 *
 * The chirp mass: M_chirp = (m1*m2)^{3/5} / (m1+m2)^{1/5}
 * The strain amplitude: h ~ (4/D_L) * (G*M_chirp/c^2) * (pi*G*M_chirp*f/c^3)^{2/3}
 *
 * For GW150914: m1~36Msun, m2~29Msun, M_chirp~28Msun, D_L~410Mpc
 *   Peak strain ~ 1.0e-21 at f~150Hz
 *
 * @param m1       Primary mass (solar masses)
 * @param m2       Secondary mass (solar masses)
 * @param distance Luminosity distance (Mpc)
 * @param f_gw     Gravitational wave frequency (Hz)
 * @return         Strain amplitude h
 *
 * Reference: LIGO-Virgo Collaboration, PRL 116, 061102 (2016)
 */
double gravwave_chirp_strain(double m1, double m2,
                              double distance, double f_gw) {
    /* Physical constants */
    const double G_over_c2 = 1.4766e-27;   /* G/c^2 in m/kg */
    const double Msun = 1.989e30;          /* Solar mass in kg */
    const double Mpc = 3.0857e22;          /* Megaparsec in meters */

    /* Chirp mass in solar masses */
    double Mc = pow(m1 * m2, 0.6) / pow(m1 + m2, 0.2);

    /* Convert to geometric units */
    double Mc_geo = Mc * Msun * G_over_c2;  /* chirp mass in meters */
    double D_geo = distance * Mpc;            /* distance in meters */

    /* Strain amplitude (simplified formula) */
    double h = (4.0 / D_geo) * Mc_geo
             * pow(M_PI * Mc_geo * f_gw, 2.0/3.0);

    return h;
}

/**
 * @brief Compute GPS relativistic correction from Schwarzschild metric.
 *
 * GPS satellites orbit at r ~ 26,600 km from Earth's center.
 * Earth mass M ~ 5.97e24 kg, Schwarzschild radius 2GM/c^2 ~ 8.9 mm.
 *
 * Relativistic effects on GPS:
 *   1. Gravitational redshift (General Relativity): +45.7 us/day
 *      dtau/dt = sqrt(1 - 2GM/(rc^2)) ~ 1 - GM/(rc^2)
 *      Clock runs faster by GM/(rc^2)
 *   2. Time dilation from orbital velocity (Special Relativity): -7.2 us/day
 *      dtau/dt = sqrt(1 - v^2/c^2) ~ 1 - v^2/(2c^2)
 *      Clock runs slower by v^2/(2c^2)
 *
 * Net effect: +45.7 - 7.2 = +38.5 us/day
 * This corresponds to ~11.5 km/day position error if uncorrected.
 *
 * @param r_orbit  Orbital radius (meters)
 * @param M_earth  Earth mass (kg)
 * @return         Fractional frequency shift (Delta f / f)
 *
 * Application: Global Positioning System (GPS)
 *   Without GR correction, GPS would accumulate ~11 km error per day.
 *   GPS satellites broadcast at 10.23 MHz, corrected to 10.22999999543 MHz
 *   at launch to compensate for relativistic effects.
 *
 * Reference: Ashby, "Relativity in the Global Positioning System",
 *   Living Reviews in Relativity 6, 1 (2003)
 */
double gps_relativistic_correction(double r_orbit, double M_earth) {
    const double G = 6.67430e-11;   /* m^3 kg^{-1} s^{-2} */
    const double c = 2.99792458e8;  /* m/s */

    /* Gravitational blueshift factor */
    double grav_shift = G * M_earth / (r_orbit * c * c);

    /* Orbital velocity: v^2 = GM/r */
    double v2 = G * M_earth / r_orbit;
    double sr_shift = 0.5 * v2 / (c * c);

    /* Net fractional frequency shift */
    double net_shift = grav_shift - sr_shift;

    return net_shift;
}

/**
 * @brief Compute the Lense-Thirring precession (frame-dragging) rate.
 *
 * For a rotating mass (Kerr metric), a gyroscope orbiting at radius r
 * precesses due to frame-dragging:
 *   Omega_LT = (2*G*J)/(c^2 * r^3)
 *
 * where J is the angular momentum of the central body.
 *
 * Gravity Probe B measured this effect for Earth:
 *   Predicted: 39.2 milliarcseconds/year
 *   Measured: 37.2 +/- 7.2 milliarcseconds/year
 *
 * @param J        Angular momentum (kg m^2 / s)
 * @param r        Orbital radius (m)
 * @param M        Central mass (kg)
 * @return         Lense-Thirring precession rate (rad/s)
 *
 * Application: Gravity Probe B (NASA, 2004-2010)
 *   Measured both geodetic effect (6.6 arcsec/yr) and
 *   frame-dragging (0.039 arcsec/yr) confirming GR predictions.
 *
 * Reference: Everitt et al., "Gravity Probe B: Final Results
 *   of a Space Experiment to Test General Relativity",
 *   PRL 106, 221101 (2011)
 */
double lense_thirring_precession(double J, double r, double M) {
    const double G = 6.67430e-11;
    const double c = 2.99792458e8;

    /* Omega_LT = (2*G*J) / (c^2 * r^3) */
    double omega = (2.0 * G * J) / (c * c * r * r * r);

    return omega;
}