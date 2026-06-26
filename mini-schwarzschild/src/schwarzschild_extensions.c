/**
 * schwarzschild_extensions.c - Kerr, Reissner-Nordstrom, BH thermodynamics
 *
 * Reference: Wald (1984) Ch.7, Ch.12; Carroll (2004) Ch.6
 * MIT 8.962 - General Relativity
 * L8 - Advanced Topics: Kerr metric, RN metric, Penrose process, BH thermodynamics
 */

#include "schwarzschild_extensions.h"
#include "schwarzschild_defs.h"
#include <math.h>

/* ========================================================================
 * Reissner-Nordstrom Black Hole
 * ======================================================================== */

void rn_black_hole_init(double mass, double charge, ReissnerNordstromBH *rn) {
    if (!rn) return;
    rn->mass = mass;
    rn->charge = charge;
    rn->geometric_mass = geometric_mass(mass);
    double k = 8.987551787e9;
    double G = SCHW_G_N;
    double c = SCHW_C;
    double Q2_geo = G * charge * charge * k / (c * c * c * c);
    rn->geometric_charge = sqrt(Q2_geo);
    double m_geo = rn->geometric_mass;
    double Q_geo = rn->geometric_charge;
    if (Q_geo > m_geo) {
        rn->r_plus = 0.0;
        rn->r_minus = 0.0;
        rn->is_extremal = 0;
        rn->is_naked_singularity = 1;
    } else if (fabs(Q_geo - m_geo) < 1e-15 * m_geo) {
        rn->r_plus = m_geo;
        rn->r_minus = m_geo;
        rn->is_extremal = 1;
        rn->is_naked_singularity = 0;
    } else {
        rn->r_plus = m_geo + sqrt(m_geo * m_geo - Q_geo * Q_geo);
        rn->r_minus = m_geo - sqrt(m_geo * m_geo - Q_geo * Q_geo);
        rn->is_extremal = 0;
        rn->is_naked_singularity = 0;
    }
}

double rn_delta(double r, const ReissnerNordstromBH *rn) {
    if (!rn) return 0.0;
    double m = rn->geometric_mass;
    double Q = rn->geometric_charge;
    return r * r - 2.0 * m * r + Q * Q;
}

double rn_g_tt(double r, const ReissnerNordstromBH *rn) {
    if (!rn) return 0.0;
    double Delta = rn_delta(r, rn);
    return -Delta / (r * r);
}

double rn_g_rr(double r, const ReissnerNordstromBH *rn) {
    if (!rn) return 0.0;
    double Delta = rn_delta(r, rn);
    if (Delta <= 0.0) return INFINITY;
    return r * r / Delta;
}

double rn_horizon_outer(const ReissnerNordstromBH *rn) {
    if (!rn) return 0.0;
    return rn->r_plus;
}

double rn_horizon_inner(const ReissnerNordstromBH *rn) {
    if (!rn) return 0.0;
    return rn->r_minus;
}

double rn_extremal_charge(double mass) {
    double m_geo = geometric_mass(mass);
    double Q2_max = m_geo * m_geo;
    double k = 8.987551787e9;
    double G = SCHW_G_N;
    double c = SCHW_C;
    return sqrt(Q2_max * c * c * c * c / (G * k));
}

double rn_hawking_temperature(const ReissnerNordstromBH *rn) {
    if (!rn || rn->is_naked_singularity) return 0.0;
    double hbar = 1.054571817e-34;
    double kB = 1.380649e-23;
    double c = SCHW_C;
    double r_plus = rn->r_plus;
    double r_minus = rn->r_minus;
    if (fabs(r_plus - r_minus) < 1e-15) return 0.0;
    double kappa = (r_plus - r_minus) / (2.0 * r_plus * r_plus);
    return hbar * c * kappa / (2.0 * M_PI * kB);
}

double rn_electric_field(double r, const ReissnerNordstromBH *rn) {
    if (!rn || r <= 0.0) return 0.0;
    double Q = rn->charge;
    double k = 8.987551787e9;
    return k * Q / (r * r);
}

void rn_metric_covariant(const SchwarzschildPoint *p, const ReissnerNordstromBH *rn,
                           SchwarzschildSymmetric4x4 *metric) {
    if (!p || !rn || !metric) return;
    double r = p->r;
    double st = sin(p->theta);
    double Delta = rn_delta(r, rn);
    metric->g00 = -Delta / (r * r);
    metric->g01 = 0.0; metric->g02 = 0.0; metric->g03 = 0.0;
    metric->g11 = (Delta != 0.0) ? r * r / Delta : INFINITY;
    metric->g12 = 0.0; metric->g13 = 0.0;
    metric->g22 = r * r;
    metric->g23 = 0.0;
    metric->g33 = r * r * st * st;
}

double rn_kretschmann(double r, const ReissnerNordstromBH *rn) {
    if (!rn || r <= 0.0) return 0.0;
    double m = rn->geometric_mass;
    double Q = rn->geometric_charge;
    double r2 = r * r, r4 = r2 * r2, r6 = r4 * r2, r8 = r4 * r4;
    double m2 = m * m, Q2 = Q * Q;
    double K = (48.0 * m2 * r2 - 96.0 * m * Q2 * r + 56.0 * Q2 * Q2) / r8;
    return K;
}

/* ========================================================================
 * Kerr Black Hole
 * ======================================================================== */

void kerr_black_hole_init(double mass, double spin, KerrBH *kerr) {
    if (!kerr) return;
    kerr->mass = mass;
    kerr->spin = spin;
    kerr->geometric_mass = geometric_mass(mass);
    kerr->geometric_spin = spin / (SCHW_C * mass);
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    double a2 = a * a;
    double disc = m * m - a2;
    if (disc < 0.0) {
        kerr->r_plus = 0.0;
        kerr->r_minus = 0.0;
    } else {
        kerr->r_plus = m + sqrt(disc);
        kerr->r_minus = m - sqrt(disc);
    }
    kerr->r_ergo_equator = 2.0 * m;
    double Z1 = 1.0 + cbrt(1.0 - a2/(m*m)) * (cbrt(1.0 + a/m) + cbrt(1.0 - a/m));
    double Z2 = sqrt(3.0 * a2/(m*m) + Z1*Z1);
    kerr->isco_prograde = m * (3.0 + Z2 - sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2)));
    kerr->isco_retrograde = m * (3.0 + Z2 + sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2)));
}

double kerr_delta(double r, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    return r * r - 2.0 * m * r + a * a;
}

double kerr_sigma(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double a = kerr->geometric_spin;
    double ct = cos(theta);
    return r * r + a * a * ct * ct;
}

double kerr_g_tt(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double Sigma = kerr_sigma(r, theta, kerr);
    return -(1.0 - 2.0 * m * r / Sigma);
}

double kerr_g_tphi(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    double Sigma = kerr_sigma(r, theta, kerr);
    double st = sin(theta);
    return -2.0 * m * a * r * st * st / Sigma;
}

double kerr_g_rr(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double Delta = kerr_delta(r, kerr);
    double Sigma = kerr_sigma(r, theta, kerr);
    if (Delta <= 0.0) return INFINITY;
    return Sigma / Delta;
}

double kerr_g_theta_theta(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    return kerr_sigma(r, theta, kerr);
}

double kerr_g_phi_phi(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    double Sigma = kerr_sigma(r, theta, kerr);
    double st = sin(theta);
    double term1 = r * r + a * a;
    double term2 = 2.0 * m * a * a * r * st * st / Sigma;
    return (term1 + term2) * st * st;
}

double kerr_horizon_outer(const KerrBH *kerr) {
    return kerr ? kerr->r_plus : 0.0;
}

double kerr_horizon_inner(const KerrBH *kerr) {
    return kerr ? kerr->r_minus : 0.0;
}

double kerr_ergosphere_radius(double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    double ct = cos(theta);
    double disc = m * m - a * a * ct * ct;
    if (disc < 0.0) return m;
    return m + sqrt(disc);
}

double kerr_isco_prograde(const KerrBH *kerr) {
    return kerr ? kerr->isco_prograde : 0.0;
}

double kerr_isco_retrograde(const KerrBH *kerr) {
    return kerr ? kerr->isco_retrograde : 0.0;
}

double kerr_frame_dragging_omega(double r, double theta, const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    double Sigma = kerr_sigma(r, theta, kerr);
    double term = r * r + a * a;
    double denom = term * Sigma + 2.0 * m * a * a * r * sin(theta) * sin(theta);
    if (fabs(denom) < 1e-30) return 0.0;
    return 2.0 * m * a * r / denom;
}

double kerr_photon_sphere_prograde(const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    return 2.0 * m * (1.0 + cos(2.0/3.0 * acos(-a/m)));
}

double kerr_photon_sphere_retrograde(const KerrBH *kerr) {
    if (!kerr) return 0.0;
    double m = kerr->geometric_mass;
    double a = kerr->geometric_spin;
    return 2.0 * m * (1.0 + cos(2.0/3.0 * acos(a/m)));
}

void kerr_metric_covariant(const SchwarzschildPoint *p, const KerrBH *kerr,
                             SchwarzschildSymmetric4x4 *metric) {
    if (!p || !kerr || !metric) return;
    double r = p->r, theta = p->theta;
    metric->g00 = kerr_g_tt(r, theta, kerr);
    metric->g01 = 0.0;
    metric->g02 = 0.0;
    metric->g03 = kerr_g_tphi(r, theta, kerr);
    metric->g11 = kerr_g_rr(r, theta, kerr);
    metric->g12 = 0.0;
    metric->g13 = 0.0;
    metric->g22 = kerr_g_theta_theta(r, theta, kerr);
    metric->g23 = 0.0;
    metric->g33 = kerr_g_phi_phi(r, theta, kerr);
}

void kerr_metric_contravariant(const SchwarzschildPoint *p, const KerrBH *kerr,
                                 SchwarzschildSymmetric4x4 *inv_metric) {
    if (!p || !kerr || !inv_metric) return;
    double r = p->r, theta = p->theta;
    double g_tt = kerr_g_tt(r, theta, kerr);
    double g_tphi = kerr_g_tphi(r, theta, kerr);
    double g_rr = kerr_g_rr(r, theta, kerr);
    double g_thth = kerr_g_theta_theta(r, theta, kerr);
    double g_phiphi = kerr_g_phi_phi(r, theta, kerr);
    double det_tphi = g_tt * g_phiphi - g_tphi * g_tphi;
    if (fabs(det_tphi) < 1e-30) return;
    inv_metric->g00 = g_phiphi / det_tphi;
    inv_metric->g01 = 0.0; inv_metric->g02 = 0.0;
    inv_metric->g03 = -g_tphi / det_tphi;
    inv_metric->g11 = (g_rr != 0.0) ? 1.0 / g_rr : INFINITY;
    inv_metric->g12 = 0.0; inv_metric->g13 = 0.0;
    inv_metric->g22 = (g_thth != 0.0) ? 1.0 / g_thth : INFINITY;
    inv_metric->g23 = 0.0;
    inv_metric->g33 = g_tt / det_tphi;
}

/* ========================================================================
 * Black Hole Thermodynamics
 * ======================================================================== */

double schwarzschild_bekenstein_hawking_entropy(double mass) {
    if (mass <= 0.0) return 0.0;
    double A = schwarzschild_horizon_area(mass);
    double hbar = 1.054571817e-34;
    double kB = 1.380649e-23;
    double c = SCHW_C;
    double G = SCHW_G_N;
    return kB * A * c * c * c / (4.0 * G * hbar);
}

double schwarzschild_hawking_temperature(double mass) {
    if (mass <= 0.0) return 0.0;
    double hbar = 1.054571817e-34;
    double kB = 1.380649e-23;
    double c = SCHW_C;
    double G = SCHW_G_N;
    return hbar * c * c * c / (8.0 * M_PI * G * mass * kB);
}

double schwarzschild_hawking_luminosity(double mass) {
    if (mass <= 0.0) return 0.0;
    double c = SCHW_C;
    double G = SCHW_G_N;
    double hbar = 1.054571817e-34;
    return hbar * c * c * c * c * c * c / (15360.0 * M_PI * G * G * mass * mass);
}

double schwarzschild_evaporation_time(double mass) {
    if (mass <= 0.0) return 0.0;
    double c = SCHW_C;
    double G = SCHW_G_N;
    double hbar = 1.054571817e-34;
    return 5120.0 * M_PI * G * G * mass * mass * mass / (hbar * c * c * c * c);
}

double schwarzschild_first_law_dM_dA(double mass) {
    if (mass <= 0.0) return 0.0;
    double T = schwarzschild_hawking_temperature(mass);
    double kB = 1.380649e-23;
    double c = SCHW_C;
    double G = SCHW_G_N;
    double kappa = schwarzschild_surface_gravity(mass);
    double dA_dM = 32.0 * M_PI * G * G * mass / (c * c * c * c);
    return dA_dM;
}

double schwarzschild_area_increase(double mass_initial, double mass_added) {
    double A_i = schwarzschild_horizon_area(mass_initial);
    double A_f = schwarzschild_horizon_area(mass_initial + mass_added);
    return A_f - A_i;
}

/* ========================================================================
 * Penrose Process and Energy Extraction
 * ======================================================================== */

double kerr_penrose_max_efficiency(double spin_param) {
    double m = spin_param;
    double a = fmin(fabs(spin_param), m * 0.999999);
    if (fabs(a) < 1e-15) return 0.0;
    double a_over_m = a / m;
    double max_eff = 1.0 - 1.0 / sqrt(2.0);
    double factor = a_over_m;
    return max_eff * factor * factor * (2.0 - factor);
}

double kerr_blandford_znajek_power(double mass, double spin,
                                     double magnetic_field) {
    if (mass <= 0.0) return 0.0;
    double m_geo = geometric_mass(mass);
    double a = spin / (SCHW_C * mass);
    double a_over_m = a / m_geo;
    double r_h = m_geo + sqrt(m_geo * m_geo - a * a);
    double c = SCHW_C;
    double mu0 = 4.0 * M_PI * 1e-7;
    double B_sq_over_mu0 = magnetic_field * magnetic_field / mu0;
    return (1.0 / 8.0) * B_sq_over_mu0 * r_h * r_h * a_over_m * a_over_m * c;
}

/* ========================================================================
 * Quasinormal Modes
 * ======================================================================== */

void schwarzschild_qnm_fundamental_l2_n0(double mass,
                                            double *freq_Hz,
                                            double *damping_time_s) {
    if (!freq_Hz || !damping_time_s || mass <= 0.0) return;
    double m_geo = geometric_mass(mass);
    double omega_real = 0.37367 / m_geo;
    double omega_imag = 0.08896 / m_geo;
    *freq_Hz = omega_real * SCHW_C / (2.0 * M_PI);
    *damping_time_s = 1.0 / (omega_imag * SCHW_C);
}

double kerr_qnm_ringdown_frequency(double mass, double spin) {
    if (mass <= 0.0) return 0.0;
    double m_geo = geometric_mass(mass);
    double a = fabs(spin) / (SCHW_C * mass);
    double a_over_m = fmin(a / m_geo, 0.9999);
    double f_lmn = 1.0 / (2.0 * M_PI * m_geo)
                   * (1.0 - 0.63 * pow(1.0 - a_over_m, 0.3));
    return f_lmn * SCHW_C;
}

double schwarzschild_qnm_quality_factor(double mass) {
    if (mass <= 0.0) return 0.0;
    double m_geo = geometric_mass(mass);
    double omega_real = 0.37367 / m_geo;
    double omega_imag = 0.08896 / m_geo;
    double tau = 1.0 / omega_imag;
    double f = omega_real / (2.0 * M_PI);
    return M_PI * f * tau;
}

double schwarzschild_ringdown_strain(double mass, double distance,
                                       double radiated_fraction) {
    if (mass <= 0.0 || distance <= 0.0) return 0.0;
    double omega_real = 0.37367 / geometric_mass(mass);
    double tau = 1.0 / (0.08896 / geometric_mass(mass));
    double E_rad = radiated_fraction * mass * SCHW_C * SCHW_C;
    double h0 = sqrt(E_rad * SCHW_G_N / (SCHW_C * SCHW_C * SCHW_C * SCHW_C * distance * distance));
    return h0;
}
