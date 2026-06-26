/**
 * @file    gw_waveforms.c
 * @brief   GW Waveform Generation — Inspiral, IMR, Ringdown, Higher Modes
 *
 * L5 — Time-domain inspiral (TaylorT4), frequency-domain SPA
 * L6 — IMR phenomenological model, ringdown (QNM), higher harmonics
 * L7 — Amplitude scaling, characteristic strain
 */

#include "gw_waveforms.h"
#include "gw_binary.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * L5 — Time-Domain Inspiral Waveform via TaylorT4
 *
 * Integrates the TaylorT4 ODE for f(t) and Phi(t), then
 * computes h_{+,x}(t) using the Newtonian quadrupole amplitude.
 * ================================================================ */

int gw_waveform_inspiral_td(GwStrainSeries *s,
    const GwBinaryParams *params,
    double f_start, double f_end,
    double dt, size_t max_samples) {
    if (!s || max_samples == 0) return -1;
    if (f_start >= f_end || dt <= 0.0) return -1;

    size_t n = 0;
    double f = f_start;
    double Phi = params->phi_ref;
    double t  = 0.0;

    double Mt = params->M_total;
    double eta = params->eta;
    double mu  = params->m1 * params->m2 / Mt;
    double D_L = params->D_L;
    double iota = params->iota;
    double Mc  = params->chirp_mass;

    while (f < f_end && n < max_samples) {
        /* Store current state */
        if (n < s->n) {
            s->t[n] = t;

            double a = gw_orbital_separation(Mt, f);
            double Phi_t = Phi;

            /* h(t) from Newtonian quadrupole */
            double fac = 2.0 * GW_G * GW_G * Mt * mu
                       / (GW_C*GW_C*GW_C*GW_C * D_L * a);

            s->hp[n] = -fac * (1.0 + cos(iota)*cos(iota)) * cos(2.0 * Phi_t);
            s->hx[n] = -fac * 2.0 * cos(iota) * sin(2.0 * Phi_t);
        }

        /* Advance via TaylorT4 */
        double df_dt, dphi_dt;
        gw_taylor_t4_rhs(f, Mt, eta, &df_dt, &dphi_dt);

        f   += df_dt * dt;
        Phi += dphi_dt * dt;
        t   += dt;
        n++;
    }

    return (int)n;
}

/* ================================================================
 * L5 — Frequency-Domain Waveform (Stationary Phase Approximation)
 *
 * SPA: h(f) = A(f) e^{i Psi(f)}
 *
 * A(f) = (1/D_L) sqrt(5/24) pi^{-2/3} (G Mc / c^3)^{5/6} f^{-7/6}
 *
 * Psi(f) = TaylorF2 phase
 * ================================================================ */

void gw_waveform_fd_spa(double *h_re, double *h_im,
    const double *f_arr, int nf,
    const GwBinaryParams *params,
    double t_c, double phi_c) {
    double Mc   = params->chirp_mass;
    double D_L  = params->D_L;
    double M_total = params->M_total;
    double eta  = params->eta;
    double iota = params->iota;

    double Mc_si = Mc * GW_G / (GW_C * GW_C * GW_C);  /* geometric chirp mass [s] */

    for (int i = 0; i < nf; i++) {
        double f = f_arr[i];
        if (f <= 0.0) {
            h_re[i] = 0.0;
            h_im[i] = 0.0;
            continue;
        }

        /* Amplitude: A(f) ~ (1/D_L) * Mc^{5/6} * f^{-7/6} */
        double amp = (GW_C / D_L)
                     * sqrt(5.0 / 24.0) * pow(GW_PI, -2.0/3.0)
                     * pow(Mc_si, 5.0/6.0) * pow(f, -7.0/6.0);
        amp *= sqrt(1.0 + 6.0*cos(iota)*cos(iota) + cos(iota)*cos(iota)*cos(iota)*cos(iota));

        /* Phase: TaylorF2 */
        double psi_f = gw_taylor_f2_phase(f, M_total, eta, t_c, phi_c);

        h_re[i] = amp * cos(psi_f);
        h_im[i] = amp * sin(psi_f);
    }
}

/* ================================================================
 * L6 — IMR Phenomenological Waveform (simplified IMRPhenomD-like)
 *
 * Inspiral: SPA with TaylorF2 phase
 * Merger-Ringdown: Lorentzian amplitude + frequency-dependent phase
 * ================================================================ */

void gw_waveform_imr_fd(double *h_re, double *h_im,
    const double *f_arr, int nf,
    const GwBinaryParams *params,
    double t_c, double phi_c) {
    double Mc   = params->chirp_mass;
    double D_L  = params->D_L;
    double M_total = params->M_total;
    double eta  = params->eta;
    double iota = params->iota;

    double f_isco = gw_isco_frequency(M_total);
    double Mc_si = Mc * GW_G / (GW_C * GW_C * GW_C);

    for (int i = 0; i < nf; i++) {
        double f = f_arr[i];
        if (f <= 0.0) {
            h_re[i] = 0.0; h_im[i] = 0.0;
            continue;
        }

        /* Amplitude */
        double amp;
        double psi_f = gw_taylor_f2_phase(f, M_total, eta, t_c, phi_c);

        if (f < f_isco) {
            /* Inspiral: SPA amplitude */
            amp = (GW_C / D_L)
                  * sqrt(5.0 / 24.0) * pow(GW_PI, -2.0/3.0)
                  * pow(Mc_si, 5.0/6.0) * pow(f, -7.0/6.0);
        } else {
            /* Merger-ringdown: Lorentzian */
            double f_rd, tau;
            double M_f, a_f;
            gw_final_mass_spin(&M_f, &a_f, params->m1, params->m2);
            gw_qnm_220_mode(&f_rd, &tau, M_f, a_f);

            double amp_isco = (GW_C / D_L)
                              * sqrt(5.0 / 24.0) * pow(GW_PI, -2.0/3.0)
                              * pow(Mc_si, 5.0/6.0) * pow(f_isco, -7.0/6.0);

            double delta_f = f - f_rd;
            double denom = 1.0 + 4.0 * tau*tau * delta_f*delta_f * 4.0*GW_PI*GW_PI;
            amp = amp_isco * pow(f_isco / f_rd, 2.0) / sqrt(denom);

            /* Phase correction for merger-ringdown */
            psi_f += 2.0 * GW_PI * f * (t_c + 0.01); /* crude time shift */
        }

        amp *= sqrt(1.0 + 6.0*cos(iota)*cos(iota)
                    + cos(iota)*cos(iota)*cos(iota)*cos(iota));

        h_re[i] = amp * cos(psi_f);
        h_im[i] = amp * sin(psi_f);
    }
}

/* ================================================================
 * L6 — Ringdown Waveform (Damped Sinusoid)
 * ================================================================ */

void gw_waveform_ringdown(GwStrainSeries *s,
    double M_f, double a_f, double A_rd, double phi_rd,
    double dt, size_t n) {
    double f_rd, tau;
    gw_qnm_220_mode(&f_rd, &tau, M_f, a_f);

    double omega = 2.0 * GW_PI * f_rd;

    for (size_t i = 0; i < n; i++) {
        double t = i * dt;
        s->t[i]  = t;
        double env = A_rd * exp(-t / tau);
        s->hp[i] = env * cos(omega * t + phi_rd);
        s->hx[i] = env * cos(omega * t + phi_rd) * 0.5; /* typical ratio */
    }
}

/* ================================================================
 * L6 — Kerr QNM Frequencies (l=m=2, n=0)
 *
 * Berti, Cardoso & Will (2006) fitting formula:
 *
 *   M omega_220 = f_22(a) + i * g_22(a)
 *
 *   f_22(a) = 1.5251 - 1.1568 (1-a)^{0.1292}
 *   Q_22(a) = 0.7000 + 1.4187 (1-a)^{-0.4990}
 *
 * Then: omega_R = f_22 / M, omega_I = -omega_R / (2 Q_22)
 *       tau = 1 / |omega_I|
 * ================================================================ */

void gw_qnm_220_mode(double *f_rd, double *tau, double M_f, double a_f) {
    if (a_f > 0.999) a_f = 0.999;
    if (a_f < 0.0)   a_f = 0.0;

    double one_minus_a = 1.0 - a_f;
    if (one_minus_a < 1e-6) one_minus_a = 1e-6;

    double f22 = 1.5251 - 1.1568 * pow(one_minus_a, 0.1292);
    double Q22 = 0.7000 + 1.4187 * pow(one_minus_a, -0.4990);

    /* Geometric mass: M_geom = G M / c^2 [m]
     * Angular frequency: omega = f22 * c^3 / (G M) [rad/s]
     * Frequency: f = omega / (2 pi) */
    double omega_R = f22 * GW_C * GW_C * GW_C / (GW_G * M_f);
    *f_rd = omega_R / (2.0 * GW_PI);

    /* Damping time: tau = 2 Q / omega_R */
    *tau = 2.0 * Q22 / omega_R;
}

/* ================================================================
 * L6 — Final Mass and Spin (Barausse et al. 2012 fitting)
 *
 * Simplified fitting formula for aligned-spin BBH mergers.
 * More accurate: use NRSurrogate or EOB.
 * ================================================================ */

void gw_final_mass_spin(double *M_f, double *a_f, double m1, double m2) {
    double M = m1 + m2;
    double eta = m1 * m2 / (M * M);

    /* Radiated energy fraction: E_rad / M ~ 0.05 for equal-mass
     * Fitting: E_rad/M ~ 4 eta^2 (leading order) + corrections */
    double E_rad_frac = 4.0 * eta * eta * (1.0 - 2.0*eta);

    *M_f = M * (1.0 - E_rad_frac);

    /* Final spin: a_f ~ 2 sqrt(3) eta - 3.87 eta^2 + ... */
    *a_f = 2.0 * sqrt(3.0) * eta - 3.87 * eta * eta;
    if (*a_f > 0.998) *a_f = 0.998;
    if (*a_f < 0.0)   *a_f = 0.0;
}

/* ================================================================
 * L6 — Higher Mode Amplitudes (Leading PN Order)
 *
 * For non-precessing binaries, the higher harmonics (l,m) relative
 * to the dominant (2,2) mode scale with eta and the mass ratio.
 *
 * h_{22} dominant
 * h_{33} ~ eta * (amplitude factor)
 * h_{44} ~ eta^2 * (amplitude factor)
 * h_{21} ~ delta_m * (amplitude factor) where delta_m = (m1-m2)/M
 * ================================================================ */

int gw_higher_mode_amplitude(double *amp_lm, int l, int m, double eta) {
    /* Relative amplitude w.r.t (2,2) at leading PN order.
     * Values from Blanchet et al. (2008). */
    double delta_m = sqrt(1.0 - 4.0*eta);  /* = (m1-m2)/M for m1>=m2 */

    if (l == 2 && m == 2) {
        *amp_lm = 1.0; return 0;
    } else if (l == 3 && m == 3) {
        *amp_lm = -0.75 * delta_m * sqrt(3.0/5.0); return 0;
    } else if (l == 4 && m == 4) {
        *amp_lm = (4.0/9.0) * (1.0 - 3.0*eta) * sqrt(2.0/7.0); return 0;
    } else if (l == 2 && m == 1) {
        *amp_lm = 0.0; /* vanishes for non-precessing at leading order */ return 0;
    } else if (l == 3 && m == 2) {
        *amp_lm = (1.0/3.0) * delta_m; return 0;
    } else {
        *amp_lm = 0.0;
        return -1;  /* mode not available at this order */
    }
}

/* ================================================================
 * L6 — Waveform with Higher Modes
 * ================================================================ */

void gw_waveform_higher_modes(double *h_re, double *h_im,
    const double *f_arr, int nf,
    const GwBinaryParams *params, double iota,
    double t_c, double phi_c) {
    double Mc   = params->chirp_mass;
    double D_L  = params->D_L;
    double M_total = params->M_total;
    double eta  = params->eta;
    double Mc_si = Mc * GW_G / (GW_C * GW_C * GW_C);

    for (int i = 0; i < nf; i++) {
        h_re[i] = 0.0; h_im[i] = 0.0;
    }

    /* Dominant (2,2) mode: phase = 2*Psi_{22} */
    for (int i = 0; i < nf; i++) {
        double f = f_arr[i];
        if (f <= 0.0) continue;

        double amp = (GW_C / D_L)
                     * sqrt(5.0 / 24.0) * pow(GW_PI, -2.0/3.0)
                     * pow(Mc_si, 5.0/6.0) * pow(f, -7.0/6.0);

        double Psi22 = gw_taylor_f2_phase(f, M_total, eta, t_c, phi_c);

        /* Include (2,2), (3,3), (4,4), (2,1) contributions */
        /* (l,m) = (2,2): dominant */
        double Y22 = sqrt(5.0/(4.0*GW_PI)) * 0.5 * (1.0 + cos(iota)*cos(iota)); /* approximate */
        h_re[i] += amp * Y22 * cos(2.0*Psi22);
        h_im[i] += amp * Y22 * sin(2.0*Psi22);

        /* (3,3): m=3 harmonic at frequency (3/2)f (approximate mapping) */
        double f33 = 1.5 * f;
        if (f33 > 0) {
            double amp33_ratio;
            gw_higher_mode_amplitude(&amp33_ratio, 3, 3, eta);
            double amp33 = amp * fabs(amp33_ratio)
                           * pow(f33 / f, -7.0/6.0) * 0.5;
            double Psi33 = gw_taylor_f2_phase(f33, M_total, eta, t_c, phi_c);
            h_re[i] += amp33 * cos(3.0*Psi33);
            h_im[i] += amp33 * sin(3.0*Psi33);
        }
    }
}

/* ================================================================
 * L7 — Strain Amplitude & Characteristic Strain
 * ================================================================ */

double gw_strain_amplitude_fd(double Mc, double D_L, double f) {
    if (f <= 0.0 || D_L <= 0.0) return 0.0;
    double Mc_si = Mc * GW_G / (GW_C * GW_C * GW_C);
    return (GW_C / D_L)
           * sqrt(5.0 / 24.0) * pow(GW_PI, -2.0/3.0)
           * pow(Mc_si, 5.0/6.0) * pow(f, -7.0/6.0);
}

double gw_characteristic_strain(double Mc, double D_L, double f) {
    double amp = gw_strain_amplitude_fd(Mc, D_L, f);
    return 2.0 * f * amp;
}
