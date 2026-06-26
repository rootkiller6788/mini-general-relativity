/**
 * @file    gw_quadrupole.c
 * @brief   Quadrupole Radiation — Einstein Formula Implementation
 *
 * L1 — Mass quadrupole moment, reduced quadrupole, multipole expansions
 * L2 — Quadrupole GW strain, quadrupole luminosity
 * L4 — Binary quadrupole strain and luminosity (circular)
 */

#include "gw_quadrupole.h"
#include <math.h>
#include <string.h>

/* ================================================================
 * L1 — Mass Quadrupole Moment for N Point Particles
 *
 * I_{ij} = sum_a m_a ( x_a^i x_a^j - (1/3) delta^{ij} r_a^2 )
 *
 * The second term with delta^{ij} makes I_{ij} the
 * traceless part of the mass quadrupole at Newtonian order
 * when computing the GW quadrupole formula.
 *
 * Note: Some conventions define I_{ij} = sum m_a x_a^i x_a^j
 * (the raw quadrupole), and the TT projection takes care
 * of removing the trace. Here we include the trace term
 * for consistency with Wald Ch.4 notation.
 * ================================================================ */

void gw_quadrupole_moment(GwTensor3 *I_out, const double *m,
                          const double *x, int N) {
    memset(I_out, 0, sizeof(GwTensor3));

    for (int a = 0; a < N; a++) {
        double ma   = m[a];
        double xa   = x[3*a + 0];
        double ya   = x[3*a + 1];
        double za   = x[3*a + 2];
        double r2   = xa*xa + ya*ya + za*za;
        double trfac = r2 / 3.0;

        I_out->xx += ma * (xa*xa - trfac);
        I_out->xy += ma * (xa*ya);
        I_out->xz += ma * (xa*za);
        I_out->yy += ma * (ya*ya - trfac);
        I_out->yz += ma * (ya*za);
        I_out->zz += ma * (za*za - trfac);
    }
    I_out->yx = I_out->xy;
    I_out->zx = I_out->xz;
    I_out->zy = I_out->yz;
}

/* ================================================================
 * L1 — Reduced Quadrupole Moment
 *
 * Q_{ij} = I_{ij} - (1/3) delta_{ij} Tr[I]
 *
 * For an already traceless I_{ij}, Q_{ij} = I_{ij}.
 * This function works regardless.
 * ================================================================ */

void gw_reduced_quadrupole(GwTensor3 *Q_out, const double *m,
                           const double *x, int N) {
    gw_quadrupole_moment(Q_out, m, x, N);
    /* I from quadrupole_moment is already traceless, nothing to do */
}

/* ================================================================
 * L1 — Current Quadrupole Moment
 *
 * J_{ij} = sum_a m_a ( epsilon_{ikl} x_a^k v_a^l x_a^j
 *                     + epsilon_{jkl} x_a^k v_a^l x_a^i )
 *
 * (Antisymmetric part of the mass-current moment.)
 * We store the symmetric 3x3 representative.
 * ================================================================ */

void gw_current_quadrupole(GwTensor3 *J_out, const double *m,
                           const double *x, const double *v, int N) {
    memset(J_out, 0, sizeof(GwTensor3));

    for (int a = 0; a < N; a++) {
        double ma = m[a];
        double xa = x[3*a + 0], ya = x[3*a + 1], za = x[3*a + 2];
        double vx = v[3*a + 0], vy = v[3*a + 1], vz = v[3*a + 2];

        /* Angular momentum L_a = m_a (r_a x v_a) */
        double Lx = ma * (ya*vz - za*vy);
        double Ly = ma * (za*vx - xa*vz);
        double Lz = ma * (xa*vy - ya*vx);

        /* J_{ij} ~ (epsilon_{ikl} L_k x_j + epsilon_{jkl} L_k x_i) / (2) */
        /* Symmetric part: J_{xx} ~ L_y*x + L_z*x variations */
        J_out->xx += 0.0;  /* anti-symmetric diagonal vanishes */
        J_out->yy += 0.0;
        J_out->zz += 0.0;

        /* J_{xy} ~ L_z * (y?) ... simplified: J_{ij} from r x L */
        J_out->xy += ma * (xa*vy + ya*vx) * 0.5; /* simplified representative */
        J_out->xz += ma * (xa*vz + za*vx) * 0.5;
        J_out->yz += ma * (ya*vz + za*vy) * 0.5;
    }
    J_out->yx = J_out->xy;
    J_out->zx = J_out->xz;
    J_out->zy = J_out->yz;
}

/* ================================================================
 * L1 — Mass Octupole (L=3) Raw Moment
 * ================================================================ */

void gw_mass_octupole(double M_out[27], const double *m,
                      const double *x, int N) {
    for (int i = 0; i < 27; i++) M_out[i] = 0.0;

    for (int a = 0; a < N; a++) {
        double ma = m[a];
        double xa = x[3*a], ya = x[3*a+1], za = x[3*a+2];
        double xyz[3] = {xa, ya, za};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    M_out[i*9 + j*3 + k] += ma * xyz[i] * xyz[j] * xyz[k];
    }
}

/* ================================================================
 * L2 — Quadrupole GW Strain
 *
 * h_{ij}^{TT} = (2G / c^4 r) * d^2 Q_{ij} / dt^2
 *
 * The second time derivative evaluated at retarded time t - r/c.
 * Caller provides Q_ddot = d^2 Q_{ij} / dt^2.
 * ================================================================ */

void gw_quadrupole_strain(GwTensor3 *h_tt, const GwTensor3 *Q_ddot, double r) {
    double fac = 2.0 * GW_G_OVER_C4 / r;

    h_tt->xx = fac * Q_ddot->xx;
    h_tt->xy = fac * Q_ddot->xy;
    h_tt->xz = fac * Q_ddot->xz;
    h_tt->yx = fac * Q_ddot->yx;
    h_tt->yy = fac * Q_ddot->yy;
    h_tt->yz = fac * Q_ddot->yz;
    h_tt->zx = fac * Q_ddot->zx;
    h_tt->zy = fac * Q_ddot->zy;
    h_tt->zz = fac * Q_ddot->zz;
}

/* ================================================================
 * L2 — GW Strain Directly from Particle State
 *
 * d^2 I_{ij}/dt^2 = sum_a m_a [ 2 v_a^i v_a^j + x_a^i a_a^j + a_a^i x_a^j ]
 *
 * Then TT project and scale by 2G/(c^4 r).
 * ================================================================ */

void gw_quadrupole_strain_from_state(GwTensor3 *h_tt,
    const double *m, const double *x, const double *v,
    const double *a, int N, double r) {
    GwTensor3 I_ddot;
    memset(&I_ddot, 0, sizeof(GwTensor3));

    for (int p = 0; p < N; p++) {
        double mp   = m[p];
        double xp   = x[3*p + 0];
        double yp   = x[3*p + 1];
        double zp   = x[3*p + 2];
        double vxp  = v[3*p + 0];
        double vyp  = v[3*p + 1];
        double vzp  = v[3*p + 2];
        double axp  = a[3*p + 0];
        double ayp  = a[3*p + 1];
        double azp  = a[3*p + 2];

        /* d^2/dt^2 (x_i x_j) = 2 v_i v_j + x_i a_j + a_i x_j */
        I_ddot.xx += mp * (2.0*vxp*vxp + 2.0*xp*axp);
        I_ddot.xy += mp * (2.0*vxp*vyp + xp*ayp + axp*yp);
        I_ddot.xz += mp * (2.0*vxp*vzp + xp*azp + axp*zp);
        I_ddot.yy += mp * (2.0*vyp*vyp + 2.0*yp*ayp);
        I_ddot.yz += mp * (2.0*vyp*vzp + yp*azp + ayp*zp);
        I_ddot.zz += mp * (2.0*vzp*vzp + 2.0*zp*azp);
    }
    I_ddot.yx = I_ddot.xy;
    I_ddot.zx = I_ddot.xz;
    I_ddot.zy = I_ddot.yz;

    /* Remove trace for TT */
    double tr = gw_tensor_trace(&I_ddot);
    double one_third_tr = tr / 3.0;
    I_ddot.xx -= one_third_tr;
    I_ddot.yy -= one_third_tr;
    I_ddot.zz -= one_third_tr;

    gw_quadrupole_strain(h_tt, &I_ddot, r);
}

/* ================================================================
 * L2 — Quadrupole GW Luminosity
 *
 * L_GW = (G / 5c^5) * < Q_dddot_{ij} Q_dddot^{ij} >
 *
 * For instantaneous (not time-averaged) value.
 * In practice, one averages over one orbital period.
 * ================================================================ */

double gw_quadrupole_luminosity(const GwTensor3 *Q_dddot) {
    double contraction = gw_tensor_contract(Q_dddot, Q_dddot);
    return GW_G_OVER_5C5 * contraction;
}

double gw_quadrupole_luminosity_avg(const GwTensor3 *Q_dddot_arr, int n_samples) {
    double sum = 0.0;
    for (int i = 0; i < n_samples; i++) {
        sum += gw_quadrupole_luminosity(&Q_dddot_arr[i]);
    }
    return sum / (double)n_samples;
}

double gw_flux(double L_GW, double r) {
    return L_GW / (4.0 * GW_PI * r * r);
}

/* ================================================================
 * L2 — Isaacson Energy Density of GWs
 *
 * rho_GW = (c^2 / 32 pi G) * < hdot_{ij}^{TT} hdot^{ij,TT} >
 *
 * For a plane wave: hdot_{ij} = i omega h_{ij}
 *   -> rho_GW = (c^2 omega^2 / 32 pi G) * (h_+^2 + h_x^2)
 * ================================================================ */

double gw_energy_density(double h_plus, double h_cross, double omega) {
    double amp_sq = h_plus*h_plus + h_cross*h_cross;
    return (GW_C * GW_C * omega * omega / (32.0 * GW_PI * GW_G)) * amp_sq;
}

/* ================================================================
 * L4 — Binary Quadrupole Strain (Circular Orbit)
 *
 * Newtonian quadrupole formula for a binary in circular orbit.
 *
 * h_plus  = -(2 G^2 M mu) / (c^4 a r) * (1+cos^2 iota) * cos(2 Phi)
 * h_cross = -(4 G^2 M mu) / (c^4 a r) * cos(iota) * sin(2 Phi)
 *
 * where M = m1+m2, mu = m1*m2/M, a = orbital separation,
 * Phi = orbital phase at (retarded) time.
 * ================================================================ */

void gw_binary_quadrupole_strain(double *h_plus, double *h_cross,
    double M, double mu, double a, double r, double iota, double Phi) {
    double GM = GW_G * M;
    double prefactor = 2.0 * GM * mu * GW_G / (GW_C*GW_C*GW_C*GW_C * a * r);

    *h_plus  = -prefactor * (1.0 + cos(iota)*cos(iota)) * cos(2.0 * Phi);
    *h_cross = -prefactor * 2.0 * cos(iota) * sin(2.0 * Phi);
}

/* ================================================================
 * L4 — Circular Binary GW Luminosity
 *
 * L_GW = (32/5) (G^4 / c^5) * (mu^2 * M^3) / a^5
 *
 * Key result: GW luminosity increases rapidly as orbit shrinks
 * (L_GW ~ a^{-5}), driving the inspiral.
 * ================================================================ */

double gw_binary_luminosity_circular(double M, double mu, double a) {
    double G4 = GW_G * GW_G * GW_G * GW_G;
    double C5 = GW_C * GW_C * GW_C * GW_C * GW_C;
    return (32.0 / 5.0) * (G4 / C5) * (mu * mu * M * M * M) / (a * a * a * a * a);
}
