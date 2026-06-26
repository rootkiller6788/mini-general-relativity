/**
 * @file    gw_core.c
 * @brief   GW Core — Tensor operations, antenna patterns, polarization
 *
 * L2 — TT projection, antenna pattern functions
 * L3 — 3x3 symmetric tensor algebra, eigenvalue computation
 * L4 — PN parameter, orbital separation, ISCO frequency
 */

#include "gw_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * L2 — Tensor Operations
 * ================================================================ */

void gw_tensor_zero(GwTensor3 *T) {
    memset(T, 0, sizeof(GwTensor3));
}

double gw_tensor_trace(const GwTensor3 *T) {
    return T->xx + T->yy + T->zz;
}

void gw_tensor_add(GwTensor3 *C, const GwTensor3 *A, const GwTensor3 *B) {
    C->xx = A->xx + B->xx; C->xy = A->xy + B->xy; C->xz = A->xz + B->xz;
    C->yx = A->yx + B->yx; C->yy = A->yy + B->yy; C->yz = A->yz + B->yz;
    C->zx = A->zx + B->zx; C->zy = A->zy + B->zy; C->zz = A->zz + B->zz;
}

void gw_tensor_sub(GwTensor3 *C, const GwTensor3 *A, const GwTensor3 *B) {
    C->xx = A->xx - B->xx; C->xy = A->xy - B->xy; C->xz = A->xz - B->xz;
    C->yx = A->yx - B->yx; C->yy = A->yy - B->yy; C->yz = A->yz - B->yz;
    C->zx = A->zx - B->zx; C->zy = A->zy - B->zy; C->zz = A->zz - B->zz;
}

void gw_tensor_scale(GwTensor3 *T, double alpha) {
    T->xx *= alpha; T->xy *= alpha; T->xz *= alpha;
    T->yx *= alpha; T->yy *= alpha; T->yz *= alpha;
    T->zx *= alpha; T->zy *= alpha; T->zz *= alpha;
}

double gw_tensor_contract(const GwTensor3 *A, const GwTensor3 *B) {
    return A->xx*B->xx + A->xy*B->xy + A->xz*B->xz
         + A->yx*B->yx + A->yy*B->yy + A->yz*B->yz
         + A->zx*B->zx + A->zy*B->zy + A->zz*B->zz;
}

void gw_tensor_tt_project(GwTensor3 *h_tt, const GwTensor3 *S,
                          double nx, double ny, double nz) {
    /* Transverse projector P_{ij} = delta_{ij} - n_i n_j */
    double P[3][3];
    P[0][0] = 1.0 - nx*nx; P[0][1] =     - nx*ny; P[0][2] =     - nx*nz;
    P[1][0] =     - ny*nx; P[1][1] = 1.0 - ny*ny; P[1][2] =     - ny*nz;
    P[2][0] =     - nz*nx; P[2][1] =     - nz*ny; P[2][2] = 1.0 - nz*nz;

    /* S matrix as 2D array */
    double Sm[3][3] = {
        {S->xx, S->xy, S->xz},
        {S->yx, S->yy, S->yz},
        {S->zx, S->zy, S->zz}
    };

    /* Lambda_{ij,kl} = P_{ik} P_{jl} - (1/2) P_{ij} P_{kl}
     * Apply: h_tt_{ij} = sum_{k,l} Lambda_{ij,kl} S_{kl}
     *         = sum_{k,l} P_{ik} P_{jl} S_{kl} - (1/2) P_{ij} sum_{k,l} P_{kl} S_{kl}
     *
     * First term: (P S P^T)_{ij} where S symmetric, so (P S P)_{ij}
     * Second term: (1/2) P_{ij} * Tr(P S)
     */

    /* Compute P*S: PS_{ik} = sum_l P_{il} S_{lk} */
    double PS[3][3];
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            PS[i][k] = 0.0;
            for (int l = 0; l < 3; l++) {
                PS[i][k] += P[i][l] * Sm[l][k];
            }
        }
    }

    /* (P*S)*P^T = (P*S)*P since P symmetric: PSP_{ij} = sum_k PS_{ik} P_{kj} */
    double PSP[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            PSP[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                PSP[i][j] += PS[i][k] * P[k][j];
            }
        }
    }

    /* Trace of P*S */
    double tr_PS = PS[0][0] + PS[1][1] + PS[2][2];

    /* h_tt_{ij} = PSP_{ij} - (1/2) P_{ij} * tr_PS */
    h_tt->xx = PSP[0][0] - 0.5 * P[0][0] * tr_PS;
    h_tt->xy = PSP[0][1] - 0.5 * P[0][1] * tr_PS;
    h_tt->xz = PSP[0][2] - 0.5 * P[0][2] * tr_PS;
    h_tt->yx = PSP[1][0] - 0.5 * P[1][0] * tr_PS;
    h_tt->yy = PSP[1][1] - 0.5 * P[1][1] * tr_PS;
    h_tt->yz = PSP[1][2] - 0.5 * P[1][2] * tr_PS;
    h_tt->zx = PSP[2][0] - 0.5 * P[2][0] * tr_PS;
    h_tt->zy = PSP[2][1] - 0.5 * P[2][1] * tr_PS;
    h_tt->zz = PSP[2][2] - 0.5 * P[2][2] * tr_PS;
}

/* ================================================================
 * L3 — 3x3 Symmetric Tensor Eigenvalue Decomposition
 *
 * Analytical method based on the cubic characteristic equation:
 *   det(T - lambda I) = -lambda^3 + I1 lambda^2 - I2 lambda + I3 = 0
 *
 * where I1 = Tr(T), I2 = (Tr(T)^2 - Tr(T^2))/2, I3 = det(T)
 *
 * Solved via trigonometric method for three real roots
 * (guaranteed for real symmetric matrices).
 * ================================================================ */

int gw_tensor_eigensystem(const GwTensor3 *T,
                          double eigenvals[3], double eigenvecs[9]) {
    /* Compute invariants */
    double Tmat[3][3] = {
        {T->xx, T->xy, T->xz},
        {T->yx, T->yy, T->yz},
        {T->zx, T->zy, T->zz}
    };

    double I1 = T->xx + T->yy + T->zz;

    double T2_trace = T->xx*T->xx + T->xy*T->xy + T->xz*T->xz
                    + T->yx*T->yx + T->yy*T->yy + T->yz*T->yz
                    + T->zx*T->zx + T->zy*T->zy + T->zz*T->zz;
    double I2 = 0.5 * (I1*I1 - T2_trace);

    double det = T->xx*(T->yy*T->zz - T->yz*T->zy)
               - T->xy*(T->yx*T->zz - T->yz*T->zx)
               + T->xz*(T->yx*T->zy - T->yy*T->zx);
    double I3 = det;

    /* Shift to zero-mean: p = I2 - I1^2/3, q = (2 I1^3)/27 - I1 I2/3 + I3 */
    double p = I2 - I1*I1/3.0;
    double q = (2.0*I1*I1*I1)/27.0 - I1*I2/3.0 + I3;

    /* For symmetric matrices: p <= 0, three real roots */
    if (p > 0) return -1;  /* Should not happen for symmetric real */

    double phi;
    double sqrt_mp3 = sqrt(-p / 3.0);
    double arg = q / (2.0 * sqrt_mp3 * sqrt_mp3 * sqrt_mp3);

    /* Clamp to [-1, 1] for numerical safety */
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;

    phi = acos(arg) / 3.0;
    double shift = I1 / 3.0;

    /* Three real eigenvalues (descending order) */
    double lam1 = 2.0 * sqrt_mp3 * cos(phi) + shift;
    double lam2 = 2.0 * sqrt_mp3 * cos(phi + 2.0*GW_PI/3.0) + shift;
    double lam3 = 2.0 * sqrt_mp3 * cos(phi + 4.0*GW_PI/3.0) + shift;

    /* Sort descending */
    double lams[3] = {lam1, lam2, lam3};
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            if (lams[i] < lams[j]) {
                double tmp = lams[i]; lams[i] = lams[j]; lams[j] = tmp;
            }
        }
    }
    eigenvals[0] = lams[0];
    eigenvals[1] = lams[1];
    eigenvals[2] = lams[2];

    /* Compute eigenvectors via inverse iteration for each eigenvalue.
     * For brevity, compute using the cofactor method.
     * For each eigenvalue lambda, solve (T - lambda I) v = 0.
     * Use the two largest rows of the adjugate matrix. */
    for (int ie = 0; ie < 3; ie++) {
        double lam = eigenvals[ie];

        /* Build (T - lam I) */
        double M[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                M[i][j] = Tmat[i][j];
        M[0][0] -= lam; M[1][1] -= lam; M[2][2] -= lam;

        /* Compute cofactors of first row */
        double c0 = M[1][1]*M[2][2] - M[1][2]*M[2][1];
        double c1 = M[1][2]*M[2][0] - M[1][0]*M[2][2];
        double c2 = M[1][0]*M[2][1] - M[1][1]*M[2][0];

        double norm = sqrt(c0*c0 + c1*c1 + c2*c2);
        if (norm < 1e-30) {
            /* Fallback: use second row */
            c0 = M[0][2]*M[2][1] - M[0][1]*M[2][2];
            c1 = M[0][0]*M[2][2] - M[0][2]*M[2][0];
            c2 = M[0][1]*M[2][0] - M[0][0]*M[2][1];
            norm = sqrt(c0*c0 + c1*c1 + c2*c2);
        }
        if (norm < 1e-30) {
            /* Fallback: identity-ish */
            c0 = (ie==0) ? 1.0 : 0.0;
            c1 = (ie==1) ? 1.0 : 0.0;
            c2 = (ie==2) ? 1.0 : 0.0;
            norm = 1.0;
        }

        eigenvecs[ie*3 + 0] = c0 / norm;
        eigenvecs[ie*3 + 1] = c1 / norm;
        eigenvecs[ie*3 + 2] = c2 / norm;
    }

    return 0;
}

/* ================================================================
 * L2 — GW Strain Memory Management
 * ================================================================ */

int gw_strain_alloc(GwStrainSeries *s, size_t n) {
    if (n == 0) return -1;
    s->n = n;
    s->t  = (double*)calloc(n, sizeof(double));
    s->hp = (double*)calloc(n, sizeof(double));
    s->hx = (double*)calloc(n, sizeof(double));
    if (!s->t || !s->hp || !s->hx) {
        gw_strain_free(s);
        return -1;
    }
    return 0;
}

void gw_strain_free(GwStrainSeries *s) {
    free(s->t);  s->t  = NULL;
    free(s->hp); s->hp = NULL;
    free(s->hx); s->hx = NULL;
    s->n = 0;
}

double gw_detector_strain(double F_plus, double F_cross,
                          double h_plus, double h_cross) {
    return F_plus * h_plus + F_cross * h_cross;
}

/* ================================================================
 * L2 — Antenna Pattern Functions
 * ================================================================ */

void gw_antenna_pattern(double *F_plus, double *F_cross,
                        double theta, double phi, double psi) {
    double cos_theta = cos(theta);
    (void)theta; /* sin_theta not needed in final formula */
    double cos_2phi  = cos(2.0 * phi);
    double sin_2phi  = sin(2.0 * phi);
    double cos_2psi  = cos(2.0 * psi);
    double sin_2psi  = sin(2.0 * psi);
    double one_plus_cos2 = 1.0 + cos_theta * cos_theta;

    *F_plus  = 0.5 * one_plus_cos2 * cos_2phi * cos_2psi
             - cos_theta * sin_2phi * sin_2psi;
    *F_cross = 0.5 * one_plus_cos2 * cos_2phi * sin_2psi
             + cos_theta * sin_2phi * cos_2psi;
}

void gw_antenna_pattern_lisa(double *F_plus, double *F_cross,
                             double theta, double phi, double psi) {
    /* LISA: 60-degree arms.
     * Approximate pattern (averaged over LISA orbit): similar form
     * but with different angular dependence.
     */
    double cos_theta = cos(theta);
    double cos_2phi  = cos(2.0 * phi);
    double sin_2phi  = sin(2.0 * phi);
    double cos_2psi  = cos(2.0 * psi);
    double sin_2psi  = sin(2.0 * psi);
    double fac = 0.5 * (1.0 + cos_theta * cos_theta);

    /* LISA has sqrt(3)/2 factor from 60 deg opening vs 90 deg */
    double pref = 0.8660254037844386;  /* sqrt(3)/2 */

    *F_plus  = pref * (fac * cos_2phi * cos_2psi - cos_theta * sin_2phi * sin_2psi);
    *F_cross = pref * (fac * cos_2phi * sin_2psi + cos_theta * sin_2phi * cos_2psi);
}

/* ================================================================
 * L2 — Polarization Tensors
 * ================================================================ */

void gw_pol_tensor_plus(GwTensor3 *e_plus) {
    e_plus->xx =  1.0;
    e_plus->xy =  0.0;
    e_plus->xz =  0.0;
    e_plus->yx =  0.0;
    e_plus->yy = -1.0;
    e_plus->yz =  0.0;
    e_plus->zx =  0.0;
    e_plus->zy =  0.0;
    e_plus->zz =  0.0;
}

void gw_pol_tensor_cross(GwTensor3 *e_cross) {
    e_cross->xx = 0.0;
    e_cross->xy = 1.0;
    e_cross->xz = 0.0;
    e_cross->yx = 1.0;
    e_cross->yy = 0.0;
    e_cross->yz = 0.0;
    e_cross->zx = 0.0;
    e_cross->zy = 0.0;
    e_cross->zz = 0.0;
}

void gw_pol_rotate(GwTensor3 *ep_out, GwTensor3 *ex_out, double psi) {
    double c2 = cos(2.0 * psi);
    double s2 = sin(2.0 * psi);

    GwTensor3 ep, ex;
    gw_pol_tensor_plus(&ep);
    gw_pol_tensor_cross(&ex);

    /* e'^+ = c2 * e^+ + s2 * e^x */
    gw_tensor_scale(&ep, c2);
    GwTensor3 tmp;
    gw_tensor_zero(&tmp);
    tmp = ex;
    gw_tensor_scale(&tmp, s2);
    gw_tensor_add(ep_out, &ep, &tmp);

    /* e'^x = -s2 * e^+ + c2 * e^x */
    gw_pol_tensor_plus(&ep);
    gw_pol_tensor_cross(&ex);
    gw_tensor_scale(&ep, -s2);
    gw_tensor_scale(&ex, c2);
    gw_tensor_add(ex_out, &ep, &ex);
}

void gw_strain_decompose(const GwTensor3 *h, double *h_plus, double *h_cross) {
    /* For propagation along z-hat:
     * h_+ = (h_xx - h_yy) / 2
     * h_x = h_xy
     */
    *h_plus  = 0.5 * (h->xx - h->yy);
    *h_cross = h->xy;
}

/* ================================================================
 * L4 — Fundamental GW Scales
 * ================================================================ */

double gw_pn_parameter(double M_total, double f_gw) {
    double x = GW_PI * GW_G * M_total * f_gw / (GW_C * GW_C * GW_C);
    return pow(x, 2.0/3.0);
}

double gw_orbital_separation(double M_total, double f_gw) {
    double omega = GW_PI * f_gw;
    return pow(GW_G * M_total / (omega * omega), 1.0/3.0);
}

double gw_isco_frequency(double M_total) {
    /* Schwarzschild ISCO: r_ISCO = 6 G M / c^2
     * Kepler frequency: omega^2 = G M / r^3
     * f_ISCO = omega / (2 pi) = c^3 / (6^{3/2} 2 pi G M)
     * Actually: f_gw = 2 * f_orb for quadrupole, and at ISCO
     * f_gw ~ c^3 / (6^{3/2} pi G M) (approximately)
     */
    double six_pow = pow(6.0, 1.5);  /* 6^{3/2} ~ 14.697 */
    return GW_C * GW_C * GW_C / (six_pow * GW_PI * GW_G * M_total);
}

double gw_schwarzschild_radius(double M) {
    return 2.0 * GW_G * M / (GW_C * GW_C);
}

double gw_lightring_frequency(double M_total) {
    double three_pow = pow(3.0, 1.5);  /* 3^{3/2} ~ 5.196 */
    return GW_C * GW_C * GW_C / (three_pow * 2.0 * GW_PI * GW_G * M_total);
}
