/**
 * @file    test_core.c
 * @brief   Unit tests for GW core: tensors, antenna patterns, polarization
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "gw_core.h"

static int passed = 0, failed = 0;
#define T(n)  printf("  TEST: %s ... ", n)
#define P()   do { printf("PASS\n"); passed++; } while(0)
#define F(m)  do { printf("FAIL: %s\n", m); failed++; } while(0)
#define C(c,m) do { if(!(c)){F(m);return;} } while(0)

static void test_tensor_ops(void) {
    T("tensor operations");
    GwTensor3 A = {1,2,3, 2,4,5, 3,5,6};
    GwTensor3 B, C;

    gw_tensor_zero(&B);
    C(B.xx == 0 && B.yy == 0 && B.zz == 0, "zero");

    double tr = gw_tensor_trace(&A);
    C(fabs(tr - 11.0) < 1e-10, "trace");

    gw_tensor_scale(&A, 2.0);
    C(fabs(gw_tensor_trace(&A) - 22.0) < 1e-10, "scale");

    double contr = gw_tensor_contract(&A, &A);
    C(contr > 0, "self-contract positive");
    P();
}

static void test_tt_projection(void) {
    T("TT gauge projection");

    /* Create a plus-polarized wave propagating in z-direction */
    GwTensor3 S = {1,0,0, 0,-1,0, 0,0,0};
    GwTensor3 h_tt;

    gw_tensor_tt_project(&h_tt, &S, 0.0, 0.0, 1.0);

    /* For n = (0,0,1): P = diag(1,1,0)
     * P S P = S (since S only has xx,yy,xy components)
     * tr(P S) = 1 + (-1) + 0 = 0
     * h_tt = PSP - 0.5*P*0 = PSP = S
     */
    C(fabs(h_tt.xx - 1.0) < 1e-10, "tt_xx");
    C(fabs(h_tt.yy + 1.0) < 1e-10, "tt_yy");
    C(fabs(h_tt.zz) < 1e-10, "tt_zz");
    C(fabs(h_tt.xy) < 1e-10, "tt_xy");
    P();
}

static void test_antenna_pattern(void) {
    T("antenna pattern functions");

    double Fp, Fc;
    /* Overhead source: theta=0, phi=0, psi=0 */
    gw_antenna_pattern(&Fp, &Fc, 0.0, 0.0, 0.0);
    /* F_plus = 0.5*(1+1)*1*1 - 0 = 1.0
     * F_cross = 0.5*(1+1)*1*0 + 0 = 0.0 */
    C(fabs(Fp - 1.0) < 1e-10, "F_plus overhead");
    C(fabs(Fc) < 1e-10, "F_cross overhead");

    /* Horizon source: theta=pi/2, phi=0, psi=0 */
    gw_antenna_pattern(&Fp, &Fc, GW_PI/2.0, 0.0, 0.0);
    /* F_plus = 0.5 * (1+0) * 1 * 1 - 0 = 0.5 */
    C(fabs(Fp - 0.5) < 1e-10, "F_plus horizon");
    P();
}

static void test_polarization_tensors(void) {
    T("polarization tensors");
    GwTensor3 ep, ex;
    gw_pol_tensor_plus(&ep);
    gw_pol_tensor_cross(&ex);

    C(ep.xx == 1.0 && ep.yy == -1.0, "e_plus diagonal");
    C(ex.xy == 1.0 && ex.yx == 1.0, "e_cross off-diagonal");

    /* Decompose */
    GwTensor3 h = {2.0, 3.0, 0.0, 3.0, -4.0, 0.0, 0.0, 0.0, 0.0};
    double hp, hx;
    gw_strain_decompose(&h, &hp, &hx);
    C(fabs(hp - 3.0) < 1e-10, "h_plus = (xx-yy)/2 = 3");
    C(fabs(hx - 3.0) < 1e-10, "h_cross = xy = 3");
    P();
}

static void test_pn_parameter(void) {
    T("PN parameter and scales");
    double M = 60.0 * GW_MSUN;
    double f = 50.0;

    double x = gw_pn_parameter(M, f);
    C(x > 0 && x < 1.0, "PN parameter in (0,1)");

    double a = gw_orbital_separation(M, f);
    C(a > 0, "orbital separation positive");

    double f_isco = gw_isco_frequency(M);
    C(f_isco > 0, "ISCO frequency positive");

    double Rs = gw_schwarzschild_radius(M);
    C(Rs > 0, "Schwarzschild radius positive");
    C(a > Rs, "orbit outside Schwarzschild radius");
    P();
}

static void test_strain_series(void) {
    T("strain series allocation");
    GwStrainSeries s = {0};
    int ret = gw_strain_alloc(&s, 100);
    C(ret == 0, "allocation");
    C(s.n == 100, "n set");
    C(s.t != NULL && s.hp != NULL && s.hx != NULL, "buffers allocated");
    gw_strain_free(&s);
    C(s.n == 0, "freed");
    C(s.t == NULL, "t freed");
    P();
}

int main(void) {
    printf("=== gw_core Tests ===\n\n");
    test_tensor_ops();
    test_tt_projection();
    test_antenna_pattern();
    test_polarization_tensors();
    test_pn_parameter();
    test_strain_series();
    printf("\nResults: %d passed, %d failed\n", passed, failed);
    return (failed > 0) ? 1 : 0;
}
