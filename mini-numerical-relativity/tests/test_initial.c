/**
 * @file test_initial.c
 * @brief Tests for initial data: Schwarzschild, Brill-Lindquist, Bowen-York.
 *
 * Knowledge: L4 (constraints), L6 (canonical initial data sets)
 */

#include "../include/nr_initial.h"
#include "../include/nr_adm.h"
#include "../include/nr_grid.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0, tests_passed = 0;

#define TEST(n) do { tests_run++; printf("  TEST: %s ... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define CHECK(c) do { if(!(c)){printf("FAILED at %s:%d: %s\n",__FILE__,__LINE__,#c);return;} } while(0)

static void test_schwarzschild_isotropic(void) {
    TEST("Schwarzschild isotropic initial data");
    int N = 12;
    nr_grid_t grid;
    nr_grid_init(&grid, N, N, N, -4.0, 4.0, -4.0, 4.0, -4.0, 4.0, 3);

    nr_adm_state_t *state = nr_adm_alloc(N, N, N, 3);
    CHECK(state != NULL);

    double M = 1.0;
    nr_init_schwarzschild_isotropic(state, &grid, M);

    /* Check conformal factor at r=2: ψ = 1 + M/(2r) = 1 + 1/4 = 1.25 */
    double x = 2.0, y = 0.0, z = 0.0;
    double psi_expected = 1.0 + M / (2.0 * sqrt(x*x + y*y + z*z));

    /* Get γ_xx at a point near (2,0,0) */
    int ng = 3;
    int i = (int)((x - grid.xmin) / grid.dx) + ng;
    int j = (int)((y - grid.ymin) / grid.dx) + ng;
    int k = (int)((z - grid.zmin) / grid.dx) + ng;

    if (i >= ng && i < N+ng && j >= ng && j < N+ng && k >= ng && k < N+ng) {
        double gxx = state->gamma->comp[0]->data[
            i*state->gamma->comp[0]->stride_x
            + j*state->gamma->comp[0]->stride_y + k];
        double psi_numerical = pow(gxx, 0.25);
        printf("  ψ_num=%f, ψ_expected=%f\n", psi_numerical, psi_expected);
        CHECK(fabs(psi_numerical - psi_expected) < 0.1);
    }

    /* Extrinsic curvature should be zero */
    double Kxx = state->K->comp[0]->data[
        (N/2+ng)*state->K->comp[0]->stride_x
        + (N/2+ng)*state->K->comp[0]->stride_y + (N/2+ng)];
    CHECK(fabs(Kxx) < 1e-10);

    nr_adm_free(state);
    nr_grid_free(&grid);
    PASS();
}

static void test_brill_lindquist(void) {
    TEST("Brill-Lindquist binary BH data");
    int N = 12;
    nr_grid_t grid;
    nr_grid_init(&grid, N, N, N, -4.0, 4.0, -4.0, 4.0, -4.0, 4.0, 3);

    nr_adm_state_t *state = nr_adm_alloc(N, N, N, 3);
    CHECK(state != NULL);

    double masses[2] = {0.5, 0.5};
    double pos_x[2] = {1.0, -1.0};
    double pos_y[2] = {0.0, 0.0};
    double pos_z[2] = {0.0, 0.0};

    nr_init_brill_lindquist(state, &grid, 2, masses, pos_x, pos_y, pos_z);

    /* The 3-metric should be conformally flat with ψ > 1 */
    int ng = 3;
    double gxx = state->gamma->comp[0]->data[
        (N/2+ng)*state->gamma->comp[0]->stride_x
        + (N/2+ng)*state->gamma->comp[0]->stride_y + (N/2+ng)];
    CHECK(gxx > 1.0);  /* ψ^4 > 1 for nontrivial data */

    /* Off-diagonal should be zero */
    double gxy = state->gamma->comp[1]->data[
        (N/2+ng)*state->gamma->comp[1]->stride_x
        + (N/2+ng)*state->gamma->comp[1]->stride_y + (N/2+ng)];
    CHECK(fabs(gxy) < 1e-10);

    nr_adm_free(state);
    nr_grid_free(&grid);
    PASS();
}

static void test_schwarzschild_trumpet(void) {
    TEST("Schwarzschild trumpet BSSN data");
    int N = 12;
    nr_grid_t grid;
    nr_grid_init(&grid, N, N, N, -4.0, 4.0, -4.0, 4.0, -4.0, 4.0, 3);

    nr_bssn_state_t *state = nr_bssn_alloc(N, N, N, 3);
    CHECK(state != NULL);

    nr_init_schwarzschild_trumpet(state, &grid, 1.0);

    /* At the puncture (origin), lapse should be near 0 */
    int ng = 3;
    double alpha_center = state->alpha->data[
        (N/2+ng)*state->alpha->stride_x
        + (N/2+ng)*state->alpha->stride_y + (N/2+ng)];
    printf("  alpha(0) = %f\n", alpha_center);
    CHECK(alpha_center < 0.5);  /* Singularity avoidance */

    /* Far away, lapse should approach 1 (small domain may not reach 1) */
    double alpha_far = state->alpha->data[
        (N-1+ng)*state->alpha->stride_x
        + (N/2+ng)*state->alpha->stride_y + (N/2+ng)];
    printf("  alpha(r_max) = %f\n", alpha_far);
    CHECK(alpha_far > 0.5);  /* Should be recovering toward 1 */

    nr_bssn_free(state);
    nr_grid_free(&grid);
    PASS();
}

int main(void) {
    printf("=== test_initial ===\n");
    test_schwarzschild_isotropic();
    test_brill_lindquist();
    test_schwarzschild_trumpet();
    printf("\nResults: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
