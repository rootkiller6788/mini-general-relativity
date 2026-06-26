/**
 * @file test_grid.c
 * @brief Tests for grid infrastructure: allocation, access, coordinates.
 *
 * Knowledge covered: L1 (grid function), L5 (memory layout)
 */

#include "../include/nr_grid.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define CHECK(cond) do { if (!(cond)) { printf("FAILED at %s:%d: %s\n", __FILE__, __LINE__, #cond); return; } } while(0)

static void test_alloc_free(void) {
    TEST("gf_alloc / gf_free");
    nr_gf_t *gf = nr_gf_alloc(10, 10, 10, 2);
    CHECK(gf != NULL);
    CHECK(gf->data != NULL);
    CHECK(gf->nx == 10);
    CHECK(gf->ng == 2);
    double *p = nr_gf_at(gf, 5, 5, 5);
    CHECK(p != NULL);
    *p = 3.14;
    CHECK(*nr_gf_at(gf, 5, 5, 5) == 3.14);
    nr_gf_free(gf);
    PASS();
}

static void test_set_and_copy(void) {
    TEST("set_all / copy / l2norm");
    nr_gf_t *a = nr_gf_alloc(8, 8, 8, 1);
    nr_gf_t *b = nr_gf_alloc(8, 8, 8, 1);
    CHECK(a && b);

    nr_gf_set_all(a, 5.0);
    CHECK(fabs(nr_gf_l2norm(a) - 5.0) < 1e-10);

    nr_gf_copy(b, a);
    CHECK(fabs(nr_gf_l2norm(b) - 5.0) < 1e-10);

    /* Modify a, check b unchanged */
    *nr_gf_at(a, 2, 2, 2) = 0.0;
    CHECK(fabs(*nr_gf_at(b, 2, 2, 2) - 5.0) < 1e-10);

    double maxv = nr_gf_maxabs(a);
    CHECK(maxv >= 5.0);

    nr_gf_free(a); nr_gf_free(b);
    PASS();
}

static void test_grid_init(void) {
    TEST("grid_init / grid_r");
    nr_grid_t grid;
    nr_grid_init(&grid, 20, 20, 20, -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, 2);

    CHECK(fabs(grid.dx - 10.0/19.0) < 1e-10);
    CHECK(fabs(grid.xmin - (-5.0)) < 1e-10);
    CHECK(fabs(grid.xmax - 5.0) < 1e-10);

    /* Check coordinate arrays — even number of points, center not exactly 0 */
    int center = 20/2 + 2; /* interior center index */
    double expected = -5.0 + 10 * (10.0/19.0); /* actual coordinate */
    CHECK(fabs(nr_grid_x(&grid, center) - expected) < 1e-10);
    CHECK(fabs(nr_grid_y(&grid, center) - expected) < 1e-10);
    CHECK(fabs(nr_grid_z(&grid, center) - expected) < 1e-10);

    /* Radial coordinate */
    double r = nr_grid_r(&grid, 20+2, 10+2, 10+2);
    /* x = 5, y = 0, z = 0 → r = 5 */
    CHECK(fabs(r - 5.0) < 0.6);

    /* Spherical coordinates */
    double rr, st, ct, sp, cp;
    nr_grid_spherical(&grid, 22, 12, 12, &rr, &st, &ct, &sp, &cp);
    /* Should have sensible values */
    CHECK(rr > 0);

    nr_grid_free(&grid);
    PASS();
}

static void test_tensor_vector_alloc(void) {
    TEST("sym_tensor3 / vector3 alloc");
    nr_sym_tensor3_t *t = nr_sym_tensor3_alloc(10, 10, 10, 2);
    nr_vector3_t *v = nr_vector3_alloc(10, 10, 10, 2);
    CHECK(t != NULL);
    CHECK(v != NULL);
    for (int c = 0; c < 6; c++) CHECK(t->comp[c] != NULL);
    for (int c = 0; c < 3; c++) CHECK(v->comp[c] != NULL);

    /* Set and check a component */
    nr_gf_set_all(t->comp[0], 3.0);  /* γ_xx */
    CHECK(fabs(nr_gf_l2norm(t->comp[0]) - 3.0) < 1e-10);

    nr_gf_set_all(v->comp[0], 1.5);  /* β^x */
    CHECK(fabs(nr_gf_l2norm(v->comp[0]) - 1.5) < 1e-10);

    nr_sym_tensor3_free(t);
    nr_vector3_free(v);
    PASS();
}

int main(void) {
    printf("=== test_grid ===\n");
    test_alloc_free();
    test_set_and_copy();
    test_grid_init();
    test_tensor_vector_alloc();
    printf("\nResults: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
