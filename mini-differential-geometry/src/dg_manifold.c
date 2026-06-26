#include "dg_manifold.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

Manifold manifold_create(int dim, const char *name) {
    Manifold m;
    m.dim = dim;
    m.num_charts = 0;
    m.is_compact = false;
    m.is_oriented = true;
    m.has_boundary = false;
    if (name) {
        strncpy(m.name, name, DG_MAX_NAME - 1);
        m.name[DG_MAX_NAME - 1] = '\0';
    } else { m.name[0] = '\0'; }
    for (int i = 0; i < DG_MAX_CHARTS; i++) { m.charts[i].dim = 0; }
    return m;
}

int manifold_add_chart(Manifold *m, const Chart *c) {
    if (!m || !c) return -1;
    if (m->num_charts >= DG_MAX_CHARTS) return -1;
    if (c->dim != m->dim) return -1;
    m->charts[m->num_charts] = *c;
    m->num_charts++;
    return m->num_charts - 1;
}

void manifold_destroy(Manifold *m) {
    if (!m) return;
    m->num_charts = 0;
    m->dim = 0;
}

Chart chart_create(int dim, const char coord_names[][DG_MAX_COORD_NAME],
                   const double domain_min[], const double domain_max[]) {
    Chart c;
    c.dim = dim;
    for (int i = 0; i < dim && i < DG_MAX_DIM; i++) {
        if (coord_names[i][0] != '\0') {
            strncpy(c.coord_names[i], coord_names[i], DG_MAX_COORD_NAME - 1);
            c.coord_names[i][DG_MAX_COORD_NAME - 1] = '\0';
        } else { c.coord_names[i][0] = '\0'; }
        c.domain_min[i] = domain_min ? domain_min[i] : -1e100;
        c.domain_max[i] = domain_max ? domain_max[i] : 1e100;
        c.is_periodic[i] = false;
        c.period[i] = 0.0;
    }
    return c;
}

bool chart_contains_point(const Chart *c, const double x[]) {
    if (!c || !x) return false;
    for (int i = 0; i < c->dim && i < DG_MAX_DIM; i++) {
        if (c->is_periodic[i]) continue;
        if (x[i] < c->domain_min[i] || x[i] > c->domain_max[i]) return false;
    }
    return true;
}

/* Transition map callback storage */
typedef struct { void (*func)(const double*, double*, int); int active; } TransEntry;
static TransEntry trans_map[DG_MAX_CHARTS * DG_MAX_CHARTS];

int manifold_set_transition_map(Manifold *m, int src, int dst,
                                 transition_func_t func) {
    if (!m || !func) return -1;
    if (src < 0 || src >= m->num_charts || dst < 0 || dst >= m->num_charts) return -1;
    if (src == dst) return 0;
    int key = src * DG_MAX_CHARTS + dst;
    trans_map[key].func = func;
    trans_map[key].active = 1;
    return 0;
}

int manifold_transition_jacobian(const Manifold *m, int src, int dst,
                                  const double x_src[], double jacobian[]) {
    if (!m || !x_src || !jacobian) return -1;
    int dim = m->dim;
    if (src == dst) {
        for (int i = 0; i < dim * dim; i++) jacobian[i] = 0.0;
        for (int i = 0; i < dim; i++) jacobian[i * dim + i] = 1.0;
        return 0;
    }
    int key = src * DG_MAX_CHARTS + dst;
    if (!trans_map[key].active) return -1;
    double h = 1e-6;
    double f0[DG_MAX_DIM], fp[DG_MAX_DIM], xp[DG_MAX_DIM];
    for (int i = 0; i < dim; i++) xp[i] = x_src[i];
    trans_map[key].func(x_src, f0, dim);
    for (int j = 0; j < dim; j++) {
        for (int i = 0; i < dim; i++) xp[i] = x_src[i];
        xp[j] += h;
        trans_map[key].func(xp, fp, dim);
        for (int i = 0; i < dim; i++)
            jacobian[i * dim + j] = (fp[i] - f0[i]) / h;
    }
    return 0;
}

bool manifold_atlas_is_complete(const Manifold *m) {
    return m && m->num_charts >= 1;
}

TangentVector tangent_vector_create(int dim, const double components[],
                                     int chart_index) {
    TangentVector tv;
    tv.dim = dim;
    tv.chart_index = chart_index;
    for (int i = 0; i < dim && i < DG_MAX_DIM; i++)
        tv.components[i] = components ? components[i] : 0.0;
    return tv;
}

CotangentVector cotangent_vector_create(int dim, const double components[],
                                         int chart_index) {
    CotangentVector cv;
    cv.dim = dim;
    cv.chart_index = chart_index;
    for (int i = 0; i < dim && i < DG_MAX_DIM; i++)
        cv.components[i] = components ? components[i] : 0.0;
    return cv;
}

double cotangent_act_on_tangent(const CotangentVector *omega,
                                 const TangentVector *X) {
    if (!omega || !X) return 0.0;
    int dim = omega->dim < X->dim ? omega->dim : X->dim;
    if (dim > DG_MAX_DIM) dim = DG_MAX_DIM;
    double result = 0.0;
    for (int i = 0; i < dim; i++)
        result += omega->components[i] * X->components[i];
    return result;
}