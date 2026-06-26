/**
 * @file nr_bssn.h
 * @brief BSSN (Baumgarte-Shapiro-Shibata-Nakamura) formulation
 *
 * The BSSN formulation is the standard evolution system used in modern
 * numerical relativity. It reformulates the ADM equations via a conformal
 * decomposition that dramatically improves numerical stability.
 *
 * Conformal decomposition (BSSN variables):
 *   φ     = (1/12) ln(det γ_{ij})        conformal factor
 *   γ̃_{ij} = e^{-4φ} γ_{ij}              conformal metric (det = 1)
 *   K     = γ^{ij} K_{ij}                trace of extrinsic curvature
 *   Ã_{ij} = e^{-4φ} (K_{ij} - γ_{ij}K/3) trace-free conformal curvature
 *   Γ̃^i  = γ̃^{jk} Γ̃^i_{jk} = -∂_j γ̃^{ij}  conformal connection functions
 *
 * Gauge conditions:
 *   1+log slicing:  ∂_t α = -2αK + β^i ∂_i α
 *   Gamma-driver:   ∂_t β^i = (3/4) B^i
 *                   ∂_t B^i = ∂_t Γ̃^i - η B^i
 *
 * References:
 *   - Shibata & Nakamura, PRD 52, 5428 (1995)
 *   - Baumgarte & Shapiro, PRD 59, 024007 (1998)
 *   - Baumgarte & Shapiro (2010), Ch. 11
 *   - Alcubierre (2008), Ch. 6
 *   - Campanelli et al., PRL 96, 111101 (2006)
 *
 * Knowledge Coverage:
 *   L1: BSSN variables φ, γ̃_{ij}, K, Ã_{ij}, Γ̃^i
 *   L2: Conformal decomposition, constraint addition, numerical stability
 *   L3: Conformal transformations, trace-free decomposition
 *   L4: BSSN evolution equations
 *   L5: Finite-difference RHS, Kreiss-Oliger dissipation
 */

#ifndef NR_BSSN_H
#define NR_BSSN_H

#include "nr_grid.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
 * L1: BSSN State Variables
 *
 * The BSSN system evolves 17 dynamical variables:
 *   1 (φ) + 5 (γ̃_{ij}, det=1 → independent) + 1 (K) + 5 (Ã_{ij}, trace-free)
 *   + 3 (Γ̃^i) + 1 (α) + 1 (β^i=3, B^i=3) = 17
 * =========================================================================== */

/**
 * L1: BSSN dynamical state on a 3D grid.
 *
 * Evolved variables:
 *   phi       — conformal factor φ
 *   gt[6]     — conformal metric γ̃_{ij} (6 stored, 5 independent, det enforced)
 *   trK       — trace of extrinsic curvature K
 *   At[6]     — trace-free conformal extrinsic curvature Ã_{ij} (6 stored, 5 independent)
 *   Gt[3]     — conformal connection functions Γ̃^i
 *
 * Gauge variables:
 *   alpha     — lapse α
 *   beta[3]   — shift β^i
 *   B[3]      — auxiliary variable for Gamma-driver
 *
 * Diagnostic (not evolved):
 *   chi       — alternative conformal factor χ = e^{-4φ}
 */
typedef struct {
    nr_gf_t      *phi;        /**< Conformal factor φ = (1/12) ln(det γ_{ij}) */
    nr_gf_t      *chi;        /**< Alternative: χ = exp(-4φ) */
    nr_sym_tensor3_t *gt;     /**< Conformal metric γ̃_{ij} (6 comp stored) */
    nr_gf_t      *trK;        /**< Trace of extrinsic curvature K */
    nr_sym_tensor3_t *At;     /**< Trace-free conformal curvature Ã_{ij} */
    nr_vector3_t *Gt;         /**< Conformal connection Γ̃^i */
    nr_gf_t      *alpha;      /**< Lapse function α */
    nr_vector3_t *beta;       /**< Shift vector β^i */
    nr_vector3_t *B;          /**< Auxiliary Gamma-driver variable B^i */
} nr_bssn_state_t;

/**
 * BSSN right-hand side (time derivatives).
 * Same structure as state but holds ∂_t of each variable.
 */
typedef struct {
    nr_gf_t      *dt_phi;
    nr_gf_t      *dt_chi;
    nr_sym_tensor3_t *dt_gt;
    nr_gf_t      *dt_trK;
    nr_sym_tensor3_t *dt_At;
    nr_vector3_t *dt_Gt;
    nr_gf_t      *dt_alpha;
    nr_vector3_t *dt_beta;
    nr_vector3_t *dt_B;
} nr_bssn_rhs_t;

/* ===========================================================================
 * L1: BSSN State Allocation
 * =========================================================================== */

/**
 * Allocate a full BSSN state (17 evolved + 3 gauge auxiliary grid functions).
 * Complexity: O(N^3) memory.
 */
nr_bssn_state_t* nr_bssn_alloc(int nx, int ny, int nz, int ng);

/** Allocate BSSN RHS matching the given state dimensions */
nr_bssn_rhs_t* nr_bssn_rhs_alloc(int nx, int ny, int nz, int ng);

/** Free a BSSN state */
void nr_bssn_free(nr_bssn_state_t *state);

/** Free a BSSN RHS structure */
void nr_bssn_rhs_free(nr_bssn_rhs_t *rhs);

/**
 * Set flat (Minkowski) initial BSSN data:
 *   φ = 0, γ̃_{ij} = δ_{ij}, K = 0, Ã_{ij} = 0, Γ̃^i = 0
 *   α = 1, β^i = 0, B^i = 0
 * Complexity: O(N^3).
 */
void nr_bssn_set_minkowski(nr_bssn_state_t *state);

/* ===========================================================================
 * L2: ADM → BSSN Conversion
 * =========================================================================== */

/**
 * Convert ADM variables (γ_{ij}, K_{ij}) to BSSN variables.
 *
 * φ = (1/12) ln(det γ)
 * γ̃_{ij} = e^{-4φ} γ_{ij}
 * K = γ^{ij} K_{ij}
 * Ã_{ij} = e^{-4φ} (K_{ij} - (1/3) γ_{ij} K)
 * Γ̃^i = -∂_j γ̃^{ij}
 *
 * Complexity: O(N^3).
 * Reference: Baumgarte & Shapiro (2010), Box 11.1
 */
void nr_bssn_from_adm(const nr_sym_tensor3_t *gamma,
                       const nr_sym_tensor3_t *K,
                       nr_bssn_state_t *bssn, double dx);

/**
 * Convert BSSN variables back to ADM variables.
 *
 * γ_{ij} = e^{4φ} γ̃_{ij}
 * K_{ij} = e^{4φ} (Ã_{ij} + (1/3) γ̃_{ij} K)
 *
 * Complexity: O(N^3).
 * Reference: Baumgarte & Shapiro (2010), Box 11.1
 */
void nr_adm_from_bssn(const nr_bssn_state_t *bssn,
                       nr_sym_tensor3_t *gamma,
                       nr_sym_tensor3_t *K);

/* ===========================================================================
 * L4: BSSN Evolution Equations (RHS)
 *
 * The BSSN evolution system with standard gauge:
 *
 * ∂_t φ = -(1/6) α K + β^k ∂_k φ + (1/6) ∂_k β^k
 * ∂_t γ̃_{ij} = -2α Ã_{ij} + β^k ∂_k γ̃_{ij}
 *              + γ̃_{ik} ∂_j β^k + γ̃_{kj} ∂_i β^k - (2/3) γ̃_{ij} ∂_k β^k
 * ∂_t K = -D^i D_i α + α(Ã_{ij} Ã^{ij} + K²/3) + β^i ∂_i K
 * ∂_t Ã_{ij} = [e^{-4φ}(-D_i D_j α + α R_{ij})]^{TF}
 *              + α(K Ã_{ij} - 2 Ã_{ik} Ã^k_j)
 *              + β^k ∂_k Ã_{ij} + Ã_{ik} ∂_j β^k + Ã_{kj} ∂_i β^k
 *              - (2/3) Ã_{ij} ∂_k β^k
 * ∂_t Γ̃^i = -2 Ã^{ij} ∂_j α
 *            + 2α(Γ̃^i_{jk} Ã^{kj} + 6 Ã^{ij} ∂_j φ - (2/3) γ̃^{ij} ∂_j K)
 *            + γ̃^{jk} ∂_j ∂_k β^i + (1/3) γ̃^{ik} ∂_k ∂_j β^j
 *            + β^j ∂_j Γ̃^i - Γ̃^j ∂_j β^i + (2/3) Γ̃^i ∂_j β^j
 *
 * Reference: Baumgarte & Shapiro (2010), Eqs. (11.4.9)-(11.4.11), (11.4.3), (11.4.5)
 * =========================================================================== */

/**
 * Compute inverse conformal metric γ̃^{ij} from γ̃_{ij} at a grid point.
 * Returns determinant.
 *
 * Complexity: O(1).
 */
double nr_bssn_inv_conformal_metric(const nr_sym_tensor3_t *gt,
                                     int i, int j, int k,
                                     double inv_gt[3][3]);

/**
 * Compute the BSSN RHS at a single grid point.
 * Fills all time derivatives in rhs.
 *
 * Complexity: O(1).
 */
void nr_bssn_rhs_point(const nr_bssn_state_t *state,
                        int i, int j, int k, double dx,
                        nr_bssn_rhs_t *rhs);

/**
 * Compute BSSN RHS over the entire interior grid.
 *
 * Complexity: O(N^3).
 */
void nr_bssn_rhs(const nr_bssn_state_t *state, double dx,
                 nr_bssn_rhs_t *rhs);

/* ===========================================================================
 * L2: BSSN Constraint Addition
 *
 * The BSSN system replaces ∂_t Γ̃^i by the momentum constraint:
 * Instead of evolving Γ̃^i from its definition, add the momentum constraint.
 *
 * ∂_t Γ̃^i = [RHS from definition] + (constraint_terms)
 *
 * This "constraint damping" is essential for numerical stability.
 * =========================================================================== */

/**
 * Enforce algebraic constraints on BSSN variables:
 *   det(γ̃_{ij}) = 1
 *   tr(Ã_{ij}) = 0
 *
 * At each interior point, rescale γ̃_{ij} and subtract trace from Ã_{ij}.
 *
 * Complexity: O(N^3).
 * Reference: Baumgarte & Shapiro (2010), Sec. 11.4
 */
void nr_bssn_enforce_constraints(nr_bssn_state_t *state);

/**
 * Enforce the Gamma constraint Γ̃^i = -∂_j γ̃^{ij}.
 * Recomputes Γ̃^i from the current γ̃_{ij} at all interior points.
 * Useful after enforcing det(γ̃) = 1.
 *
 * Complexity: O(N^3).
 */
void nr_bssn_recompute_Gamma(nr_bssn_state_t *state, double dx);

/* ===========================================================================
 * L2: Constraint Monitoring
 *
 * The BSSN system preserves three types of constraints:
 *   1. Hamiltonian constraint H
 *   2. Momentum constraint M_i
 *   3. Gamma constraint: G^i = Γ̃^i + ∂_j γ̃^{ij} = 0
 * =========================================================================== */

/**
 * Compute the Hamiltonian constraint in BSSN variables.
 *
 * H = R + (2/3)K² - Ã_{ij}Ã^{ij} = 0
 *
 * R = R̃ + R_φ, where:
 *   R̃ = Ricci scalar of γ̃_{ij}
 *   R_φ = -8e^{-4φ} D̃² e^{φ} (or equivalent conformal expression)
 *
 * Complexity: O(1) per point.
 */
double nr_bssn_hamiltonian_constraint(const nr_bssn_state_t *state,
                                       int i, int j, int k, double dx);

/**
 * Compute the momentum constraint in BSSN variables.
 * M^i = D̃_j Ã^{ij} + 6Ã^{ij} D̃_j φ - (2/3) D̃^i K = 0
 *
 * Returns 3-vector M^i.
 *
 * Complexity: O(1) per point.
 */
void nr_bssn_momentum_constraint(const nr_bssn_state_t *state,
                                  int i, int j, int k, double dx,
                                  double Mi[3]);

/**
 * Compute the Gamma constraint at a point.
 *
 * G^i = Γ̃^i + ∂_j γ̃^{ij}
 *
 * Complexity: O(1) per point.
 */
void nr_bssn_gamma_constraint(const nr_bssn_state_t *state,
                               int i, int j, int k, double dx,
                               double Gi[3]);

/**
 * Compute L2 norms of all three constraints over the interior grid.
 *
 * Complexity: O(N^3).
 */
void nr_bssn_constraint_norms(const nr_bssn_state_t *state, double dx,
                               double *norm_H, double *norm_M,
                               double *norm_G);

/* ===========================================================================
 * L5: Kreiss-Oliger Dissipation
 *
 * High-frequency numerical noise is suppressed by adding a dissipation term:
 *   ∂_t u → ∂_t u + ε (-1)^{r/2} h^r D_+^r D_-^r u
 *
 * where r is the dissipation order (typically 6 for 4th-order FD).
 * =========================================================================== */

/**
 * Apply Kreiss-Oliger dissipation to all BSSN evolved variables.
 *
 * @param state      Current state (modified in place)
 * @param eps        Dissipation strength (typ. 0.1-0.5)
 * @param order      Dissipation order (typ. 6 for 4th-order FD)
 * @param dx         Grid spacing
 *
 * Complexity: O(N^3).
 * Reference: Kreiss & Oliger (1973); Baumgarte & Shapiro (2010), Sec. 6.4
 */
void nr_bssn_kreiss_oliger(nr_bssn_state_t *state, double eps,
                            int order, double dx);

#ifdef __cplusplus
}
#endif

#endif /* NR_BSSN_H */
