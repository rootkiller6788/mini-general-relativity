/**
 * @file nr_adm.h
 * @brief ADM (Arnowitt-Deser-Misner) 3+1 decomposition
 *
 * The ADM formalism foliates spacetime into spatial hypersurfaces Σ_t
 * parameterized by a time coordinate t. This enables Einstein's equations
 * to be written as a Cauchy initial value problem.
 *
 * Spacetime metric in 3+1 form:
 *   ds² = -α² dt² + γ_{ij} (dx^i + β^i dt)(dx^j + β^j dt)
 *
 * where:
 *   α     = lapse function (proper time between hypersurfaces)
 *   β^i   = shift vector (spatial coordinate threading)
 *   γ_{ij} = 3-metric on spatial hypersurface Σ_t
 *
 * References:
 *   - Arnowitt, Deser, Misner, "Dynamics of General Relativity" (1962)
 *   - Wald, "General Relativity" (1984), Appendix E
 *   - Baumgarte & Shapiro (2010), Ch. 2
 *   - Gourgoulhon, "3+1 Formalism in General Relativity" (2012)
 *
 * Knowledge Coverage:
 *   L1: lapse α, shift β^i, 3-metric γ_{ij}, extrinsic curvature K_{ij}
 *   L2: Foliation of spacetime, Cauchy problem
 *   L3: Riemannian geometry on Σ_t, covariant derivative D_i
 *   L4: ADM evolution equations, Hamiltonian & momentum constraints
 *   L5: Finite-difference implementation of ADM RHS
 */

#ifndef NR_ADM_H
#define NR_ADM_H

#include "nr_grid.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
 * L1: ADM State Variables
 *
 * The complete state of the gravitational field on a spatial hypersurface
 * is described by the 3-metric γ_{ij} and the extrinsic curvature K_{ij}.
 * Together they form the canonical pair (γ_{ij}, K_{ij}) in the Hamiltonian
 * formulation of GR.
 * =========================================================================== */

/**
 * ADM dynamical variables on a 3D grid.
 *
 * γ_{ij}: 3-metric — symmetric, positive-definite, 6 components
 * K_{ij}: extrinsic curvature — symmetric, 6 components
 *
 * The gauge variables (lapse α and shift β^i) are freely specifiable
 * and determine the foliation and spatial coordinates respectively.
 */
typedef struct {
    nr_sym_tensor3_t *gamma;  /**< 3-metric γ_{ij} (6 symmetric components) */
    nr_sym_tensor3_t *K;      /**< Extrinsic curvature K_{ij} */
    nr_gf_t          *alpha;  /**< Lapse function α */
    nr_vector3_t     *beta;   /**< Shift vector β^i (3 components) */
} nr_adm_state_t;

/* ===========================================================================
 * L1: ADM State Allocation and Management
 * =========================================================================== */

/**
 * Allocate a full ADM state on a 3D grid.
 * Includes γ_{ij} (6 comp), K_{ij} (6 comp), α (1 comp), β^i (3 comp).
 * Total: 16 grid functions.
 * Complexity: O(N^3) memory.
 */
nr_adm_state_t* nr_adm_alloc(int nx, int ny, int nz, int ng);

/** Free all memory associated with an ADM state */
void nr_adm_free(nr_adm_state_t *state);

/**
 * Set Minkowski initial data: γ_{ij} = δ_{ij}, K_{ij} = 0, α = 1, β^i = 0.
 * Complexity: O(N^3).
 */
void nr_adm_set_minkowski(nr_adm_state_t *state);

/* ===========================================================================
 * L3: 3-Covariant Derivative
 *
 * On a spatial slice Σ_t, the covariant derivative D_i is compatible
 * with the 3-metric: D_k γ_{ij} = 0.
 *
 * Christoffel symbols of the 3-metric:
 *   Γ^i_{jk} = (1/2) γ^{il} (∂_j γ_{lk} + ∂_k γ_{lj} - ∂_l γ_{jk})
 * =========================================================================== */

/**
 * Compute 3-Christoffel symbol Γ^i_{jk} at a grid point.
 *
 * The Christoffel symbols are computed from finite-difference derivatives
 * of the 3-metric γ_{ij} and the inverse 3-metric γ^{ij}.
 *
 * @param gamma     3-metric tensor
 * @param i, j, k   Grid indices
 * @param dx        Grid spacing
 * @param Gamma_i_jk[3][3][3]  Output: Γ^i_{jk} (27 scalars, note symmetry Γ^i_{jk} = Γ^i_{kj})
 *
 * Complexity: O(1) per call, used in loops over grid.
 * Reference: Wald (1984), Eq. (3.1.30)
 */
void nr_adm_christoffel(const nr_sym_tensor3_t *gamma,
                         int i, int j, int k, double dx,
                         double Gamma[3][3][3]);

/**
 * Compute the 3-Ricci tensor R_{ij} at a grid point.
 *
 * R_{ij} = ∂_k Γ^k_{ij} - ∂_i Γ^k_{kj} + Γ^k_{kl} Γ^l_{ij} - Γ^k_{il} Γ^l_{kj}
 *
 * @param gamma     3-metric tensor
 * @param i, j, k   Grid indices
 * @param dx        Grid spacing
 * @param Rij[3][3] Output: 3-Ricci tensor (symmetric)
 *
 * Complexity: O(1).
 * Reference: Wald (1984), Eq. (3.4.3)
 */
void nr_adm_ricci(const nr_sym_tensor3_t *gamma,
                  int i, int j, int k, double dx,
                  double Rij[3][3]);

/**
 * Compute the 3-Ricci scalar R = γ^{ij} R_{ij} at a grid point.
 *
 * R = γ^{xx}R_{xx} + γ^{yy}R_{yy} + γ^{zz}R_{zz}
 *   + 2(γ^{xy}R_{xy} + γ^{xz}R_{xz} + γ^{yz}R_{yz})
 *
 * @param Rij   3-Ricci tensor (from nr_adm_ricci)
 * @param inv_gamma  Inverse 3-metric γ^{ij}
 * @return       Ricci scalar R
 *
 * Complexity: O(1).
 * Reference: Wald (1984), Eq. (3.4.5)
 */
double nr_adm_ricci_scalar(const double Rij[3][3],
                           const double inv_gamma[3][3]);

/**
 * Compute the inverse 3-metric γ^{ij} from γ_{ij} at a point.
 *
 * For a symmetric 2x2 (upper-left, for debug) or 3x3 matrix:
 * Uses analytic 3x3 symmetric matrix inversion.
 *
 * @param gamma_ij[6]  3-metric components {xx, xy, xz, yy, yz, zz} (in)
 * @param inv_gamma[3][3] Inverse metric (out)
 * @return Determinant of γ_{ij}
 *
 * Complexity: O(1).
 */
double nr_adm_inverse_metric(const double gamma_ij[6],
                             double inv_gamma[3][3]);

/* ===========================================================================
 * L4: Einstein Constraints
 *
 * On each spatial hypersurface, the initial data must satisfy the
 * Hamiltonian and momentum constraints. These are elliptic equations
 * that restrict the allowed (γ_{ij}, K_{ij}) pairs.
 * =========================================================================== */

/**
 * L4: Hamiltonian constraint (scalar equation).
 *
 * H = R + K² - K_{ij} K^{ij} = 0   (vacuum)
 *
 * where R is the 3-Ricci scalar, K = tr(K_{ij}) = γ^{ij} K_{ij}.
 *
 * @param gamma  3-metric
 * @param K       Extrinsic curvature
 * @param i,j,k  Grid point
 * @param dx     Grid spacing
 * @return       Value of Hamiltonian constraint H (should be ~0)
 *
 * Reference: Wald (1984), Eq. (E.2.32)
 * Complexity: O(1).
 */
double nr_adm_hamiltonian_constraint(const nr_sym_tensor3_t *gamma,
                                      const nr_sym_tensor3_t *K,
                                      int i, int j, int k, double dx);

/**
 * L4: Momentum constraint (vector equation, 3 components).
 *
 * M_i = D_j K^j_i - D_i K = 0   (vacuum)
 *
 * Returns the 3-vector M_i of the momentum constraint.
 *
 * @param gamma  3-metric
 * @param K       Extrinsic curvature
 * @param i,j,k  Grid point
 * @param dx     Grid spacing
 * @param Mi[3]  Output: momentum constraint M_i
 *
 * Reference: Wald (1984), Eq. (E.2.33)
 * Complexity: O(1).
 */
void nr_adm_momentum_constraint(const nr_sym_tensor3_t *gamma,
                                 const nr_sym_tensor3_t *K,
                                 int i, int j, int k, double dx,
                                 double Mi[3]);

/**
 * Compute the L2 norm of both constraints over the interior grid.
 * Returns norms via pointers.
 *
 * Complexity: O(N^3).
 */
void nr_adm_constraint_norms(const nr_adm_state_t *state, double dx,
                              double *norm_H, double *norm_M);

/* ===========================================================================
 * L4: ADM Evolution Equations (Right-Hand Side)
 *
 * The ADM evolution equations advance γ_{ij} and K_{ij} in time.
 *
 * ∂_t γ_{ij} = -2α K_{ij} + D_i β_j + D_j β_i
 * ∂_t K_{ij} = -D_i D_j α + α(R_{ij} + K K_{ij} - 2K_{ik} K^k_j)
 *               + β^k D_k K_{ij} + K_{ik} D_j β^k + K_{kj} D_i β^k
 *
 * Reference: Baumgarte & Shapiro (2010), Eqs. (2.5.5)-(2.5.6)
 * =========================================================================== */

/**
 * Compute RHS of the ADM evolution equations at a single grid point.
 *
 * ∂_t γ_{ij} — stored in dt_gamma[6]
 * ∂_t K_{ij} — stored in dt_K[6]
 *
 * @param state       Current ADM state
 * @param i, j, k     Grid point indices
 * @param dx          Grid spacing
 * @param dt_gamma[6] Output: ∂_t γ_{ij}
 * @param dt_K[6]     Output: ∂_t K_{ij}
 *
 * Complexity: O(1).
 */
void nr_adm_rhs_point(const nr_adm_state_t *state,
                       int i, int j, int k, double dx,
                       double dt_gamma[6], double dt_K[6]);

/**
 * Compute the RHS of ADM evolution over the entire interior grid.
 * Fills dt_state with time derivatives.
 *
 * Complexity: O(N^3).
 */
void nr_adm_rhs(const nr_adm_state_t *state, double dx,
                nr_adm_state_t *dt_state);

/* ===========================================================================
 * L2: ADM Mass (Energy at spatial infinity)
 *
 * The ADM mass measures the total energy contained in an asymptotically
 * flat spacelike hypersurface.
 * =========================================================================== */

/**
 * Compute ADM mass via the surface integral at infinity.
 *
 * For a conformally flat metric γ_{ij} = ψ^4 δ_{ij}:
 *   M_ADM = -(1/(2π)) ∮_∞ ∂_i ψ dS^i
 *
 * For Schwarzschild in isotropic coordinates (ψ = 1 + M/(2r)):
 *   Should return M.
 *
 * @param psi    Conformal factor ψ on the grid
 * @param grid   Grid structure
 * @param r_ext  Extraction radius (should be near outer boundary)
 * @return       ADM mass estimate
 *
 * Reference: Wald (1984), Eq. (11.2.5); Baumgarte & Shapiro (2010), Eq. (3.3.100)
 * Complexity: O(N^2) — surface integral at fixed r.
 */
double nr_adm_mass(const nr_gf_t *psi, const nr_grid_t *grid, double r_ext);

/**
 * Compute the Komar mass for a stationary spacetime.
 *
 * M_Komar = (1/(4π)) ∮_S (D^i α) dS_i
 *
 * For Schwarzschild in standard coordinates: returns M.
 * For Kerr: returns M (not M_irr).
 *
 * @param alpha  Lapse function
 * @param gamma  3-metric (for proper area element)
 * @param grid   Grid structure
 * @param r_ext  Extraction surface radius
 * @return       Komar mass
 *
 * Reference: Wald (1984), Eq. (11.2.9-10)
 * Complexity: O(N^2).
 */
double nr_adm_komar_mass(const nr_gf_t *alpha,
                          const nr_sym_tensor3_t *gamma,
                          const nr_grid_t *grid, double r_ext);

#ifdef __cplusplus
}
#endif

#endif /* NR_ADM_H */
