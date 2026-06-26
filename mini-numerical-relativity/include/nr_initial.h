/**
 * @file nr_initial.h
 * @brief Initial data construction for numerical relativity
 *
 * The Einstein constraint equations (Hamiltonian + momentum) must be
 * satisfied on the initial hypersurface. We solve the conformal
 * transverse-traceless (CTT) decomposition:
 *
 *   γ_{ij} = ψ^4 γ̄_{ij}          (conformal metric)
 *   K_{ij} = ψ^{-2} Â_{ij} + (1/3) γ_{ij} K   (conformal curvature)
 *
 * The Hamiltonian constraint becomes an elliptic equation for ψ:
 *   8 Δ̄ ψ − R̄ ψ + (2/3) K^2 ψ^5 − Â_{ij} Â^{ij} ψ^{-7} = 0
 *
 * The momentum constraint becomes an elliptic equation for the vector
 * potential that generates Â^{ij}_{TT}.
 *
 * References:
 *   - Bowen & York, PRD 21, 2047 (1980)
 *   - Brill & Lindquist, PR 131, 471 (1963)
 *   - Brandt & Brügmann, PRL 78, 3606 (1997) — puncture method
 *   - Baumgarte & Shapiro (2010), Ch. 3 & 12
 *
 * Knowledge Coverage:
 *   L1: Conformal factor ψ, puncture mass, Bowen-York momentum
 *   L2: Conformal transverse-traceless decomposition
 *   L3: Elliptic equations on conformal background
 *   L4: Hamiltonian & momentum constraints in CTT form
 *   L5: Puncture method, Bowen-York solutions
 *   L6: Schwarzschild, Kerr, Brill-Lindquist, Bowen-York initial data
 */

#ifndef NR_INITIAL_H
#define NR_INITIAL_H

#include "nr_grid.h"
#include "nr_adm.h"
#include "nr_bssn.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
 * L1: Puncture Data
 *
 * The puncture method (Brandt & Brügmann 1997) writes:
 *   ψ = 1 + Σ_A m_A/(2r_A) + u
 *
 * where the first two terms are the Brill-Lindquist solution and u is a
 * regular correction solved from an elliptic equation.
 *
 * r_A = |x - x_A| is the coordinate distance from puncture A.
 * =========================================================================== */

/**
 * L6: Brill-Lindquist initial data — multiple Schwarzschild black holes
 * at the moment of time symmetry (K_{ij} = 0).
 *
 * Conformal factor: ψ = 1 + Σ_A m_A/(2r_A)
 * This solves the Hamiltonian constraint exactly.
 *
 * Each puncture has mass parameter m_A and position x_A.
 *
 * @param state    ADM state to fill
 * @param grid     Grid definition
 * @param n_punc   Number of punctures (≥ 1)
 * @param masses   Puncture mass parameters m_A (length n_punc)
 * @param pos_x    Puncture x-positions (length n_punc)
 * @param pos_y    Puncture y-positions (length n_punc)
 * @param pos_z    Puncture z-positions (length n_punc)
 *
 * Reference: Brill & Lindquist, PR 131, 471 (1963)
 * Complexity: O(N^3).
 */
void nr_init_brill_lindquist(nr_adm_state_t *state, const nr_grid_t *grid,
                              int n_punc, const double *masses,
                              const double *pos_x, const double *pos_y,
                              const double *pos_z);

/**
 * L6: Bowen-York initial data — black holes with linear momentum and spin.
 *
 * The extrinsic curvature is constructed from the Bowen-York solution:
 *   Â^{ij} = (3/(2r^2)) [P^i n^j + P^j n^i − (δ^{ij} − n^i n^j) P_k n^k]
 *            + (3/r^3) [ε^{ikl} S_l n_k n^j + ε^{jkl} S_l n_k n^i]
 *
 * where P^i is the linear momentum, S^i is the spin angular momentum,
 * and n^i is the unit normal from the puncture.
 *
 * Then the Hamiltonian constraint is solved for ψ (using Brill-Lindquist
 * as initial guess for iterative elliptic solver).
 *
 * @param state    ADM state to fill
 * @param grid     Grid definition
 * @param n_punc   Number of punctures
 * @param masses   Mass parameters
 * @param pos_x, pos_y, pos_z  Positions
 * @param Px, Py, Pz           Linear momenta
 * @param Sx, Sy, Sz           Spin angular momenta
 *
 * Reference: Bowen & York, PRD 21, 2047 (1980)
 * Complexity: O(N^3).
 */
void nr_init_bowen_york(nr_adm_state_t *state, const nr_grid_t *grid,
                         int n_punc, const double *masses,
                         const double *pos_x, const double *pos_y,
                         const double *pos_z,
                         const double *Px, const double *Py, const double *Pz,
                         const double *Sx, const double *Sy, const double *Sz);

/* ===========================================================================
 * L6: Schwarzschild Initial Data
 *
 * Schwarzschild black hole in various coordinate systems.
 * =========================================================================== */

/**
 * Schwarzschild in isotropic coordinates.
 *
 * γ_{ij} = ψ^4 δ_{ij}
 * ψ = 1 + M/(2r)
 * K_{ij} = 0 (moment of time symmetry)
 *
 * ADM mass = M.
 *
 * @param state    ADM state to fill
 * @param grid     Grid definition
 * @param M        ADM mass of the black hole
 *
 * Reference: Baumgarte & Shapiro (2010), Sec. 3.3
 * Complexity: O(N^3).
 */
void nr_init_schwarzschild_isotropic(nr_adm_state_t *state,
                                      const nr_grid_t *grid, double M);

/**
 * Schwarzschild in trumpet (1+log) coordinates.
 *
 * The maximal slicing of Schwarzschild:
 *   α = sqrt(1 - 2M/R)  (where R is areal radius)
 *   R = r (1 + M/(2r))^2  (relation to isotropic radius)
 *
 * This is the stationary end-state of 1+log sliced evolutions.
 *
 * @param state    BSSN state to fill
 * @param grid     Grid definition  
 * @param M        ADM mass
 *
 * Reference: Hannam et al., PRL 99, 241102 (2007); Baumgarte & Shapiro (2010), Sec. 8.4
 * Complexity: O(N^3).
 */
void nr_init_schwarzschild_trumpet(nr_bssn_state_t *state,
                                    const nr_grid_t *grid, double M);

/* ===========================================================================
 * L6: Kerr-Schild Initial Data
 *
 * Kerr black hole in Kerr-Schild coordinates:
 *   ds^2 = (η_{μν} + 2H l_μ l_ν) dx^μ dx^ν
 *
 * where H = M r / (r^2 + a^2 cos^2 θ) and l_μ is a null vector.
 *
 * These coordinates are horizon-penetrating (no coordinate singularity).
 * =========================================================================== */

/**
 * Kerr-Schild initial data for a spinning black hole.
 *
 * @param state    ADM state to fill
 * @param grid     Grid definition
 * @param M        Mass
 * @param a        Spin parameter (a = J/M, |a| ≤ M)
 *
 * Reference: Kerr, PRL 11, 237 (1963); Kerr & Schild (1965)
 * Complexity: O(N^3).
 */
void nr_init_kerr_schild(nr_adm_state_t *state, const nr_grid_t *grid,
                          double M, double a);

/**
 * Convert Kerr-Schild to BSSN variables.
 *
 * Kerr-Schild in BSSN form is used as initial data for spinning black
 * hole evolutions.
 *
 * @param state    BSSN state to fill
 * @param grid     Grid definition
 * @param M        Mass
 * @param a        Spin parameter
 *
 * Complexity: O(N^3).
 */
void nr_init_kerr_schild_bssn(nr_bssn_state_t *state, const nr_grid_t *grid,
                               double M, double a);

/* ===========================================================================
 * L2: Conformal Factor Solver
 *
 * Solve the Hamiltonian constraint for the conformal factor ψ:
 *   8 Δ̄ ψ − R̄ ψ + (2/3) K^2 ψ^5 − Â_{ij} Â^{ij} ψ^{-7} = 0
 *
 * with outer boundary condition ψ → 1 + M_ADM/(2r) as r → ∞.
 * =========================================================================== */

/**
 * Solve the Hamiltonian constraint via Gauss-Seidel relaxation.
 *
 * For the conformally flat case (γ̄_{ij} = δ_{ij}, R̄ = 0) with K = 0:
 *   Δ_flat ψ = -(1/8) Â_{ij} Â^{ij} ψ^{-7}
 *
 * @param psi      Initial guess for ψ (modified in place)
 * @param A2       |Â|^2 = Â_{ij} Â^{ij} on the grid
 * @param grid     Grid structure
 * @param dx       Grid spacing
 * @param tol      Convergence tolerance
 * @param max_iter Maximum iterations
 * @param omega    SOR over-relaxation parameter (1.0 = no over-relaxation)
 * @return         Number of iterations performed
 *
 * Reference: Baumgarte & Shapiro (2010), Sec. 3.2
 * Complexity: O(N^3 · iter).
 */
int nr_init_solve_psi(nr_gf_t *psi, const nr_gf_t *A2,
                       const nr_grid_t *grid, double dx,
                       double tol, int max_iter, double omega);

/**
 * Set up the puncture equation source for Hamiltonian constraint solver.
 *
 * For punctures, ψ = 1 + Σ m_A/(2r_A) + u, solve for u:
 *   Δ_flat u = -β ψ^{-7} (in interior, u → 0 at boundary)
 *
 * @param u        Regular correction (filled in place)
 * @param psi0     Brill-Lindquist initial ψ = 1 + Σ m_A/(2r_A) (not modified)
 * @param grid     Grid structure
 * @param dx       Grid spacing
 * @param tol      Convergence tolerance
 * @param max_iter Maximum iterations
 * @param omega    SOR parameter
 * @return         Number of iterations
 *
 * Reference: Brandt & Brügmann, PRL 78, 3606 (1997)
 * Complexity: O(N^3 · iter).
 */
int nr_init_solve_puncture_u(nr_gf_t *u, const nr_gf_t *psi0,
                              const nr_grid_t *grid, double dx,
                              double tol, int max_iter, double omega);

#ifdef __cplusplus
}
#endif

#endif /* NR_INITIAL_H */
