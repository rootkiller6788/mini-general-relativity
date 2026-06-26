/**
 * @file nr_horizon.h
 * @brief Apparent horizon finder for numerical relativity
 *
 * An apparent horizon is the outermost marginally outer trapped surface
 * (MOTS) on a spatial hypersurface Σ_t. It is defined by the condition
 * that the expansion Θ of outgoing null geodesics vanishes:
 *
 *   Θ = D_i s^i + K_{ij} s^i s^j − K = 0
 *
 * where s^i is the outward-pointing unit normal to the 2-surface S
 * embedded in Σ_t.
 *
 * For spherical topology, the surface is parameterized as:
 *   h(θ, φ) = r_0 + Σ_{lm} a_{lm} Y_{lm}(θ, φ)
 *
 * References:
 *   - Gundlach, PRD 57, 3480 (1998) — fast flow method
 *   - Thornburg, Living Rev. Relativity 10, 3 (2007) — review
 *   - Baumgarte & Shapiro (2010), Ch. 7
 *   - Alcubierre (2008), Ch. 5
 *
 * Knowledge Coverage:
 *   L1: Apparent horizon, MOTS, expansion Θ
 *   L2: Horizon finding as root-finding on S^2
 *   L3: Spectral decomposition in spherical harmonics
 *   L5: Newton's method on S^2, flow method
 *   L6: Schwarzschild/Kerr apparent horizon
 */

#ifndef NR_HORIZON_H
#define NR_HORIZON_H

#include "nr_grid.h"
#include "nr_adm.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum spherical harmonic order for horizon representation */
#define NR_HORIZON_LMAX 8

/**
 * L1: Horizon surface represented in spherical harmonics.
 *
 * The surface radius as a function of angle:
 *   h(θ, φ) = Σ_{l=0}^{Lmax} Σ_{m=-l}^{l} a_{lm} Y_{lm}(θ, φ)
 *
 * For axisymmetric configurations, only m=0 terms survive.
 * h(θ) = Σ_{l=0}^{Lmax} a_l P_l(cos θ)
 */
typedef struct {
    int    lmax;                      /**< Maximum l for spherical harmonic expansion */
    double coeffs[(NR_HORIZON_LMAX+1)*(NR_HORIZON_LMAX+1)];  /**< a_{lm} coefficients, row-major */
    double centroid[3];              /**< (x0, y0, z0) center of the horizon */
    bool   found;                    /**< Whether a horizon was successfully located */
    double area;                     /**< Proper area of the horizon */
    double mass_irr;                 /**< Irreducible mass: M_irr = sqrt(A/(16π)) */
    double expansion_rms;            /**< RMS of Θ on the surface (convergence measure) */
} nr_horizon_t;

/* ===========================================================================
 * L2: Expansion Computation
 *
 * The expansion Θ of a null congruence orthogonal to a 2-surface:
 *   Θ = D_i s^i + s^i s^j K_{ij} − K
 *
 * For a surface r = h(θ, φ), the unit normal in Σ_t is:
 *   s^i = (1/|∇h|) γ^{ij} ∂_j (r − h(θ, φ))
 * =========================================================================== */

/**
 * Compute the expansion Θ at a specific angular direction (θ, φ)
 * on a trial surface radius h.
 *
 * @param gamma     3-metric
 * @param K         Extrinsic curvature
 * @param grid      Grid structure
 * @param dx        Grid spacing
 * @param h         Trial horizon radius in direction (θ, φ)
 * @param theta     Polar angle θ
 * @param phi       Azimuthal angle φ
 * @param centroid  (x0, y0, z0) center for the horizon
 * @return          Expansion Θ (should be 0 on apparent horizon)
 *
 * Complexity: O(1) — interpolation at one angular direction.
 * Reference: Baumgarte & Shapiro (2010), Eq. (7.3.5)
 */
double nr_horizon_expansion(const nr_sym_tensor3_t *gamma,
                             const nr_sym_tensor3_t *K,
                             const nr_grid_t *grid, double dx,
                             double h, double theta, double phi,
                             const double centroid[3]);

/**
 * Evaluate the horizon surface radius h(θ, φ) from spherical harmonic coefficients.
 *
 * @param horizon   Horizon structure (reads coeffs, centroid)
 * @param theta     Polar angle
 * @param phi       Azimuthal angle
 * @return          Surface radius at (θ, φ)
 *
 * Complexity: O(Lmax²).
 */
double nr_horizon_evaluate(const nr_horizon_t *horizon,
                            double theta, double phi);

/**
 * Compute the RMS expansion over the trial horizon surface.
 *
 * Evaluates Θ at (N_θ × N_φ) uniformly spaced points on the sphere
 * and returns sqrt(mean(Θ²)).
 *
 * @param gamma     3-metric
 * @param K         Extrinsic curvature
 * @param grid      Grid structure
 * @param dx        Grid spacing
 * @param horizon   Trial horizon (reads coeffs, centroid)
 * @param n_theta   Number of polar samples
 * @param n_phi     Number of azimuthal samples
 * @return          RMS expansion
 *
 * Complexity: O(N_θ × N_φ).
 */
double nr_horizon_expansion_rms(const nr_sym_tensor3_t *gamma,
                                 const nr_sym_tensor3_t *K,
                                 const nr_grid_t *grid, double dx,
                                 const nr_horizon_t *horizon,
                                 int n_theta, int n_phi);

/* ===========================================================================
 * L5: Apparent Horizon Finder — Flow Method
 *
 * The flow method (Gundlach 1998) evolves a trial surface toward the
 * apparent horizon by solving the parabolic equation:
 *
 *   ∂_λ h = −ρ Θ
 *
 * where λ is a fictitious "flow time" and ρ is a positive weight function.
 * As λ → ∞, h → h_AH (the apparent horizon radius), where Θ → 0.
 *
 * =========================================================================== */

/**
 * Find apparent horizon using the fast flow method.
 *
 * Starts from an initial spherical guess at radius r0 and iteratively
 * deforms the surface until Θ ≈ 0 everywhere.
 *
 * @param gamma     3-metric
 * @param K         Extrinsic curvature
 * @param grid      Grid structure
 * @param dx        Grid spacing
 * @param r0        Initial guess radius
 * @param centroid  Center coordinates
 * @param lmax      Maximum spherical harmonic order
 * @param max_iter  Maximum flow iterations
 * @param tol       Convergence tolerance (RMS expansion)
 * @param dt_flow   Flow step size
 * @param horizon   Output horizon (filled on success)
 * @return          Number of iterations, or -1 on failure
 *
 * Reference: Gundlach, PRD 57, 3480 (1998)
 * Complexity: O(N_iter × N_angular × Lmax²).
 */
int nr_horizon_find_flow(const nr_sym_tensor3_t *gamma,
                          const nr_sym_tensor3_t *K,
                          const nr_grid_t *grid, double dx,
                          double r0, const double centroid[3],
                          int lmax, int max_iter, double tol,
                          double dt_flow, nr_horizon_t *horizon);

/* ===========================================================================
 * L6: Known Horizon Locations
 * =========================================================================== */

/**
 * Get the exact apparent horizon radius for Schwarzschild.
 *
 * In isotropic coordinates, the apparent horizon is at r = M/2.
 * In Schwarzschild areal coordinates, it's at R = 2M.
 *
 * @param M           ADM mass
 * @param isotropic   true → return in isotropic coords, false → areal coords
 * @return            Horizon radius
 *
 * Complexity: O(1).
 */
double nr_horizon_schwarzschild_radius(double M, bool isotropic);

/**
 * Get the exact apparent horizon radius for Kerr.
 *
 * In Boyer-Lindquist coordinates: r_+ = M + sqrt(M² − a²)
 *
 * @param M    Mass
 * @param a    Spin parameter (|a| ≤ M)
 * @return     Outer horizon radius in Boyer-Lindquist coordinates
 *
 * Complexity: O(1).
 */
double nr_horizon_kerr_radius(double M, double a);

/**
 * Find apparent horizon for Schwarzschild as a verification test.
 *
 * Given exact Schwarzschild data on the grid, call the flow finder
 * and verify that the found horizon matches the analytic location.
 *
 * @param M          Mass
 * @param grid       Grid structure
 * @param dx         Grid spacing
 * @param horizon    Output horizon
 * @return           0 on success, -1 on failure
 *
 * Complexity: O(N_iter × N_angular).
 */
int nr_horizon_find_schwarzschild(double M, const nr_grid_t *grid,
                                   double dx, nr_horizon_t *horizon);

/**
 * Compute the proper area of the horizon surface.
 *
 * A = ∮_S sqrt(det(h_AB)) dθ dφ
 *
 * where h_AB is the induced 2-metric on the surface S.
 *
 * @param gamma     3-metric
 * @param grid      Grid structure
 * @param dx        Grid spacing
 * @param horizon   Horizon definition (coeffs and centroid)
 * @param n_theta, n_phi  Angular integration resolution
 * @return          Proper area of the surface
 *
 * Complexity: O(N_θ × N_φ).
 */
double nr_horizon_area(const nr_sym_tensor3_t *gamma,
                        const nr_grid_t *grid, double dx,
                        const nr_horizon_t *horizon,
                        int n_theta, int n_phi);

#ifdef __cplusplus
}
#endif

#endif /* NR_HORIZON_H */
