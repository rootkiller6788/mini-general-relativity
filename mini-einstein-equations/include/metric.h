/**
 * @file metric.h
 * @brief Metric tensor, Christoffel symbols, covariant derivative.
 *
 * Reference: Wald (1984), Ch.3; Carroll (2004), Ch.3
 *            MIT 8.962 Lecture Notes, Schutz "A First Course in GR"
 *
 * Knowledge map:
 *   L1: Metric g_{mu,nu}, inverse metric g^{mu,nu}
 *   L3: Levi-Civita connection (Christoffel symbols)
 *   L3: Covariant derivative
 *   L3: Metric compatibility ∇_rho g_{mu,nu} = 0
 */

#ifndef MINI_EINSTEIN_METRIC_H
#define MINI_EINSTEIN_METRIC_H

#include "tensor.h"

/* ---------------------------------------------------------------------------
 * L1: Metric structure
 * -------------------------------------------------------------------------*/

/**
 * Spacetime metric on 4D Lorentzian manifold.
 * Contains both covariant (lower-index) and contravariant (upper-index) forms.
 */
typedef struct {
    Tensor2 g;        /**< g_{mu,nu} — covariant metric */
    Tensor2 g_inv;    /**< g^{mu,nu} — inverse/contravariant metric */
    double det_g;     /**< determinant of g_{mu,nu} (< 0 for Lorentzian) */
    int is_valid;     /**< 1 if metric is non-degenerate, 0 otherwise */
} Metric;

/* ---------------------------------------------------------------------------
 * L3: Metric operations
 * -------------------------------------------------------------------------*/

/**
 * Initialize metric to Minkowski (flat spacetime).
 * g_{mu,nu} = diag(-1, 1, 1, 1)
 */
void metric_init_minkowski(Metric *m);

/**
 * Compute matrix inverse of 4x4 metric using Gauss-Jordan elimination.
 * Also computes determinant.
 * @return 0 on success, -1 if singular (|det| < 1e-12)
 */
int metric_compute_inverse(Metric *m);

/**
 * Compute determinant of 4x4 matrix via Laplace expansion.
 * @param mat  4x4 matrix
 * @return determinant
 */
double matrix4x4_determinant(const Tensor2 mat);

/**
 * Compute cofactor matrix (adjugate / det). Used internally.
 */
void matrix4x4_cofactor(const Tensor2 mat, Tensor2 cof);

/**
 * Check signature of metric (count positive/negative eigenvalues).
 * For a physical Lorentzian metric: 1 negative + 3 positive eigenvalues.
 * @return 1 if Lorentzian, 0 otherwise
 */
int metric_check_lorentz_signature(const Metric *m);

/* ---------------------------------------------------------------------------
 * L3: Christoffel symbols (Levi-Civita connection)
 * -------------------------------------------------------------------------*/

/**
 * Compute Christoffel symbols Gamma^{lambda}_{mu,nu} from metric and its
 * first partial derivatives:
 *
 *   Gamma^{lambda}_{mu,nu} = (1/2) * g^{lambda,sigma} *
 *       (d_g_{sigma,mu}/dx^{nu} + d_g_{sigma,nu}/dx^{mu} - d_g_{mu,nu}/dx^{sigma})
 *
 * @param m         metric with valid g_inv
 * @param dg        partial derivatives: dg[sigma][mu][nu] = ∂_{sigma} g_{mu,nu}
 *                  Note: first index is the derivative direction.
 * @param Gamma     output: Gamma^{lambda}_{mu,nu}
 */
void metric_christoffel(const Metric *m, const Tensor3 dg, Tensor3 Gamma);

/**
 * Compute Christoffel symbols in terms of their symmetry.
 * On a torsion-free connection: Gamma^{lambda}_{mu,nu} = Gamma^{lambda}_{nu,mu}
 * This version verifies the symmetry of the result.
 */
void metric_christoffel_symmetric(const Metric *m, const Tensor3 dg, Tensor3 Gamma);

/* ---------------------------------------------------------------------------
 * L3: Covariant derivative
 * -------------------------------------------------------------------------*/

/**
 * Covariant derivative of rank-1 contravariant vector:
 *   ∇_{mu} V^{nu} = ∂_{mu} V^{nu} + Gamma^{nu}_{mu,sigma} V^{sigma}
 *
 * @param Gamma     Christoffel symbols
 * @param V         vector field V^{nu}
 * @param dV        partial derivatives: dV[mu][nu] = ∂_{mu} V^{nu}
 * @param covV      output: ∇_{mu} V^{nu}
 */
void covariant_derivative_vector(const Tensor3 Gamma,
                                  const Vector4 V,
                                  const Tensor2 dV,
                                  Tensor2 covV);

/**
 * Covariant derivative of rank-1 covariant vector (1-form):
 *   ∇_{mu} omega_{nu} = ∂_{mu} omega_{nu} - Gamma^{sigma}_{mu,nu} omega_{sigma}
 *
 * @param Gamma     Christoffel symbols
 * @param omega     1-form omega_{nu}
 * @param d_omega   partial derivatives: d_omega[mu][nu] = ∂_{mu} omega_{nu}
 * @param cov_omega output: ∇_{mu} omega_{nu}
 */
void covariant_derivative_oneform(const Tensor3 Gamma,
                                   const Vector4 omega,
                                   const Tensor2 d_omega,
                                   Tensor2 cov_omega);

/**
 * Covariant derivative of rank-2 covariant tensor:
 *   ∇_{rho} T_{mu,nu} = ∂_{rho} T_{mu,nu}
 *     - Gamma^{sigma}_{rho,mu} T_{sigma,nu}
 *     - Gamma^{sigma}_{rho,nu} T_{mu,sigma}
 *
 * @param Gamma     Christoffel symbols
 * @param T         rank-2 tensor T_{mu,nu}
 * @param dT        partial derivatives: dT[rho][mu][nu] = ∂_{rho} T_{mu,nu}
 * @param covT      output: ∇_{rho} T_{mu,nu} — rank-3 tensor
 */
void covariant_derivative_tensor2(const Tensor3 Gamma,
                                   const Tensor2 T,
                                   const Tensor3 dT,
                                   Tensor3 covT);

/**
 * Verify metric compatibility: ∇_{rho} g_{mu,nu} = 0
 * This is the defining property of the Levi-Civita connection.
 * @return max absolute component of ∇ g (should be ≈ 0)
 */
double metric_compatibility_check(const Metric *m, const Tensor3 Gamma);

/**
 * Compute d'Alembertian (wave operator) on scalar field:
 *   □ phi = g^{mu,nu} * ∇_{mu} ∇_{nu} phi
 *         = g^{mu,nu} * (∂_{mu} ∂_{nu} phi - Gamma^{sigma}_{mu,nu} ∂_{sigma} phi)
 *
 * @param m        metric
 * @param Gamma    Christoffel symbols
 * @param dphi     first derivatives ∂_{mu} phi
 * @param ddphi    second derivatives ∂_{mu} ∂_{nu} phi
 * @return □ phi
 */
double dalembertian_scalar(const Metric *m, const Tensor3 Gamma,
                            const Vector4 dphi, const Tensor2 ddphi);

#endif /* MINI_EINSTEIN_METRIC_H */
