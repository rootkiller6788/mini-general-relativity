/**
 * @file curvature.h
 * @brief Riemann curvature tensor, Ricci tensor, Einstein tensor, Weyl tensor.
 *
 * Reference: Wald (1984), Ch.3; Carroll (2004), Ch.3
 *            MIT 8.962 General Relativity
 *
 * Knowledge map:
 *   L1: Riemann tensor R^{rho}_{sigma,mu,nu}
 *   L1: Ricci tensor R_{mu,nu}, Ricci scalar R
 *   L1: Einstein tensor G_{mu,nu} = R_{mu,nu} - (1/2)R g_{mu,nu}
 *   L1: Weyl tensor (conformal curvature)
 *   L3: Bianchi identities
 *   L4: Einstein field equations
 */

#ifndef MINI_EINSTEIN_CURVATURE_H
#define MINI_EINSTEIN_CURVATURE_H

#include "metric.h"

/* ---------------------------------------------------------------------------
 * L1: Curvature tensor structures
 * -------------------------------------------------------------------------*/

/** Full curvature information at a spacetime point */
typedef struct {
    Tensor4 Riemann;     /**< R^{rho}_{sigma,mu,nu} */
    Tensor2 Ricci;       /**< R_{mu,nu} */
    Tensor2 Einstein;    /**< G_{mu,nu} */
    double  RicciScalar; /**< R = g^{mu,nu} R_{mu,nu} */
    double  Kretschmann; /**< R_{mu,nu,rho,sigma} R^{mu,nu,rho,sigma} */
} Curvature;

/* ---------------------------------------------------------------------------
 * L1: Riemann tensor computation
 * -------------------------------------------------------------------------*/

/**
 * Compute Riemann curvature tensor from Christoffel symbols and their
 * first partial derivatives:
 *
 *   R^{rho}_{sigma,mu,nu} = ∂_{mu} Gamma^{rho}_{nu,sigma}
 *                         - ∂_{nu} Gamma^{rho}_{mu,sigma}
 *                         + Gamma^{rho}_{mu,lambda} * Gamma^{lambda}_{nu,sigma}
 *                         - Gamma^{rho}_{nu,lambda} * Gamma^{lambda}_{mu,sigma}
 *
 * @param Gamma     Christoffel symbols Gamma^{rho}_{mu,nu}
 * @param dGamma    partial derivatives: dGamma[lambda][rho][mu][nu]
 *                  = ∂_{lambda} Gamma^{rho}_{mu,nu}
 * @param R_riemann output Riemann tensor
 */
void compute_riemann_tensor(const Tensor3 Gamma,
                             const Tensor4 dGamma,
                             Tensor4 R_riemann);

/**
 * Compute fully covariant Riemann tensor:
 *   R_{rho,sigma,mu,nu} = g_{rho,lambda} * R^{lambda}_{sigma,mu,nu}
 */
void riemann_covariant(const Metric *m, const Tensor4 R_riemann_up,
                       Tensor4 R_riemann_dn);

/* ---------------------------------------------------------------------------
 * L1: Ricci tensor and scalar
 * -------------------------------------------------------------------------*/

/**
 * Compute Ricci tensor from Riemann:
 *   R_{mu,nu} = R^{rho}_{mu,rho,nu}
 */
void compute_ricci_tensor(const Tensor4 R_riemann, Tensor2 R_ricci);

/**
 * Compute Ricci scalar:
 *   R = g^{mu,nu} * R_{mu,nu}
 */
double compute_ricci_scalar(const Metric *m, const Tensor2 R_ricci);

/* ---------------------------------------------------------------------------
 * L1: Einstein tensor
 * -------------------------------------------------------------------------*/

/**
 * Compute Einstein tensor:
 *   G_{mu,nu} = R_{mu,nu} - (1/2) * R * g_{mu,nu}
 */
void compute_einstein_tensor(const Metric *m, const Tensor2 R_ricci,
                              double R_scalar, Tensor2 G);

/* ---------------------------------------------------------------------------
 * L1: Weyl tensor (conformal curvature)
 * -------------------------------------------------------------------------*/

/**
 * Compute Weyl tensor (the trace-free part of Riemann):
 *
 *   C_{rho,sigma,mu,nu} = R_{rho,sigma,mu,nu}
 *     - (1/2) * (g_{rho,mu} R_{sigma,nu} - g_{rho,nu} R_{sigma,mu}
 *                - g_{sigma,mu} R_{rho,nu} + g_{sigma,nu} R_{rho,mu})
 *     + (R/6) * (g_{rho,mu} g_{sigma,nu} - g_{rho,nu} g_{sigma,mu})
 *
 * The Weyl tensor vanishes identically in 3 or fewer dimensions.
 * In vacuum (R_{mu,nu}=0) the Weyl tensor equals the Riemann tensor.
 */
void compute_weyl_tensor(const Metric *m, const Tensor4 R_dn,
                          const Tensor2 R_ricci, double R_scalar,
                          Tensor4 Weyl);

/* ---------------------------------------------------------------------------
 * L3: Curvature invariants
 * -------------------------------------------------------------------------*/

/**
 * Compute Kretschmann scalar:
 *   K = R_{mu,nu,rho,sigma} * R^{mu,nu,rho,sigma}
 *
 * For Schwarzschild: K = 48 * M^2 / r^6
 * Diverges at r=0 (physical singularity).
 */
double compute_kretschmann_scalar(const Tensor4 R_riemann_dn, const Metric *m);

/**
 * Compute Chern-Pontryagin scalar (pseudoscalar):
 *   *R R = (1/2) * epsilon^{mu,nu,rho,sigma} * R_{mu,nu,alpha,beta}
 *          * R_{rho,sigma}^{alpha,beta}
 */
double compute_chern_pontryagin(const Tensor4 R_riemann_dn);

/* ---------------------------------------------------------------------------
 * L4: Bianchi identities
 * -------------------------------------------------------------------------*/

/**
 * Algebraic (first) Bianchi identity:
 *   R^{rho}_{[sigma,mu,nu]} = 0
 *
 * Checks the cyclic sum: R^rho_sigma_mu_nu + R^rho_mu_nu_sigma + R^rho_nu_sigma_mu = 0
 * @return max absolute value of the cyclic sum (should be ≈ 0)
 */
double check_algebraic_bianchi(const Tensor4 R_riemann);

/**
 * Differential (second) Bianchi identity:
 *   ∇_{[lambda} R_{rho,sigma],mu,nu} = 0
 *
 * Simplified: ∇_{rho} G^{rho}_{nu} = 0 (contracted form — energy conservation)
 * @param m        metric
 * @param Gamma    Christoffel symbols
 * @param G        Einstein tensor
 * @param dG       partial derivatives ∂_{rho} G^{mu}_{nu}
 * @return max |∇_{rho} G^{rho}_{nu}| (should be ≈ 0)
 */
double check_contracted_bianchi(const Metric *m, const Tensor3 Gamma,
                                 const Tensor2 G, const Tensor3 dG);

/* ---------------------------------------------------------------------------
 * L4: Einstein tensor trace-reversed form
 * -------------------------------------------------------------------------*/

/**
 * Trace-reversed Einstein tensor (useful in initial value formulation):
 *   \bar{G}_{mu,nu} = G_{mu,nu} - (1/2) * g_{mu,nu} * G^{alpha}_{alpha}
 */
void compute_trace_reversed_einstein(const Metric *m, const Tensor2 G,
                                      Tensor2 G_bar);

/**
 * Compute full curvature structure at a point from metric and its first
 * two derivatives. Bundles all the above computations.
 */
void compute_full_curvature(const Metric *m, const Tensor3 dg,
                             const Tensor4 ddg, Curvature *C);

#endif /* MINI_EINSTEIN_CURVATURE_H */
