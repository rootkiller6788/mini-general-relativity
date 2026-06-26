/**
 * @file tensor.h
 * @brief Tensor algebra for 4-dimensional Lorentzian spacetime.
 *
 * Reference: Wald (1984), Ch.2-3; Carroll (2004), Ch.1-2
 * MIT 8.962 General Relativity
 *
 * All tensors are defined on a 4D pseudo-Riemannian manifold M
 * with signature (-,+,+,+) (East-coast convention).
 *
 * Knowledge map:
 *   L1: Scalar, Vector, Rank-2, Rank-3, Rank-4 tensor typedefs
 *   L3: Tensor contraction, symmetrization, antisymmetrization
 */

#ifndef MINI_EINSTEIN_TENSOR_H
#define MINI_EINSTEIN_TENSOR_H

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * L1: Core tensor type definitions
 * -------------------------------------------------------------------------*/

/** Rank-0 tensor (scalar field) */
typedef double Scalar;

/** Rank-1 tensor (vector / covector), 4 components */
typedef double Vector4[4];

/** Rank-2 tensor, 4x4 components, g_{mu,nu} layout */
typedef double Tensor2[4][4];

/** Rank-3 tensor, 4x4x4 components, Gamma^{lambda}_{mu,nu} layout */
typedef double Tensor3[4][4][4];

/** Rank-4 tensor, 4x4x4x4 components, R^{rho}_{sigma,mu,nu} layout */
typedef double Tensor4[4][4][4][4];

/* ---------------------------------------------------------------------------
 * L1: Minkowski metric (flat spacetime reference)
 * -------------------------------------------------------------------------*/

/** Minkowski metric eta_{mu,nu} = diag(-1, +1, +1, +1) */
static const double MINKOWSKI_METRIC[4][4] = {
    {-1.0,  0.0,  0.0,  0.0},
    { 0.0,  1.0,  0.0,  0.0},
    { 0.0,  0.0,  1.0,  0.0},
    { 0.0,  0.0,  0.0,  1.0}
};

/* ---------------------------------------------------------------------------
 * L1: Kronecker delta and Levi-Civita symbol
 * -------------------------------------------------------------------------*/

/**
 * Kronecker delta delta^{mu}_{nu}.
 * @param mu  row index (0..3)
 * @param nu  column index (0..3)
 * @return 1.0 if mu==nu else 0.0
 */
double kronecker_delta(int mu, int nu);

/**
 * 4D Levi-Civita symbol epsilon_{mu,nu,rho,sigma} (purely covariant).
 * Weight = +1 for even permutations, -1 for odd, 0 otherwise.
 * @return +1.0, -1.0, or 0.0
 */
double levi_civita_4d(int mu, int nu, int rho, int sigma);

/* ---------------------------------------------------------------------------
 * L3: Tensor operations — contraction
 * -------------------------------------------------------------------------*/

/**
 * Contract rank-2 tensor with metric to raise/lower indices:
 *   T^mu_nu = g^{mu,alpha} * T_{alpha,nu}
 * Raises first index of T2.
 * @param g_up    inverse metric (upper indices)
 * @param T_dn    tensor (lower indices)
 * @param result  output: T^{mu}_{nu}
 */
void raise_first_index(const Tensor2 g_up, const Tensor2 T_dn, Tensor2 result);

/**
 * Contract rank-2 tensor with metric to lower first index:
 *   T_{mu,nu} = g_{mu,alpha} * T^{alpha}_{nu}
 * @param g_dn    metric (lower indices)
 * @param T_up    tensor (upper indices)
 * @param result  output: T_{mu,nu}
 */
void lower_first_index(const Tensor2 g_dn, const Tensor2 T_up, Tensor2 result);

/**
 * Lower a vector: v_μ = g_{μ,ν} v^ν
 * @param g_dn    metric (lower indices)
 * @param v_up    contravariant vector v^ν
 * @param v_dn    output covariant vector v_μ
 */
void lower_vector(const Tensor2 g_dn, const Vector4 v_up, Vector4 v_dn);

/**
 * Raise a vector: v^μ = g^{μ,ν} v_ν
 * @param g_up    inverse metric
 * @param v_dn    covariant vector v_ν
 * @param v_up    output contravariant vector v^μ
 */
void raise_vector(const Tensor2 g_up, const Vector4 v_dn, Vector4 v_up);

/**
 * Full double-trace of two rank-2 tensors: A^{mu,nu} * B_{mu,nu}
 * @param A  first tensor (upper indices)
 * @param B  second tensor (lower indices)
 * @return scalar trace sum_{mu,nu} A^{mu,nu} * B_{mu,nu}
 */
double double_contraction(const Tensor2 A, const Tensor2 B);

/**
 * Contract rank-3 with rank-1: C^lambda = Gamma^{lambda}_{mu,nu} * U^{mu} * U^{nu}
 * Appears in geodesic equation: d^2 x^lambda/dtau^2 + Gamma^lambda_{mu,nu} U^mu U^nu = 0
 * @param Gamma  Christoffel symbols (upper first index)
 * @param U      four-velocity
 * @param result output C[0..3]
 */
void contract_gamma_uu(const Tensor3 Gamma, const Vector4 U, Vector4 result);

/**
 * Riemann tensor contraction to Ricci: R_{mu,nu} = R^{rho}_{mu,rho,nu}
 * @param R_riemann  Riemann tensor R^rho_sigma_mu_nu
 * @param R_ricci    output Ricci tensor R_{mu,nu}
 */
void riemann_to_ricci(const Tensor4 R_riemann, Tensor2 R_ricci);

/**
 * Ricci scalar: R = g^{mu,nu} * R_{mu,nu}
 * @param g_up    inverse metric
 * @param R_ricci Ricci tensor (lower indices)
 * @return Ricci scalar R
 */
double ricci_scalar_contract(const Tensor2 g_up, const Tensor2 R_ricci);

/* ---------------------------------------------------------------------------
 * L3: Tensor operations — symmetrization
 * -------------------------------------------------------------------------*/

/**
 * Symmetrize rank-2 tensor in place: T_{(mu,nu)} = (T_{mu,nu}+T_{nu,mu})/2
 * @param T  tensor to symmetrize (modified in place)
 */
void symmetrize_tensor2(Tensor2 T);

/**
 * Antisymmetrize rank-2 tensor in place: T_{[mu,nu]} = (T_{mu,nu}-T_{nu,mu})/2
 * @param T  tensor to antisymmetrize (modified in place)
 */
void antisymmetrize_tensor2(Tensor2 T);

/**
 * Compute trace of rank-2 tensor: T = T^{mu}_{mu}
 * @param T  tensor (mixed indices, i.e. one up one down)
 * @return trace
 */
double trace_tensor2(const Tensor2 T);

/* ---------------------------------------------------------------------------
 * L3: Tensor operations — algebra utilities
 * -------------------------------------------------------------------------*/

/**
 * Copy tensor T2: dst = src
 */
void copy_tensor2(const Tensor2 src, Tensor2 dst);

/**
 * Scale tensor T2: T *= alpha
 */
void scale_tensor2(Tensor2 T, double alpha);

/**
 * Add two tensors: C = A + B
 */
void add_tensor2(const Tensor2 A, const Tensor2 B, Tensor2 C);

/**
 * Subtract two tensors: C = A - B
 */
void subtract_tensor2(const Tensor2 A, const Tensor2 B, Tensor2 C);

/**
 * Set tensor2 to identity (Kronecker delta): T[mu][nu] = delta_{mu,nu}
 */
void set_identity_tensor2(Tensor2 T);

/**
 * Zero a rank-2 tensor
 */
void zero_tensor2(Tensor2 T);

/**
 * Zero a rank-3 tensor
 */
void zero_tensor3(Tensor3 T);

/**
 * Zero a rank-4 tensor
 */
void zero_tensor4(Tensor4 T);

/**
 * Frobenius norm of rank-2 tensor: sqrt(sum T_{mu,nu}^2)
 */
double norm_tensor2(const Tensor2 T);

/**
 * Check if tensor2 is symmetric to within tolerance
 * @return 1 if symmetric, 0 otherwise
 */
int is_symmetric_tensor2(const Tensor2 T, double tol);

/**
 * Matrix multiplication of two 4x4 tensors: C^{mu}_{nu} = A^{mu}_{alpha} * B^{alpha}_{nu}
 * Both inputs treated as mixed-index tensors.
 */
void matmul_tensor2(const Tensor2 A, const Tensor2 B, Tensor2 C);

/**
 * Outer product of two vectors: T^{mu,nu} = U^{mu} * V^{nu}
 */
void outer_product_vector4(const Vector4 U, const Vector4 V, Tensor2 T);

#endif /* MINI_EINSTEIN_TENSOR_H */
