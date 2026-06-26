/**
 * @file tensor.c
 * @brief Tensor algebra implementation for 4D Lorentzian spacetime.
 *
 * Reference: Wald (1984), Ch.2-3; Carroll (2004), Ch.1-2
 *
 * All operations work on double-precision arrays representing
 * rank-0 through rank-4 tensors in 4-dimensional spacetime.
 *
 * Complexity note: Most operations are O(1) since dimension is fixed at 4.
 * Contractions are O(D^2) to O(D^4) but D=4 is constant.
 */

#include "tensor.h"
#include <math.h>
#include <string.h>

/* ========================================================================
 * L1: Kronecker delta and Levi-Civita symbol
 * ========================================================================*/

double kronecker_delta(int mu, int nu)
{
    return (mu == nu) ? 1.0 : 0.0;
}

/**
 * 4D Levi-Civita symbol ε_{μνρσ}.
 *
 * Definition: ε_{0123} = +1, totally antisymmetric.
 * Even permutation → +1, odd → -1, repeated index → 0.
 *
 * This is the pure number (not the tensor density), i.e.,
 * it does NOT include sqrt(-g).
 *
 * Ref: Carroll (2004) Appendix B; Wald (1984) Appendix B.
 */
double levi_civita_4d(int mu, int nu, int rho, int sigma)
{
    /* All indices must be distinct */
    if (mu == nu || mu == rho || mu == sigma ||
        nu == rho || nu == sigma || rho == sigma)
        return 0.0;

    /* Permutation parity via inversion count */
    int idx[4] = {mu, nu, rho, sigma};
    int inversions = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            if (idx[i] > idx[j]) inversions++;
        }
    }
    return (inversions % 2 == 0) ? 1.0 : -1.0;
}

/* ========================================================================
 * L3: Tensor contraction operations
 * ========================================================================*/

/**
 * Raise first index: T^μ_ν = g^{μ,α} * T_{α,ν}
 *
 * Knowledge point: Index raising/lowering via metric.
 * This is the fundamental operation that maps between tangent
 * and cotangent spaces in a pseudo-Riemannian manifold.
 */
void raise_first_index(const Tensor2 g_up, const Tensor2 T_dn, Tensor2 result)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                sum += g_up[mu][alpha] * T_dn[alpha][nu];
            }
            result[mu][nu] = sum;
        }
    }
}

/**
 * Lower first index: T_{μ,ν} = g_{μ,α} * T^{α}_{ν}
 */
void lower_first_index(const Tensor2 g_dn, const Tensor2 T_up, Tensor2 result)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                sum += g_dn[mu][alpha] * T_up[alpha][nu];
            }
            result[mu][nu] = sum;
        }
    }
}

/**
 * Lower a vector: v_μ = g_{μ,ν} v^ν
 *
 * Knowledge point: Index lowering maps vectors between tangent and
 * cotangent spaces. This is the most basic metric operation.
 */
void lower_vector(const Tensor2 g_dn, const Vector4 v_up, Vector4 v_dn)
{
    for (int mu = 0; mu < 4; mu++) {
        double sum = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            sum += g_dn[mu][nu] * v_up[nu];
        }
        v_dn[mu] = sum;
    }
}

/**
 * Raise a vector: v^μ = g^{μ,ν} v_ν
 */
void raise_vector(const Tensor2 g_up, const Vector4 v_dn, Vector4 v_up)
{
    for (int mu = 0; mu < 4; mu++) {
        double sum = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            sum += g_up[mu][nu] * v_dn[nu];
        }
        v_up[mu] = sum;
    }
}

/**
 * Double contraction A^{μ,ν} * B_{μ,ν}.
 *
 * Knowledge point: Full trace — the canonical pairing between
 * a tensor and its dual. Used in action principles.
 */
double double_contraction(const Tensor2 A, const Tensor2 B)
{
    double sum = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            sum += A[mu][nu] * B[mu][nu];
        }
    }
    return sum;
}

/**
 * Contract Christoffel symbols with four-velocity:
 *   C^λ = Γ^λ_{μ,ν} * U^μ * U^ν
 *
 * Knowledge point: Gamma-U-U contraction — the nonlinear term
 * in the geodesic equation.
 */
void contract_gamma_uu(const Tensor3 Gamma, const Vector4 U, Vector4 result)
{
    for (int lambda = 0; lambda < 4; lambda++) {
        double sum = 0.0;
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                sum += Gamma[lambda][mu][nu] * U[mu] * U[nu];
            }
        }
        result[lambda] = sum;
    }
}

/**
 * Riemann to Ricci contraction: R_{μ,ν} = R^ρ_{μ,ρ,ν}
 *
 * Knowledge point: The Ricci tensor as a trace of the Riemann tensor.
 * This is the only non-vanishing contraction of Riemann (due to its
 * algebraic symmetries — all others vanish or are proportional).
 */
void riemann_to_ricci(const Tensor4 R_riemann, Tensor2 R_ricci)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int rho = 0; rho < 4; rho++) {
                sum += R_riemann[rho][mu][rho][nu];
            }
            R_ricci[mu][nu] = sum;
        }
    }
}

/**
 * Ricci scalar: R = g^{μ,ν} * R_{μ,ν}
 *
 * Knowledge point: The curvature scalar — the simplest
 * diffeomorphism-invariant curvature quantity.
 */
double ricci_scalar_contract(const Tensor2 g_up, const Tensor2 R_ricci)
{
    double R = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            R += g_up[mu][nu] * R_ricci[mu][nu];
        }
    }
    return R;
}

/* ========================================================================
 * L3: Symmetrization / antisymmetrization
 * ========================================================================*/

/**
 * Symmetrize rank-2 tensor: T_{(μ,ν)} = (T_{μ,ν} + T_{ν,μ}) / 2
 *
 * Knowledge point: Symmetric part of a tensor.
 * The metric g_{μ,ν}, Ricci R_{μ,ν}, Einstein G_{μ,ν}, and
 * stress-energy T_{μ,ν} are all symmetric rank-2 tensors.
 */
void symmetrize_tensor2(Tensor2 T)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            double sym = 0.5 * (T[mu][nu] + T[nu][mu]);
            T[mu][nu] = sym;
            T[nu][mu] = sym;
        }
    }
}

/**
 * Antisymmetrize rank-2 tensor: T_{[μ,ν]} = (T_{μ,ν} - T_{ν,μ}) / 2
 *
 * Knowledge point: Antisymmetric part.
 * The Faraday tensor F_{μ,ν} is antisymmetric.
 */
void antisymmetrize_tensor2(Tensor2 T)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            double anti = 0.5 * (T[mu][nu] - T[nu][mu]);
            T[mu][nu] =  anti;
            T[nu][mu] = -anti;
        }
    }
    /* Diagonal must be zero for antisymmetric tensor */
    for (int mu = 0; mu < 4; mu++) {
        T[mu][mu] = 0.0;
    }
}

/**
 * Trace of rank-2 tensor: T^μ_μ
 *
 * Knowledge point: The trace is invariant under change of basis.
 * Einstein equation trace: -R + 4Λ = κ T^μ_μ
 */
double trace_tensor2(const Tensor2 T)
{
    double trace = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        trace += T[mu][mu];
    }
    return trace;
}

/* ========================================================================
 * L3: Tensor utility functions
 * ========================================================================*/

void copy_tensor2(const Tensor2 src, Tensor2 dst)
{
    memcpy(dst, src, sizeof(Tensor2));
}

void scale_tensor2(Tensor2 T, double alpha)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] *= alpha;
        }
    }
}

void add_tensor2(const Tensor2 A, const Tensor2 B, Tensor2 C)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            C[mu][nu] = A[mu][nu] + B[mu][nu];
        }
    }
}

void subtract_tensor2(const Tensor2 A, const Tensor2 B, Tensor2 C)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            C[mu][nu] = A[mu][nu] - B[mu][nu];
        }
    }
}

void set_identity_tensor2(Tensor2 T)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] = (mu == nu) ? 1.0 : 0.0;
        }
    }
}

void zero_tensor2(Tensor2 T)
{
    memset(T, 0, sizeof(Tensor2));
}

void zero_tensor3(Tensor3 T)
{
    memset(T, 0, sizeof(Tensor3));
}

void zero_tensor4(Tensor4 T)
{
    memset(T, 0, sizeof(Tensor4));
}

/**
 * Frobenius norm of rank-2 tensor.
 *
 * Knowledge point: The Frobenius norm provides a natural scalar
 * measure of tensor magnitude, used in convergence checks.
 */
double norm_tensor2(const Tensor2 T)
{
    double sum_sq = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double v = T[mu][nu];
            sum_sq += v * v;
        }
    }
    return sqrt(sum_sq);
}

/**
 * Symmetry check for rank-2 tensor.
 *
 * Knowledge point: Physical rank-2 tensors (g, R, G, T) must be symmetric.
 * This is a consistency check for numerical computations.
 */
int is_symmetric_tensor2(const Tensor2 T, double tol)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            if (fabs(T[mu][nu] - T[nu][mu]) > tol) return 0;
        }
    }
    return 1;
}

/**
 * Matrix multiplication: C = A * B (both treated as 4×4 matrices)
 *
 * Knowledge point: The composition of two linear maps.
 * Used to verify: g^{μ,α} * g_{α,ν} = δ^μ_ν.
 */
void matmul_tensor2(const Tensor2 A, const Tensor2 B, Tensor2 C)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                sum += A[mu][alpha] * B[alpha][nu];
            }
            C[mu][nu] = sum;
        }
    }
}

/**
 * Outer product: T^{μ,ν} = U^μ * V^ν
 *
 * Knowledge point: Dyadic product — a rank-2 tensor built from
 * two vectors. Used in constructing stress-energy of dust.
 */
void outer_product_vector4(const Vector4 U, const Vector4 V, Tensor2 T)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] = U[mu] * V[nu];
        }
    }
}
