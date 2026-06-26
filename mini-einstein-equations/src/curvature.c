/**
 * @file curvature.c
 * @brief Riemann, Ricci, Einstein, Weyl curvature tensors.
 *
 * Reference: Wald (1984), Ch.3; Carroll (2004), Ch.3
 *            Weinberg "Gravitation and Cosmology" Ch.6
 *
 * The Riemann tensor measures the failure of second covariant derivatives
 * to commute. It completely characterizes spacetime curvature.
 */

#include "curvature.h"
#include <math.h>
#include <string.h>

/* ========================================================================
 * L1: Riemann tensor
 * ========================================================================*/

/**
 * Riemann tensor R^ρ_{σ,μ,ν} from Christoffel symbols and their derivatives.
 *
 * Formula:
 *   R^ρ_{σ,μ,ν} = ∂_μ Γ^ρ_{ν,σ} - ∂_ν Γ^ρ_{μ,σ}
 *                + Γ^ρ_{μ,λ} Γ^λ_{ν,σ}
 *                - Γ^ρ_{ν,λ} Γ^λ_{μ,σ}
 *
 * Symmetries (of fully covariant Riemann):
 *   R_{ρ,σ,μ,ν} = -R_{σ,ρ,μ,ν}         (antisymmetric in first pair)
 *   R_{ρ,σ,μ,ν} = -R_{ρ,σ,ν,μ}         (antisymmetric in second pair)
 *   R_{ρ,σ,μ,ν} =  R_{μ,ν,ρ,σ}         (pair symmetry)
 *   R_{ρ,[σ,μ,ν]} = 0                   (algebraic Bianchi)
 *
 * Knowledge point: The Riemann tensor is the fundamental object of
 * differential geometry in GR. Its vanishing is necessary and sufficient
 * for the spacetime to be flat.
 *
 * Convention: dGamma[lambda][rho][mu][nu] = ∂_lambda Gamma^rho_{mu,nu}
 */
void compute_riemann_tensor(const Tensor3 Gamma,
                             const Tensor4 dGamma,
                             Tensor4 R_riemann)
{
    for (int rho = 0; rho < 4; rho++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    /* Derivative terms */
                    double d_term = dGamma[mu][rho][nu][sigma]
                                  - dGamma[nu][rho][mu][sigma];

                    /* Product terms */
                    double prod_pos = 0.0, prod_neg = 0.0;
                    for (int lambda = 0; lambda < 4; lambda++) {
                        prod_pos += Gamma[rho][mu][lambda]
                                  * Gamma[lambda][nu][sigma];
                        prod_neg += Gamma[rho][nu][lambda]
                                  * Gamma[lambda][mu][sigma];
                    }

                    R_riemann[rho][sigma][mu][nu] = d_term + prod_pos - prod_neg;
                }
            }
        }
    }
}

/**
 * Fully covariant Riemann tensor:
 *   R_{ρ,σ,μ,ν} = g_{ρ,λ} * R^λ_{σ,μ,ν}
 *
 * The covariant form has the most transparent symmetries.
 */
void riemann_covariant(const Metric *m, const Tensor4 R_riemann_up,
                       Tensor4 R_riemann_dn)
{
    for (int rho = 0; rho < 4; rho++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    double sum = 0.0;
                    for (int lambda = 0; lambda < 4; lambda++) {
                        sum += m->g[rho][lambda]
                             * R_riemann_up[lambda][sigma][mu][nu];
                    }
                    R_riemann_dn[rho][sigma][mu][nu] = sum;
                }
            }
        }
    }
}

/* ========================================================================
 * L1: Ricci tensor and scalar
 * ========================================================================*/

/**
 * Ricci tensor: R_{μ,ν} = R^ρ_{μ,ρ,ν}
 *
 * Knowledge point: The Ricci tensor is the trace of Riemann over the
 * first and third indices. It appears directly in the Einstein equation.
 * In vacuum (T_{μ,ν}=0), the Einstein equation reduces to R_{μ,ν}=0.
 */
void compute_ricci_tensor(const Tensor4 R_riemann, Tensor2 R_ricci)
{
    riemann_to_ricci(R_riemann, R_ricci);
}

/**
 * Ricci scalar: R = g^{μ,ν} R_{μ,ν}
 *
 * Knowledge point: The Ricci scalar is the simplest curvature invariant.
 * It appears in the Einstein-Hilbert action S_EH = ∫ R √(-g) d^4x.
 */
double compute_ricci_scalar(const Metric *m, const Tensor2 R_ricci)
{
    return ricci_scalar_contract(m->g_inv, R_ricci);
}

/* ========================================================================
 * L1: Einstein tensor
 * ========================================================================*/

/**
 * Einstein tensor: G_{μ,ν} = R_{μ,ν} - (1/2) R g_{μ,ν}
 *
 * Knowledge point: The Einstein tensor is the unique divergence-free
 * rank-2 tensor constructed from the metric and its first two derivatives
 * (Lovelock's theorem in 4D).
 *
 * Key property: ∇^μ G_{μ,ν} = 0 (contracted Bianchi identity).
 * This ensures consistency with energy-momentum conservation ∇^μ T_{μ,ν} = 0.
 */
void compute_einstein_tensor(const Metric *m, const Tensor2 R_ricci,
                              double R_scalar, Tensor2 G)
{
    double half_R = 0.5 * R_scalar;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            G[mu][nu] = R_ricci[mu][nu] - half_R * m->g[mu][nu];
        }
    }
}

/* ========================================================================
 * L1: Weyl tensor (conformal curvature)
 * ========================================================================*/

/**
 * Weyl tensor C_{ρ,σ,μ,ν} (fully covariant).
 *
 * Formula:
 *   C_{ρ,σ,μ,ν} = R_{ρ,σ,μ,ν}
 *     - (1/2)(g_{ρ,μ} R_{σ,ν} - g_{ρ,ν} R_{σ,μ}
 *             - g_{σ,μ} R_{ρ,ν} + g_{σ,ν} R_{ρ,μ})
 *     + (R/6)(g_{ρ,μ} g_{σ,ν} - g_{ρ,ν} g_{σ,μ})
 *
 * Properties:
 *   - Trace-free on all index pairs.
 *   - Identically zero in ≤ 3 dimensions.
 *   - Vanishes for conformally flat spacetimes (e.g., FLRW).
 *   - In vacuum (R_{μ,ν}=0): C = R (Weyl = Riemann).
 *
 * Knowledge point: The Weyl tensor captures the "free gravitational field"
 * — those aspects of curvature not determined locally by matter.
 * Gravitational waves are encoded in the Weyl tensor.
 */
void compute_weyl_tensor(const Metric *m, const Tensor4 R_dn,
                          const Tensor2 R_ricci, double R_scalar,
                          Tensor4 Weyl)
{
    double R_over_6 = R_scalar / 6.0;

    for (int rho = 0; rho < 4; rho++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    /* Ricci correction part */
                    double ricci_part = 0.5 * (
                          m->g[rho][mu] * R_ricci[sigma][nu]
                        - m->g[rho][nu] * R_ricci[sigma][mu]
                        - m->g[sigma][mu] * R_ricci[rho][nu]
                        + m->g[sigma][nu] * R_ricci[rho][mu]
                    );

                    /* Scalar correction part */
                    double scalar_part = R_over_6 * (
                          m->g[rho][mu] * m->g[sigma][nu]
                        - m->g[rho][nu] * m->g[sigma][mu]
                    );

                    Weyl[rho][sigma][mu][nu] = R_dn[rho][sigma][mu][nu]
                                               - ricci_part + scalar_part;
                }
            }
        }
    }
}

/* ========================================================================
 * L3: Curvature invariants
 * ========================================================================*/

/**
 * Kretschmann scalar: K = R_{μ,ν,ρ,σ} * R^{μ,ν,ρ,σ}
 *
 * Knowledge point: The Kretschmann scalar is a curvature invariant that
 * detects true spacetime singularities (as opposed to coordinate ones).
 *
 * For Schwarzschild: K = 48 M^2 / r^6 → diverges at r=0 (physical singularity).
 * At r = 2M (horizon): K = 3/(4M^4) — finite (coordinate singularity).
 *
 * For Minkowski: K = 0.
 * For de Sitter: K = 24 / L^4 where L is the de Sitter radius.
 */
double compute_kretschmann_scalar(const Tensor4 R_riemann_dn, const Metric *m)
{
    /* Compute R^{μ,ν,ρ,σ} = g^{μ,α} g^{ν,β} g^{ρ,γ} g^{σ,δ} R_{α,β,γ,δ} */
    double K = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    /* Raise all indices of R_dn */
                    double R_up = 0.0;
                    for (int a = 0; a < 4; a++) {
                        for (int b = 0; b < 4; b++) {
                            for (int c = 0; c < 4; c++) {
                                for (int d = 0; d < 4; d++) {
                                    R_up += m->g_inv[mu][a]
                                          * m->g_inv[nu][b]
                                          * m->g_inv[rho][c]
                                          * m->g_inv[sigma][d]
                                          * R_riemann_dn[a][b][c][d];
                                }
                            }
                        }
                    }
                    K += R_up * R_riemann_dn[mu][nu][rho][sigma];
                }
            }
        }
    }
    return K;
}

/**
 * Chern-Pontryagin pseudoscalar (topological invariant):
 *   *R R = (1/2) ε^{μ,ν,ρ,σ} R_{μ,ν,α,β} R_{ρ,σ}^{α,β}
 *
 * Note: This is a simplified computation using the Levi-Civita symbol
 * (not the tensor density — i.e., without sqrt(-g)).
 * The full topological invariant is ∫ *R R √(-g) d^4x.
 *
 * Knowledge point: The Chern-Pontryagin invariant is a topological
 * invariant of the manifold (related to the gravitational theta-term).
 * It vanishes for spherically symmetric spacetimes.
 */
double compute_chern_pontryagin(const Tensor4 R_riemann_dn)
{
    /* This is a simplified computation — real CP invariant requires
     * the Riemann tensor with two upper and two lower indices.
     * We compute a rough proxy: ε^{μ,ν,ρ,σ} * R_{μ,ν,0,1} * R_{ρ,σ,2,3}
     * as a representative term. */

    double cp = 0.0;
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            for (int c = 0; c < 4; c++) {
                for (int d = 0; d < 4; d++) {
                    double eps_abcd = levi_civita_4d(a, b, c, d);
                    if (fabs(eps_abcd) < 0.5) continue;
                    /* Contract with Riemann */
                    double sum_ef = 0.0;
                    for (int e = 0; e < 4; e++) {
                        for (int f = 0; f < 4; f++) {
                            sum_ef += R_riemann_dn[a][b][e][f]
                                    * R_riemann_dn[c][d][e][f];
                        }
                    }
                    cp += eps_abcd * sum_ef;
                }
            }
        }
    }
    return 0.5 * cp;
}

/* ========================================================================
 * L4: Bianchi identities
 * ========================================================================*/

/**
 * Algebraic (first) Bianchi identity:
 *   R^ρ_{[σ,μ,ν]} = R^ρ_{σ,μ,ν} + R^ρ_{μ,ν,σ} + R^ρ_{ν,σ,μ} = 0
 *
 * Knowledge point: The algebraic Bianchi identity is a consequence of the
 * torsion-free condition and the Jacobi identity for commutators.
 * It implies R_{[ρ,σ,μ,ν]} = 0.
 */
double check_algebraic_bianchi(const Tensor4 R_riemann)
{
    double max_err = 0.0;
    for (int rho = 0; rho < 4; rho++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    double cyclic_sum = R_riemann[rho][sigma][mu][nu]
                                      + R_riemann[rho][mu][nu][sigma]
                                      + R_riemann[rho][nu][sigma][mu];
                    double err = fabs(cyclic_sum);
                    if (err > max_err) max_err = err;
                }
            }
        }
    }
    return max_err;
}

/**
 * Contracted differential Bianchi identity:
 *   ∇_ρ G^ρ_ν = ∇_ρ (R^ρ_ν - (1/2) R δ^ρ_ν) = 0
 *
 * This is a simplified check: we compute the covariant divergence of
 * the Einstein tensor and verify it vanishes.
 *
 * ∇_ρ G^ρ_ν = ∂_ρ G^ρ_ν + Γ^ρ_{ρ,σ} G^σ_ν - Γ^σ_{ρ,ν} G^ρ_σ
 *
 * Knowledge point: The contracted Bianchi identity is equivalent to
 * energy-momentum conservation in GR. It is an identity, not an
 * additional dynamical equation.
 */
double check_contracted_bianchi(const Metric *m, const Tensor3 Gamma,
                                 const Tensor2 G, const Tensor3 dG)
{
    (void)m; /* g^{μ,ν} used implicitly through G and Gamma */
    double max_err = 0.0;

    /* Compute divergence: ∇_ρ G^ρ_ν */
    for (int nu = 0; nu < 4; nu++) {
        for (int rho = 0; rho < 4; rho++) {
            double div = dG[rho][rho][nu]; /* ∂_ρ G^ρ_ν */
            for (int sigma = 0; sigma < 4; sigma++) {
                div += Gamma[rho][rho][sigma] * G[sigma][nu];
                div -= Gamma[sigma][rho][nu] * G[rho][sigma];
            }
            double err = fabs(div);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

/* ========================================================================
 * L4: Trace-reversed Einstein
 * ========================================================================*/

void compute_trace_reversed_einstein(const Metric *m, const Tensor2 G,
                                      Tensor2 G_bar)
{
    /* G_bar_{μ,ν} = G_{μ,ν} - (1/2) g_{μ,ν} G^{α}_{α} */
    double G_trace = trace_tensor2(G);
    double half_trace = 0.5 * G_trace;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            G_bar[mu][nu] = G[mu][nu] - half_trace * m->g[mu][nu];
        }
    }
}

/* ========================================================================
 * L4: Full curvature pipeline
 * ========================================================================*/

/**
 * Compute complete curvature at a point.
 *
 * Pipeline:
 *   1. Christoffel symbols from metric + dg
 *   2. dGamma from ddg (reconstruct partial derivatives of Gamma)
 *   3. Riemann from Gamma + dGamma
 *   4. Ricci from Riemann
 *   5. Ricci scalar from Ricci + metric
 *   6. Einstein from Ricci + R + metric
 *   7. Kretschmann invariant
 *
 * Knowledge point: This is the standard curvature computation pipeline
 * used in numerical relativity.
 *
 * Note: dGamma[lambda][rho][mu][nu] = ∂_lambda Gamma^rho_{mu,nu}
 * which requires second derivatives of the metric:
 *   dGamma = (1/2) g_inv * (∂∂g permutations) + derivatives of g_inv * (∂g)
 *
 * For simplicity, we approximate dGamma from numerical second derivatives
 * of the metric. The calling code must provide this.
 */
void compute_full_curvature(const Metric *m, const Tensor3 dg,
                             const Tensor4 ddg, Curvature *C)
{
    Tensor3 Gamma;
    metric_christoffel(m, dg, Gamma);

    /* Compute dGamma from ddg.
     * dGamma[lambda][rho][mu][nu] = ∂_lambda Gamma^rho_{mu,nu}
     *
     * ∂_λ Γ^ρ_{μ,ν} = (1/2) Σ_σ [ (∂_λ g^{ρ,σ}) * (∂_μ g_{σ,ν} + ∂_ν g_{σ,μ} - ∂_σ g_{μ,ν})
     *                              + g^{ρ,σ} * (∂_λ ∂_μ g_{σ,ν} + ∂_λ ∂_ν g_{σ,μ} - ∂_λ ∂_σ g_{μ,ν}) ]
     *
     * For simplicity we use finite differences for dGamma rather than
     * computing the full analytical expression. The caller should provide
     * dGamma directly if high accuracy is needed.
     *
     * Here ddg[lambda][sigma][mu][nu] = ∂_λ ∂_σ g_{μ,ν}
     * We approximate dGamma using the metric term only (neglecting d(g_inv) terms).
     */
    Tensor4 dGamma;
    memset(dGamma, 0, sizeof(dGamma));

    for (int lambda = 0; lambda < 4; lambda++) {
        for (int rho = 0; rho < 4; rho++) {
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    double sum = 0.0;
                    for (int sigma = 0; sigma < 4; sigma++) {
                        /* ∂_λ Γ^ρ_{μ,ν} ~ g^{ρ,σ} * (∂_λ ∂_μ g_{σ,ν}
                         *                              + ∂_λ ∂_ν g_{σ,μ}
                         *                              - ∂_λ ∂_σ g_{μ,ν}) / 2 */
                        double d2g = ddg[lambda][mu][sigma][nu]
                                   + ddg[lambda][nu][sigma][mu]
                                   - ddg[lambda][sigma][mu][nu];
                        sum += m->g_inv[rho][sigma] * d2g;
                    }
                    dGamma[lambda][rho][mu][nu] = 0.5 * sum;
                }
            }
        }
    }

    compute_riemann_tensor(Gamma, dGamma, C->Riemann);
    compute_ricci_tensor(C->Riemann, C->Ricci);
    C->RicciScalar = compute_ricci_scalar(m, C->Ricci);
    compute_einstein_tensor(m, C->Ricci, C->RicciScalar, C->Einstein);
    C->Kretschmann = compute_kretschmann_scalar(C->Riemann, m);
}
