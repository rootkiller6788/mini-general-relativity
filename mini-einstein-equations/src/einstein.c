/**
 * @file einstein.c
 * @brief Einstein field equations: evaluation, trace forms, action.
 *
 * Reference: Wald (1984), Ch.4, Ch.E (Appendix E)
 *            Carroll (2004), Ch.4
 *            Weinberg "Gravitation and Cosmology" Ch.7
 *
 * The Einstein field equations are the central equations of general
 * relativity, relating spacetime curvature to matter content.
 */

#include "einstein.h"
#include <math.h>
#include <string.h>

/* ========================================================================
 * L4: Einstein equation residuals and components
 * ========================================================================*/

/**
 * Einstein field equations residual:
 *   residual_{μ,ν} = G_{μ,ν} + Λ g_{μ,ν} - κ T_{μ,ν}
 *
 * A valid solution has residual = 0 (to within numerical tolerance).
 *
 * Knowledge point: The Einstein equation is a system of 10 coupled,
 * nonlinear, second-order PDEs for the metric components.
 * Only 10 of the 16 components are independent due to symmetry.
 * The Bianchi identities provide 4 constraints, leaving 6 dynamical degrees
 * of freedom (2 tensor polarizations of gravitational waves).
 */
void einstein_equation_residual(const Metric *m, const Tensor2 G,
                                 double Lambda, double kappa,
                                 const Tensor2 T, Tensor2 residual)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            residual[mu][nu] = G[mu][nu]
                             + Lambda * m->g[mu][nu]
                             - kappa * T[mu][nu];
        }
    }
}

void einstein_lhs(const Metric *m, const Tensor2 G, double Lambda,
                  Tensor2 E)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            E[mu][nu] = G[mu][nu] + Lambda * m->g[mu][nu];
        }
    }
}

void einstein_rhs(double kappa, const Tensor2 T, Tensor2 S)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            S[mu][nu] = kappa * T[mu][nu];
        }
    }
}

/**
 * Trace of Einstein equation:
 *   -R + 4Λ = κ T
 *
 * Derivation:
 *   G_{μ,ν} + Λ g_{μ,ν} = κ T_{μ,ν}
 *   Trace: G^μ_μ + 4Λ = κ T^μ_μ
 *   But G^μ_μ = R^μ_μ - (1/2)R δ^μ_μ = R - 2R = -R
 *   So: -R + 4Λ = κ T
 *
 * Knowledge point: The trace equation is algebraically simpler and
 * useful for substituting R in terms of T in the full equations.
 * For vacuum (T=0, Λ=0): R=0 (but R_{μ,ν} may still be non-zero
 * via the Weyl tensor — gravitational waves carry Ricci-flat curvature).
 */
double einstein_trace_equation(double R_scalar, double Lambda,
                                double T_trace, double kappa)
{
    return -R_scalar + 4.0 * Lambda - kappa * T_trace;
}

/**
 * Trace-reversed Einstein equation:
 *   R_{μ,ν} = κ (T_{μ,ν} - (1/2) g_{μ,ν} T)
 *
 * Derivation:
 *   Start with G_{μ,ν} = κ T_{μ,ν}
 *   G_{μ,ν} = R_{μ,ν} - (1/2) R g_{μ,ν}
 *   Substitute R = -κ T from trace equation:
 *   R_{μ,ν} - (1/2)(-κ T) g_{μ,ν} = κ T_{μ,ν}
 *   R_{μ,ν} = κ T_{μ,ν} - (κ T/2) g_{μ,ν}
 *
 * This form is often more convenient for solving the equations because
 * the Einstein tensor's nonlinearity is absorbed into the simpler Ricci
 * tensor, once the trace is known.
 *
 * Knowledge point: The trace-reversed form is the starting point for
 * the weak-field approximation, post-Newtonian expansion, and the
 * initial value formulation (ADM decomposition).
 */
void einstein_trace_reversed(const Metric *m, const Tensor2 T,
                              double T_trace, double kappa, Tensor2 R_ricci)
{
    double kappa_T_over_2 = 0.5 * kappa * T_trace;

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            R_ricci[mu][nu] = kappa * T[mu][nu]
                            - kappa_T_over_2 * m->g[mu][nu];
        }
    }
}

void einstein_from_matter(const Metric *m, const Tensor2 T,
                           double kappa, double Lambda, Tensor2 G)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            G[mu][nu] = kappa * T[mu][nu] - Lambda * m->g[mu][nu];
        }
    }
}

/* ========================================================================
 * L4: Einstein-Hilbert action
 * ========================================================================*/

/**
 * Gravitational Lagrangian density: L_G = R √(-g)
 *
 * The Einstein-Hilbert action:
 *   S_EH = (c^4 / (16πG)) ∫ R √(-g) d^4x
 *
 * Variation δS_EH = 0 with respect to g^{μ,ν} yields the vacuum
 * Einstein equations G_{μ,ν} = 0.
 *
 * Knowledge point: The Einstein-Hilbert action is the simplest
 * diffeomorphism-invariant action that yields second-order field equations.
 * Lovelock's theorem (1971) proves it's unique in 4D.
 *
 * The √(-g) factor ensures the action is a scalar (d^4x alone is not
 * invariant — it transforms with the Jacobian of the coordinate change,
 * while √(-g) transforms with the inverse Jacobian, canceling it).
 */
double einstein_hilbert_lagrangian_density(double R_scalar, double det_g)
{
    return R_scalar * sqrt(-det_g);
}

double matter_lagrangian_density_perfect_fluid(double rho, double det_g)
{
    /* For pressureless dust: L_M = -ρ √(-g) */
    return -rho * sqrt(-det_g);
}

double total_lagrangian_density(double R_scalar, double rho, double det_g)
{
    return einstein_hilbert_lagrangian_density(R_scalar, det_g)
         + matter_lagrangian_density_perfect_fluid(rho, det_g);
}

/* ========================================================================
 * L7: Vacuum equation check
 * ========================================================================*/

/**
 * Check vacuum Einstein equations: G_{μ,ν} + Λ g_{μ,ν} = 0
 *
 * For Λ = 0: G_{μ,ν} = 0 ⇔ R_{μ,ν} = 0
 * For Λ ≠ 0: G_{μ,ν} = -Λ g_{μ,ν} ⇔ R_{μ,ν} = Λ g_{μ,ν} (de Sitter-like)
 *
 * Returns the max absolute component of the residual.
 */
double einstein_vacuum_check(const Metric *m, const Tensor2 G, double Lambda)
{
    double max_err = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double residual = G[mu][nu] + Lambda * m->g[mu][nu];
            double err = fabs(residual);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

/**
 * Komar mass (approximate for Schwarzschild-like spacetimes).
 *
 * The full Komar integral:
 *   M = (1/8π) ∫_{S_∞} ∇^μ ξ^ν dS_{μ,ν}
 *
 * For Schwarzschild in Schwarzschild coordinates, with the timelike
 * Killing vector ξ = ∂_t, the only non-vanishing component of ∇^μ ξ^ν
 * at large r gives M.
 *
 * This simplified version computes the Komar mass from the derivative
 * of g_{00} at a given radius:
 *   M ≈ (r^2 / 2G) * ∂_r g_{00}  (in geometric units)
 *
 * Knowledge point: The Komar mass is a conserved quantity associated
 * with the timelike Killing vector in stationary spacetimes.
 * It equals the ADM mass for asymptotically flat spacetimes.
 */
double komar_mass(const Metric *m, const Tensor3 Gamma,
                  const Vector4 xi, double radius)
{
    /* Approximate: M ≈ -(r^2/2) * Γ^r_{t,t} for Schwarzschild-like metric.
     * ∂_r g_{00} ≈ 2M/r^2 → M ≈ (1/2) r^2 * ∂_r g_{00} */
    (void)xi;
    /* Extract from Γ^1_{0,0} ≈ M/r^2 for Schwarzschild at large r */
    double Gamma_r_tt = Gamma[1][0][0];
    double M_approx = radius * radius * Gamma_r_tt;

    (void)m;
    return M_approx;
}
