/**
 * @file stress_energy.c
 * @brief Stress-energy tensor models and energy conditions.
 *
 * Reference: Wald (1984), Ch.4; Carroll (2004), Ch.4
 *            Schutz "A First Course in GR" Ch.4
 *            Hawking & Ellis "The Large Scale Structure of Space-Time" (1973)
 *
 * The stress-energy tensor T_{μ,ν} is the source of spacetime curvature.
 * Its components represent energy density, momentum density, and stress.
 */

#include "stress_energy.h"
#include <math.h>
#include <string.h>

/* ========================================================================
 * L6: Perfect fluid stress-energy
 * ========================================================================*/

/**
 * Perfect fluid: T_{μ,ν} = (ρ + p) u_μ u_ν + p g_{μ,ν}
 *
 * where u^μ is the 4-velocity (u_μ u^μ = -1 in our convention).
 * u_μ = g_{μ,ν} u^ν.
 *
 * Knowledge point: The perfect fluid is the simplest matter model in GR.
 * It describes any fluid with no viscosity or heat conduction. The
 * equation of state p = w ρ determines the fluid type:
 *   w = 0    : dust (matter-dominated era)
 *   w = 1/3  : radiation
 *   w = -1   : cosmological constant (dark energy)
 */
void stress_energy_perfect_fluid(const Metric *m, double rho, double p,
                                  const Vector4 u, Tensor2 T)
{
    /* Compute u_μ = g_{μ,ν} u^ν */
    Vector4 u_dn;
    lower_vector(m->g, u, u_dn);

    /* Outer product u_μ u_ν */
    Tensor2 uu;
    outer_product_vector4(u_dn, u_dn, uu);

    /* T_{μ,ν} = (ρ + p) u_μ u_ν + p g_{μ,ν} */
    double rho_plus_p = rho + p;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] = rho_plus_p * uu[mu][nu] + p * m->g[mu][nu];
        }
    }
}

/**
 * Dust (pressureless fluid): T_{μ,ν} = ρ u_μ u_ν
 *
 * Knowledge point: Dust models non-relativistic matter (cold dark matter,
 * galaxies as point particles). The vanishing pressure means the matter
 * is "pressureless" — no internal random motions.
 *
 * For dust, the energy conditions simplify to ρ ≥ 0.
 */
void stress_energy_dust(const Metric *m, double rho,
                         const Vector4 u, Tensor2 T)
{
    /* u_μ = g_{μ,ν} u^ν */
    Vector4 u_dn;
    lower_vector(m->g, u, u_dn);

    outer_product_vector4(u_dn, u_dn, T);
    scale_tensor2(T, rho);
}

/**
 * Radiation: T_{μ,ν} = ρ [u_μ u_ν + (1/3) h_{μ,ν}]
 *
 * where h_{μ,ν} = g_{μ,ν} + u_μ u_ν is the spatial projection tensor.
 *
 * Knowledge point: Radiation has p = ρ/3 (traced from the Maxwell
 * stress-energy tensor being trace-free: T^μ_μ = 0).
 *
 * Radiation dominated the early universe. The CMB temperature
 * scales as T ∝ 1/a during this era.
 */
void stress_energy_radiation(const Metric *m, double rho,
                              const Vector4 u, Tensor2 T)
{
    Vector4 u_dn;
    lower_vector(m->g, u, u_dn);

    Tensor2 uu;
    outer_product_vector4(u_dn, u_dn, uu);

    /* h_{μ,ν} = g_{μ,ν} + u_μ u_ν (spatial projection) */
    /* T_{μ,ν} = ρ u_μ u_ν + (ρ/3) h_{μ,ν} = (4ρ/3) u_μ u_ν + (ρ/3) g_{μ,ν} */
    double rho_over_3 = rho / 3.0;
    double four_thirds_rho = (4.0 / 3.0) * rho;

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] = four_thirds_rho * uu[mu][nu]
                      + rho_over_3 * m->g[mu][nu];
        }
    }
}

/* ========================================================================
 * L6: Electromagnetic stress-energy
 * ========================================================================*/

/**
 * Electromagnetic stress-energy:
 *   T_{μ,ν} = F_{μ,α} F_ν^α - (1/4) g_{μ,ν} F_{α,β} F^{α,β}
 *
 * where F^{α,β} = g^{α,γ} g^{β,δ} F_{γ,δ}.
 *
 * Properties:
 *   - T is symmetric (F asymmetric, but the contraction yields symmetry)
 *   - T^μ_μ = 0 (trace-free — conformal invariance of EM)
 *   - ∇^μ T_{μ,ν} = 0 follows from Maxwell equations ∇^μ F_{μ,ν} = 0
 *
 * Knowledge point: The EM stress-energy is the source in the
 * electrovacuum Einstein equations (e.g., Reissner-Nordström black hole).
 */
void stress_energy_electromagnetic(const Metric *m, const Tensor2 F, Tensor2 T)
{
    /* Compute F^{α,β} = g^{α,γ} g^{β,δ} F_{γ,δ} */
    Tensor2 F_up;
    for (int alpha = 0; alpha < 4; alpha++) {
        for (int beta = 0; beta < 4; beta++) {
            double sum = 0.0;
            for (int gamma = 0; gamma < 4; gamma++) {
                for (int delta = 0; delta < 4; delta++) {
                    sum += m->g_inv[alpha][gamma]
                         * m->g_inv[beta][delta]
                         * F[gamma][delta];
                }
            }
            F_up[alpha][beta] = sum;
        }
    }

    /* Compute F_{α,β} F^{α,β} (scalar invariant) */
    double F2 = double_contraction(F_up, F);

    /* T_{μ,ν} = F_{μ,α} F_ν^α - (1/4) g_{μ,ν} F^2 */
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            /* F_{μ,α} F_ν^α = F_{μ,α} g^{α,β} F_{β,ν} ... wait, no.
             * F_ν^α = g^{α,β} F_{β,ν} is wrong; F_ν^α means the first index
             * is upper: F^α_ν = g^{α,β} F_{β,ν}.
             *
             * But F_{μ,α} F_ν^α = F_{μ,α} * (g^{α,β} F_{ν,β})? No.
             *
             * F_ν^α means: one index raised: F^α_ν.
             * But the expression is F_{μ,α} F_ν^α where the α on F_{μ,α}
             * is lower and on F_ν^α the α is upper.
             * So: F_ν^α = g^{α,β} F_{ν,β}? No, that would raise the second index.
             *
             * Convention: F_ν^α = g^{α,β} F_{β,ν}.
             * This means the first index of F is raised and then swapped.
             *
             * So: F_{μ,α} F_ν^α = F_{μ,α} * g^{α,β} * F_{β,ν}
             */

            double Fmu_nu = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    Fmu_nu += F[mu][alpha] * m->g_inv[alpha][beta] * F[beta][nu];
                }
            }

            T[mu][nu] = Fmu_nu - 0.25 * F2 * m->g[mu][nu];
        }
    }
}

/* ========================================================================
 * L6: Scalar field stress-energy
 * ========================================================================*/

/**
 * Minimally coupled scalar field:
 *   T_{μ,ν} = ∂_μ φ ∂_ν φ - g_{μ,ν} [(1/2) g^{α,β} ∂_α φ ∂_β φ + V(φ)]
 *
 * The equation of motion for φ is the Klein-Gordon equation:
 *   □ φ - dV/dφ = 0
 *
 * Knowledge point: Scalar fields are fundamental in cosmology (inflaton,
 * quintessence) and appear in many modified gravity theories.
 */
void stress_energy_scalar_field(const Metric *m, double phi,
                                 const Vector4 dphi, double V, Tensor2 T)
{
    (void)phi; /* phi not directly used, only dphi and potential */

    /* Kinetic term: X = (1/2) g^{α,β} ∂_α φ ∂_β φ */
    double X = 0.0;
    for (int alpha = 0; alpha < 4; alpha++) {
        for (int beta = 0; beta < 4; beta++) {
            X += m->g_inv[alpha][beta] * dphi[alpha] * dphi[beta];
        }
    }
    X *= 0.5;

    /* T_{μ,ν} = ∂_μ φ ∂_ν φ - g_{μ,ν} (X + V) */
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] = dphi[mu] * dphi[nu]
                      - m->g[mu][nu] * (X + V);
        }
    }
}

/* ========================================================================
 * L6: Cosmological constant as effective stress-energy
 * ========================================================================*/

void stress_energy_cosmological_constant(const Metric *m, double Lambda,
                                          double kappa, Tensor2 T)
{
    double factor = -Lambda / kappa;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            T[mu][nu] = factor * m->g[mu][nu];
        }
    }
}

/* ========================================================================
 * L8: Energy conditions
 * ========================================================================*/

/**
 * Check classical energy conditions.
 *
 * Four standard energy conditions (Hawking & Ellis 1973):
 *
 * WEC (Weak Energy Condition):
 *   ρ ≥ 0,  ρ + p_i ≥ 0 for all i
 *   Physical meaning: energy density is non-negative for all observers.
 *   Satisfied by all known classical matter.
 *
 * NEC (Null Energy Condition):
 *   ρ + p_i ≥ 0 for all i
 *   Physical meaning: light rays are focused by matter (always attractive).
 *   Required for the Penrose singularity theorem.
 *
 * SEC (Strong Energy Condition):
 *   ρ + Σ p_i ≥ 0,  ρ + p_i ≥ 0 for all i
 *   Physical meaning: gravity is always attractive.
 *   Violated by dark energy (cosmological constant), causing acceleration.
 *
 * DEC (Dominant Energy Condition):
 *   ρ ≥ |p_i| for all i
 *   Physical meaning: energy flux is timelike or null (no superluminal flow).
 *   Equivalent to causality of energy transport.
 *
 * Knowledge point: Energy conditions constrain possible matter types.
 * Dark energy violates SEC (explaining accelerating expansion).
 * Quantum effects (Casimir) can violate all conditions.
 */
void energy_conditions_check(double rho, const double p[3],
                              EnergyConditions *ec)
{
    ec->rho = rho;
    ec->p[0] = p[0]; ec->p[1] = p[1]; ec->p[2] = p[2];

    double sum_p = p[0] + p[1] + p[2];

    /* WEC */
    ec->WEC = (rho >= 0.0)
              && (rho + p[0] >= 0.0)
              && (rho + p[1] >= 0.0)
              && (rho + p[2] >= 0.0);

    /* NEC */
    ec->NEC = (rho + p[0] >= 0.0)
              && (rho + p[1] >= 0.0)
              && (rho + p[2] >= 0.0);

    /* SEC */
    ec->SEC = (rho + sum_p >= 0.0)
              && (rho + p[0] >= 0.0)
              && (rho + p[1] >= 0.0)
              && (rho + p[2] >= 0.0);

    /* DEC */
    ec->DEC = (rho >= fabs(p[0]))
              && (rho >= fabs(p[1]))
              && (rho >= fabs(p[2]));
}

void energy_conditions_from_tensor(const Metric *m, const Tensor2 T,
                                    const Vector4 u, EnergyConditions *ec)
{
    /* Extract energy density as measured by observer u:
     * ρ = T_{μ,ν} u^μ u^ν */
    Vector4 u_dn;
    lower_vector(m->g, u, u_dn);

    double rho = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            rho += T[mu][nu] * u[mu] * u[nu];
        }
    }

    /* Spatial stress: p_i = T_{i,i} (in orthonormal frame, approximately) */
    /* For simplicity, use diagonal components of T in coordinate basis
     * and assume near-orthonormal coordinates. A full computation would
     * project T onto a spatial triad orthogonal to u. */
    double p[3] = {T[1][1], T[2][2], T[3][3]};
    energy_conditions_check(rho, p, ec);
}

/* ========================================================================
 * L4: Stress-energy conservation
 * ========================================================================*/

/**
 * Covariant divergence of stress-energy tensor:
 *   ∇_μ T^{μ,ν} = ∂_μ T^{μ,ν} + Γ^μ_{μ,σ} T^{σ,ν} + Γ^ν_{μ,σ} T^{μ,σ}
 *
 * Energy-momentum conservation: ∇_μ T^{μ,ν} = 0.
 *
 * Knowledge point: This is the GR generalization of the conservation
 * laws ∂_μ T^{μ,ν} = 0 from special relativity. It follows identically
 * from the contracted Bianchi identity and the Einstein equation.
 */
void stress_energy_divergence(const Metric *m, const Tensor3 Gamma,
                               const Tensor2 T, const Tensor3 dT,
                               Vector4 divT)
{
    /* Compute T^{μ,ν} by raising both indices */
    Tensor2 T_up;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    sum += m->g_inv[mu][alpha]
                         * m->g_inv[nu][beta]
                         * T[alpha][beta];
                }
            }
            T_up[mu][nu] = sum;
        }
    }

    /* dT[sigma][mu][nu] is ∂_sigma T^{mu,nu} */
    for (int nu = 0; nu < 4; nu++) {
        double div = 0.0;
        for (int mu = 0; mu < 4; mu++) {
            /* Partial derivative term */
            div += dT[mu][mu][nu];

            /* Connection terms */
            for (int sigma = 0; sigma < 4; sigma++) {
                div += Gamma[mu][mu][sigma] * T_up[sigma][nu]
                     + Gamma[nu][mu][sigma] * T_up[mu][sigma];
            }
        }
        divT[nu] = div;
    }
}
