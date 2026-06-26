/**
 * @file stress_energy.h
 * @brief Stress-energy tensor models (perfect fluid, EM field, scalar field).
 *
 * Reference: Wald (1984), Ch.4; Carroll (2004), Ch.4
 *            MIT 8.962 General Relativity
 *            Schutz "A First Course in GR" Ch.4
 *
 * Knowledge map:
 *   L1: Stress-energy tensor T_{mu,nu}
 *   L4: Energy-momentum conservation ∇_{mu} T^{mu,nu} = 0
 *   L6: Perfect fluid, dust, radiation, EM field
 *   L8: Energy conditions (WEC, SEC, DEC, NEC)
 */

#ifndef MINI_EINSTEIN_STRESS_ENERGY_H
#define MINI_EINSTEIN_STRESS_ENERGY_H

#include "metric.h"

/* ---------------------------------------------------------------------------
 * L1: Energy conditions result
 * -------------------------------------------------------------------------*/

/** Results of checking classical energy conditions */
typedef struct {
    int WEC;    /**< Weak Energy Condition: rho >= 0, rho + p_i >= 0 */
    int NEC;    /**< Null Energy Condition: rho + p_i >= 0 */
    int SEC;    /**< Strong Energy Condition: rho + sum(p_i) >= 0, rho + p_i >= 0 */
    int DEC;    /**< Dominant Energy Condition: rho >= |p_i| */
    double rho; /**< Energy density */
    double p[3]; /**< Principal pressures */
} EnergyConditions;

/* ---------------------------------------------------------------------------
 * L6: Perfect fluid stress-energy
 * -------------------------------------------------------------------------*/

/**
 * Perfect fluid stress-energy tensor:
 *   T_{mu,nu} = (rho + p) * u_{mu} * u_{nu} + p * g_{mu,nu}
 *
 * where u^{mu} is the fluid 4-velocity (u_{mu} u^{mu} = -1).
 *
 * @param m      metric
 * @param rho    energy density (rest frame)
 * @param p      isotropic pressure
 * @param u      fluid 4-velocity (contravariant)
 * @param T      output stress-energy tensor
 */
void stress_energy_perfect_fluid(const Metric *m, double rho, double p,
                                  const Vector4 u, Tensor2 T);

/**
 * Dust (pressureless perfect fluid):
 *   T_{mu,nu} = rho * u_{mu} * u_{nu}
 *
 * Models cold dark matter, non-relativistic matter.
 */
void stress_energy_dust(const Metric *m, double rho,
                         const Vector4 u, Tensor2 T);

/**
 * Radiation stress-energy:
 *   T_{mu,nu} = rho * (u_{mu} * u_{nu} + (1/3) * h_{mu,nu})
 *
 * where h_{mu,nu} = g_{mu,nu} + u_{mu} * u_{nu} is the spatial metric.
 * Equation of state: p = rho / 3.
 */
void stress_energy_radiation(const Metric *m, double rho,
                              const Vector4 u, Tensor2 T);

/* ---------------------------------------------------------------------------
 * L6: Electromagnetic stress-energy
 * -------------------------------------------------------------------------*/

/**
 * Electromagnetic stress-energy tensor:
 *   T_{mu,nu} = F_{mu,alpha} * F_{nu}^{alpha}
 *             - (1/4) * g_{mu,nu} * F_{alpha,beta} * F^{alpha,beta}
 *
 * @param m  metric
 * @param F  Faraday tensor F_{mu,nu} (antisymmetric)
 * @param T  output stress-energy
 */
void stress_energy_electromagnetic(const Metric *m, const Tensor2 F, Tensor2 T);

/* ---------------------------------------------------------------------------
 * L6: Scalar field stress-energy
 * -------------------------------------------------------------------------*/

/**
 * Minimally coupled scalar field stress-energy:
 *   T_{mu,nu} = ∇_{mu} phi * ∇_{nu} phi
 *             - g_{mu,nu} * [(1/2) * ∇^{alpha} phi * ∇_{alpha} phi + V(phi)]
 *
 * @param m      metric
 * @param phi    scalar field value
 * @param dphi   covariant derivative ∇_{mu} phi (typically = ∂_{mu} phi)
 * @param V      potential energy V(phi)
 * @param T      output stress-energy
 */
void stress_energy_scalar_field(const Metric *m, double phi,
                                 const Vector4 dphi, double V, Tensor2 T);

/* ---------------------------------------------------------------------------
 * L6: Cosmological constant "stress-energy"
 * -------------------------------------------------------------------------*/

/**
 * Effective stress-energy of cosmological constant:
 *   T^{Lambda}_{mu,nu} = -(Lambda / kappa) * g_{mu,nu}
 *
 * Often written as Lambda * g_{mu,nu} appearing on the left side of
 * the Einstein equation.
 */
void stress_energy_cosmological_constant(const Metric *m, double Lambda,
                                          double kappa, Tensor2 T);

/* ---------------------------------------------------------------------------
 * L8: Energy conditions
 * -------------------------------------------------------------------------*/

/**
 * Check all four classical energy conditions given energy density and
 * principal pressures (in the rest frame of the matter).
 *
 * - WEC: rho >= 0 and rho + p_i >= 0 for all i
 * - NEC: rho + p_i >= 0 for all i
 * - SEC: rho + sum(p_i) >= 0 and rho + p_i >= 0 for all i
 * - DEC: rho >= |p_i| for all i
 *
 * @param rho  energy density
 * @param p    array of 3 principal pressures
 * @param ec   output energy conditions result
 */
void energy_conditions_check(double rho, const double p[3],
                              EnergyConditions *ec);

/**
 * Check energy conditions from a general stress-energy tensor
 * by diagonalizing it to find principal pressures.
 *
 * @param m   metric
 * @param T   stress-energy tensor T_{mu,nu}
 * @param u   local observer 4-velocity
 * @param ec  output energy conditions
 */
void energy_conditions_from_tensor(const Metric *m, const Tensor2 T,
                                    const Vector4 u, EnergyConditions *ec);

/* ---------------------------------------------------------------------------
 * L4: Stress-energy conservation
 * -------------------------------------------------------------------------*/

/**
 * Check stress-energy conservation:
 *   ∇_{mu} T^{mu,nu} = 0
 *
 * @param m      metric
 * @param Gamma  Christoffel symbols
 * @param T      stress-energy tensor
 * @param dT     partial derivatives: dT[sigma][mu][nu] = ∂_{sigma} T^{mu,nu}
 * @param divT   output: ∇_{mu} T^{mu,nu} (4-vector)
 */
void stress_energy_divergence(const Metric *m, const Tensor3 Gamma,
                               const Tensor2 T, const Tensor3 dT,
                               Vector4 divT);

#endif /* MINI_EINSTEIN_STRESS_ENERGY_H */
