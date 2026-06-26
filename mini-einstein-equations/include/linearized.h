/**
 * @file linearized.h
 * @brief Linearized gravity: metric perturbations, gravitational waves.
 *
 * Reference: Wald (1984), Ch.4.4; Carroll (2004), Ch.7
 *            MIT 8.962 General Relativity
 *            Maggiore "Gravitational Waves" Vol.1 (2008)
 *
 * Knowledge map:
 *   L4: Linearized Einstein equations: □ hbar_{mu,nu} = -16pi G / c^4 * T_{mu,nu}
 *   L4: Lorenz gauge condition: ∂^{mu} hbar_{mu,nu} = 0
 *   L6: Gravitational wave solutions (TT gauge)
 *   L7: Gravitational wave strain from binary systems
 */

#ifndef MINI_EINSTEIN_LINEARIZED_H
#define MINI_EINSTEIN_LINEARIZED_H

#include "einstein.h"

/* ---------------------------------------------------------------------------
 * L1: Metric perturbation
 * -------------------------------------------------------------------------*/

/**
 * Metric perturbation around Minkowski background:
 *   g_{mu,nu} = eta_{mu,nu} + h_{mu,nu}
 *
 * with |h_{mu,nu}| << 1.
 */
typedef struct {
    Tensor2 h;           /**< perturbation h_{mu,nu} */
    Tensor2 hbar;        /**< trace-reversed: hbar_{mu,nu} = h_{mu,nu} - (1/2) eta_{mu,nu} * h */
    double trace_h;      /**< h = eta^{mu,nu} * h_{mu,nu} */
} MetricPerturbation;

/* ---------------------------------------------------------------------------
 * L4: Linearized Einstein equation
 * -------------------------------------------------------------------------*/

/**
 * Compute trace-reversed perturbation:
 *   hbar_{mu,nu} = h_{mu,nu} - (1/2) * eta_{mu,nu} * h
 */
void compute_hbar(const Tensor2 h, double trace_h, Tensor2 hbar);

/**
 * Compute linearized Einstein tensor in Lorenz gauge:
 *
 *   G^{(1)}_{mu,nu} = -(1/2) * □ hbar_{mu,nu}
 *
 * where □ = eta^{alpha,beta} * ∂_{alpha} ∂_{beta} is the flat-spacetime
 * d'Alembertian.
 *
 * @param dd_hbar  second partial derivatives: dd_hbar[alpha][beta][mu][nu]
 *                 = ∂_{alpha} ∂_{beta} hbar_{mu,nu}
 * @param G1       output: linearized Einstein tensor
 */
void linearized_einstein_tensor(const Tensor4 dd_hbar, Tensor2 G1);

/**
 * Linearized Ricci tensor:
 *   R^{(1)}_{mu,nu} = (1/2) * [∂_{mu} ∂^{alpha} h_{nu,alpha}
 *                             + ∂_{nu} ∂^{alpha} h_{mu,alpha}
 *                             - □ h_{mu,nu}
 *                             - ∂_{mu} ∂_{nu} h]
 *
 * (In Lorenz gauge this simplifies to -(1/2) * □ h_{mu,nu}.)
 *
 * @param h      perturbation
 * @param dh     1st derivatives: dh[alpha][mu][nu] = ∂_{alpha} h_{mu,nu}
 * @param ddh    2nd derivatives: ddh[alpha][beta][mu][nu] = ∂_{alpha} ∂_{beta} h_{mu,nu}
 * @param R1     output: linearized Ricci tensor
 */
void linearized_ricci_tensor(const Tensor2 h, const Tensor3 dh,
                              const Tensor4 ddh, Tensor2 R1);

/* ---------------------------------------------------------------------------
 * L4: Lorenz gauge
 * -------------------------------------------------------------------------*/

/**
 * Lorenz gauge condition:
 *   ∂^{mu} hbar_{mu,nu} = 0
 *
 * @param d_hbar  partial derivatives: d_hbar[mu][mu][nu] = ∂^{mu} hbar_{mu,nu}
 *                (first index = contracted direction)
 * @param div     output divergence 4-vector (should be ≈ 0)
 */
void lorenz_gauge_divergence(const Tensor3 d_hbar, Vector4 div);

/**
 * Gauge transformation of metric perturbation:
 *   h_{mu,nu} -> h_{mu,nu} + ∂_{mu} xi_{nu} + ∂_{nu} xi_{mu}
 *
 * where xi^{mu} is an infinitesimal coordinate transformation.
 *
 * @param h     perturbation (modified in place)
 * @param dxi   derivatives: dxi[mu][nu] = ∂_{mu} xi_{nu}
 */
void gauge_transform_perturbation(Tensor2 h, const Tensor2 dxi);

/* ---------------------------------------------------------------------------
 * L4: Linearized Einstein equation (full form)
 * -------------------------------------------------------------------------*/

/**
 * Evaluate full linearized Einstein equation:
 *   □ hbar_{mu,nu} - ∂_{mu} ∂^{alpha} hbar_{alpha,nu}
 *   - ∂_{nu} ∂^{alpha} hbar_{mu,alpha}
 *   + eta_{mu,nu} * ∂^{alpha} ∂^{beta} hbar_{alpha,beta}
 *   = -16 * pi * G / c^4 * T_{mu,nu}
 *
 * In Lorenz gauge the last three terms vanish.
 *
 * @return max residual
 */
double linearized_einstein_equation(const MetricPerturbation *mp,
                                     const Tensor4 dd_hbar,
                                     const Tensor2 T, double kappa,
                                     Tensor2 residual);

/* ---------------------------------------------------------------------------
 * L6: Gravitational waves — TT gauge
 * -------------------------------------------------------------------------*/

/**
 * Transverse-Traceless (TT) gauge projection.
 *
 * Given a symmetric rank-2 spatial perturbation, project onto
 * the TT subspace:
 *   - Transverse: h^{TT}_{ij} k^{j} = 0 (propagation along z)
 *   - Traceless:  h^{TT}_{ii} = 0
 *
 * For a wave propagating in the z-direction, the TT components are:
 *   h_+   = (h_{xx} - h_{yy}) / 2     (plus polarization)
 *   h_x   =  h_{xy}                    (cross polarization)
 *
 * @param h_spatial  3x3 spatial perturbation matrix
 * @param h_plus     output: plus polarization amplitude
 * @param h_cross    output: cross polarization amplitude
 */
void tt_gauge_project(const double h_spatial[3][3],
                       double *h_plus, double *h_cross);

/**
 * Strain amplitude from a binary system (quadrupole formula):
 *   h ~ (2 * G / (c^4 * D)) * d^2 I_{ij} / dt^2
 *
 * where I_{ij} is the reduced quadrupole moment.
 *
 * @param M_chirp  chirp mass in kg
 * @param f_gw     gravitational wave frequency in Hz
 * @param D        distance to source in meters
 * @return strain amplitude h
 */
double gravitational_wave_strain(double M_chirp, double f_gw, double D);

/**
 * Compute quadrupole radiation power:
 *   P = (G / (5 * c^5)) * sum_{i,j} (d^3 I_{ij}/dt^3)^2
 *
 * @param d3I  array of third time derivatives of reduced quadrupole
 * @return radiated power in Watts
 */
double quadrupole_power(double d3I[3][3]);

/* ---------------------------------------------------------------------------
 * L7: Memory effect
 * -------------------------------------------------------------------------*/

/**
 * Gravitational wave memory effect: net change in metric after wave passes.
 *
 * For a burst of GWs passing an observer:
 *   Delta h_{ij} = h_{ij}(t=+∞) - h_{ij}(t=-∞)
 *
 * Nonlinear (Christodoulou) memory from energy flux of radiated GWs.
 *
 * @param h_initial    initial strain
 * @param h_final      final strain
 * @param Delta_h      output memory (difference)
 */
void gravitational_memory(const double h_initial[3][3],
                           const double h_final[3][3],
                           double Delta_h[3][3]);

#endif /* MINI_EINSTEIN_LINEARIZED_H */
