/**
 * @file linearized.c
 * @brief Linearized gravity: perturbations, gravitational waves.
 *
 * Reference: Wald (1984), Ch.4.4; Carroll (2004), Ch.7
 *            Maggiore "Gravitational Waves" Vol.1 (2008)
 *            Misner, Thorne, Wheeler "Gravitation" Ch.35-37
 *
 * Linearized gravity is the weak-field approximation:
 *   g_{μ,ν} = η_{μ,ν} + h_{μ,ν},   |h_{μ,ν}| ≪ 1
 */

#include "linearized.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ========================================================================
 * L4: Metric perturbation basics
 * ========================================================================*/

/**
 * Compute trace-reversed perturbation:
 *   h̅_{μ,ν} = h_{μ,ν} - (1/2) η_{μ,ν} h
 *   where h = η^{μ,ν} h_{μ,ν}
 *
 * The trace reversal simplifies the linearized Einstein equation
 * dramatically when combined with the Lorenz gauge condition.
 *
 * Knowledge point: The trace-reversed perturbation is the key variable
 * that decouples the linearized Einstein equations into wave equations.
 * Without it, the equations are coupled by the trace of h.
 */
void compute_hbar(const Tensor2 h, double trace_h, Tensor2 hbar)
{
    double half_trace = 0.5 * trace_h;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            hbar[mu][nu] = h[mu][nu]
                         - half_trace * MINKOWSKI_METRIC[mu][nu];
        }
    }
}

/**
 * Linearized Einstein tensor in Lorenz gauge:
 *   G^{(1)}_{μ,ν} = -(1/2) □ h̅_{μ,ν}
 *
 * where □ = η^{α,β} ∂_α ∂_β is the flat d'Alembertian.
 *
 * The full expression (before gauge fixing) is:
 *   G^{(1)}_{μ,ν} = -(1/2)[□ h̅_{μ,ν} + η_{μ,ν} ∂^α ∂^β h̅_{α,β}
 *                          - ∂_μ ∂^α h̅_{α,ν} - ∂_ν ∂^α h̅_{μ,α}]
 *
 * In Lorenz gauge (∂^μ h̅_{μ,ν} = 0), only the □ term survives.
 *
 * Knowledge point: In Lorenz gauge, each component of h̅ satisfies a
 * simple wave equation with the stress-energy as source.
 */
void linearized_einstein_tensor(const Tensor4 dd_hbar, Tensor2 G1)
{
    /* □ h̅_{μ,ν} = Σ_{α,β} η^{α,β} ∂_α ∂_β h̅_{μ,ν} */
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double box_hbar = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    box_hbar += MINKOWSKI_METRIC[alpha][beta]
                              * dd_hbar[alpha][beta][mu][nu];
                }
            }
            G1[mu][nu] = -0.5 * box_hbar;
        }
    }
}

/**
 * Linearized Ricci tensor (full expression, before gauge fixing):
 *
 *   R^{(1)}_{μ,ν} = (1/2)[∂_μ ∂^α h_{ν,α} + ∂_ν ∂^α h_{μ,α}
 *                         - □ h_{μ,ν} - ∂_μ ∂_ν h]
 *
 * After Lorenz gauge fixing: R^{(1)}_{μ,ν} = -(1/2) □ h_{μ,ν}
 *
 * Knowledge point: The linearized Ricci tensor captures the first-order
 * curvature induced by a weak gravitational field. For a Newtonian source
 * (static, weak field): R^{(1)}_{00} ≈ ∇^2 φ where φ = -GM/r is the
 * Newtonian potential — this connects GR to Newtonian gravity.
 */
void linearized_ricci_tensor(const Tensor2 h, const Tensor3 dh,
                              const Tensor4 ddh, Tensor2 R1)
{
    (void)dh; /* dh is used indirectly through ddh (second derivatives) */
    /* Compute trace h = η^{μ,ν} h_{μ,ν} */
    double trace_h = 0.0;
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            trace_h += MINKOWSKI_METRIC[a][b] * h[a][b];
        }
    }

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            /* □ h_{μ,ν} */
            double box_h = 0.0;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    box_h += MINKOWSKI_METRIC[a][b]
                           * ddh[a][b][mu][nu];
                }
            }

            /* ∂_μ ∂^α h_{ν,α} */
            double term_mu = 0.0;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    term_mu += MINKOWSKI_METRIC[a][b]
                             * ddh[mu][a][nu][b];
                }
            }

            /* ∂_ν ∂^α h_{μ,α} */
            double term_nu = 0.0;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    term_nu += MINKOWSKI_METRIC[a][b]
                             * ddh[nu][a][mu][b];
                }
            }

            /* ∂_μ ∂_ν h */
            double term_trace = ddh[mu][nu][0][0] * MINKOWSKI_METRIC[0][0];
            for (int a = 1; a < 4; a++) {
                term_trace += ddh[mu][nu][a][a]; /* η^{aa}=1 for a=1,2,3 */
            }

            R1[mu][nu] = 0.5 * (term_mu + term_nu - box_h - term_trace);
        }
    }
}

/* ========================================================================
 * L4: Lorenz gauge
 * ========================================================================*/

/**
 * Lorenz gauge divergence:
 *   ∂^μ h̅_{μ,ν} = 0
 *
 * If this holds, the linearized Einstein equations reduce to simple
 * wave equations. The gauge condition is the GR analog of the
 * Lorenz gauge in electromagnetism: ∂^μ A_μ = 0.
 *
 * Knowledge point: The Lorenz gauge decouples the equations, making
 * gravitational waves mathematically analogous to EM waves but with
 * spin-2 rather than spin-1.
 */
void lorenz_gauge_divergence(const Tensor3 d_hbar, Vector4 div)
{
    for (int nu = 0; nu < 4; nu++) {
        double sum = 0.0;
        for (int mu = 0; mu < 4; mu++) {
            for (int alpha = 0; alpha < 4; alpha++) {
                /* ∂^μ h̅_{μ,ν} = η^{μ,α} ∂_α h̅_{μ,ν} */
                sum += MINKOWSKI_METRIC[mu][alpha]
                     * d_hbar[alpha][mu][nu];
            }
        }
        div[nu] = sum;
    }
}

/**
 * Gauge transformation of perturbation:
 *   h_{μ,ν} → h_{μ,ν} + ∂_μ ξ_ν + ∂_ν ξ_μ
 *
 * where ξ^μ is an infinitesimal vector field.
 *
 * This is the linearized version of the diffeomorphism:
 *   x^μ → x^μ + ξ^μ
 *
 * Physical quantities (like the Riemann tensor) are gauge-invariant.
 *
 * Knowledge point: Gauge freedom in linearized gravity corresponds to
 * infinitesimal coordinate transformations. Physical observables must
 * be gauge-invariant. This is analogous to the gauge freedom in EM:
 * A_μ → A_μ + ∂_μ λ.
 */
void gauge_transform_perturbation(Tensor2 h, const Tensor2 dxi)
{
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            h[mu][nu] += dxi[mu][nu] + dxi[nu][mu];
        }
    }
}

/* ========================================================================
 * L4: Full linearized Einstein equation
 * ========================================================================*/

double linearized_einstein_equation(const MetricPerturbation *mp,
                                     const Tensor4 dd_hbar,
                                     const Tensor2 T, double kappa,
                                     Tensor2 residual)
{
    Tensor2 G1;
    linearized_einstein_tensor(dd_hbar, G1);

    double max_res = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            /* G^{(1)}_{μ,ν} = κ T_{μ,ν} */
            residual[mu][nu] = G1[mu][nu] - kappa * T[mu][nu];
            double r = fabs(residual[mu][nu]);
            if (r > max_res) max_res = r;
        }
    }
    (void)mp;
    return max_res;
}

/* ========================================================================
 * L6: TT gauge and gravitational waves
 * ========================================================================*/

/**
 * TT (Transverse-Traceless) gauge projection.
 *
 * For a wave propagating in the z-direction, the spatial perturbation
 * is reduced to two independent components:
 *   h_+   = (h_{xx} - h_{yy}) / 2   (plus polarization)
 *   h_×   =  h_{xy}                 (cross polarization)
 *
 * Properties of TT gauge:
 *   - h^TT_{0μ} = 0 (no temporal components)
 *   - h^TT_{zi} = 0 (transverse to propagation)
 *   - Σ_i h^TT_{ii} = 0 (traceless)
 *
 * Knowledge point: The TT gauge isolates the two physical degrees of
 * freedom of gravitational waves. This is the gauge used by GW detectors
 * (LIGO, Virgo, KAGRA) to express strain measurements.
 *
 * LIGO's first detection (GW150914) measured a peak strain of h~10^{-21}.
 */
void tt_gauge_project(const double h_spatial[3][3],
                       double *h_plus, double *h_cross)
{
    *h_plus  = 0.5 * (h_spatial[0][0] - h_spatial[1][1]);
    *h_cross = h_spatial[0][1]; /* h_{xy} */
}

/**
 * Gravitational wave strain from quadrupole formula.
 *
 * Amplitude:
 *   h ≈ (2G / (c^4 D)) * d²I/dt²
 *
 * For a binary system:
 *   h ≈ (4G / (c^4 D)) * μ * ω² * r²
 *   (where μ = reduced mass, ω = orbital frequency, r = separation)
 *
 * With chirp mass M_c = (m1*m2)^{3/5} / (m1+m2)^{1/5},
 * at frequency f_gw:
 *   h ≈ (4G^{5/3} / (c^4 D)) * (π f_gw)^{2/3} * M_c^{5/3}
 *
 * Knowledge point: The quadrupole formula (Einstein 1918) is the
 * foundation of gravitational wave astronomy. It shows that GW strain
 * falls as 1/D (unlike EM which falls as 1/D²), enabling detection
 * of sources at cosmological distances.
 */
double gravitational_wave_strain(double M_chirp, double f_gw, double D)
{
    double pi_f = M_PI * f_gw;
    double factor = 4.0 * pow(G_NEWTON, 5.0/3.0)
                    / (C_LIGHT * C_LIGHT * C_LIGHT * C_LIGHT * D);
    double M_c_53 = pow(M_chirp, 5.0/3.0);
    double f_23 = pow(pi_f, 2.0/3.0);

    return factor * f_23 * M_c_53;
}

/**
 * Quadrupole radiation power:
 *   P = (G / (5 c^5)) Σ_{i,j} ⟨d³I_{ij}/dt³ d³I_{ij}/dt³⟩
 *
 * For a binary:
 *   P = (32/5) * (G^4 / c^5) * (μ² M³ / r^5)
 *
 * where μ = m1 m2 / (m1+m2) [reduced mass], M = m1+m2 [total mass].
 *
 * Knowledge point: The quadrupole power formula quantifies the energy
 * carried away by gravitational waves. This energy loss causes binary
 * orbits to decay — the Hulse-Taylor pulsar (PSR B1913+16) provided
 * the first indirect evidence for GWs through orbital decay matching
 * the GR prediction (Nobel Prize 1993).
 */
double quadrupole_power(double d3I[3][3])
{
    double sum_sq = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sum_sq += d3I[i][j] * d3I[i][j];
        }
    }
    return G_NEWTON / (5.0 * pow(C_LIGHT, 5.0)) * sum_sq;
}

/* ========================================================================
 * L7: Gravitational memory effect
 * ========================================================================*/

/**
 * Gravitational wave memory: net strain offset after wave passage.
 *
 * Two types:
 *   Linear memory: from unbound motion of sources (e.g., hyperbolic flyby)
 *   Nonlinear (Christodoulou) memory: from energy flux of GWs themselves
 *
 * The memory Δh_{ij} is the difference between initial and final
 * asymptotic metric states.
 *
 * Knowledge point: The memory effect is a permanent change in the
 * relative separation of test masses after a GW passes. It is a
 * prediction of the nonlinearity of GR (hereditary effect).
 * Detection would be a profound test of GR's nonlinear structure.
 */
void gravitational_memory(const double h_initial[3][3],
                           const double h_final[3][3],
                           double Delta_h[3][3])
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Delta_h[i][j] = h_final[i][j] - h_initial[i][j];
        }
    }
}
