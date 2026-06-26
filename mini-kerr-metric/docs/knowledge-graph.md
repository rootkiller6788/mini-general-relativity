# Knowledge Graph — mini-kerr-metric

## L1 — Definitions (Complete ✅)

| Item | Definition | Location |
|------|-----------|----------|
| Kerr metric g_μν | ds² = -(1-2Mr/Σ)dt² - (4Mar sin²θ/Σ)dtdφ + (Σ/Δ)dr² + Σdθ² + ...dφ² | kerr_metric.h, kerr_metric.c |
| Sigma Σ | r² + a² cos²θ | kerr_metric.c:kerr_sigma |
| Delta Δ | r² - 2Mr + a² | kerr_metric.c:kerr_delta |
| A-function | (r²+a²)² - a²Δ sin²θ | kerr_metric.c:kerr_A_function |
| ω-tilde | 2Mar/Σ | kerr_metric.c:kerr_omega_tilde |
| Horizon radii r_± | M ± √(M² - a²) | kerr_horizons.c:kerr_horizon_radii |
| Ergosphere boundary | r_ergo = M + √(M² - a² cos²θ) | kerr_horizons.c:kerr_ergosphere_radius |
| Ring singularity | r=0, θ=π/2, ring of radius a | kerr_horizons.c:kerr_ring_singularity_data |
| Frame-dragging ω | -g_tφ/g_φφ | kerr_metric.c:kerr_frame_dragging_omega |
| Carter constant Q | From Killing tensor K_μν | kerr_geodesics.c:kerr_carter_constant |
| ISCO radius | Innermost stable circular orbit | kerr_geodesics.c:kerr_isco |
| Petrov Type D | Algebraically special classification | kerr_metric.h:PetrovType |
| Weyl scalars Ψ₁ | Ψ₂ = -M/(r-ia cosθ)³ | kerr_metric.c:kerr_weyl_scalars |

## L2 — Core Concepts (Complete ✅)

| Concept | Description | Implementation |
|---------|-------------|----------------|
| Frame dragging | Rotation of inertial frames by spinning mass | kerr_frame_dragging_omega |
| Ergoregion | Region where no static observers can exist | kerr_is_inside_ergosphere |
| Penrose process | Rotational energy extraction via particle decay | kerr_penrose_process |
| Superradiance | Wave amplification by rotating BH | kerr_superradiance_condition |
| Blandford-Znajek | EM energy extraction via B-fields | kerr_blandford_znajek |
| Black hole thermodynamics | Entropy, temperature, area theorem | kerr_horizons.c |
| Cosmic censorship | a ≤ M for physical BH | kerr_extremality_parameter |
| No-hair theorem | Kerr uniquely characterized by M and a | kerr_bh_classify |
| Killing vectors | ∂_t (stationary), ∂_φ (axisymmetric) | kerr_metric_bl |
| Killing tensor | Quadratic symmetry for Carter constant | kerr_carter_constant |

## L3 — Mathematical Structures (Complete ✅)

| Structure | Description | Location |
|-----------|-------------|----------|
| Boyer-Lindquist coordinates | Standard Kerr coordinates (singular at Δ=0) | kerr_metric_bl |
| Kerr-Schild coordinates | Cartesian, regular at horizon | kerr_metric_ks |
| Christoffel symbols Γ^μ_νρ | 20 non-zero components | kerr_christoffel_bl |
| Riemann tensor R^μ_νρσ | 20 independent components (vacuum = Weyl) | kerr_riemann_bl |
| Kretschmann scalar K | Curvature invariant (→ ∞ at ring sing.) | kerr_kretschmann |
| Weyl scalars Ψ₀..Ψ₄ | NP formalism (only Ψ₂ ≠ 0 for Type D) | kerr_weyl_scalars |
| Spin coefficients | 12 complex NP spin coefficients | kerr_spin_coefficients |
| Kinnersley tetrad | Null tetrad aligned with principal null dirs | kerr_kinnersley_tetrad |
| ADM 3+1 decomposition | Lapse, shift, 3-metric, extrinsic curvature | kerr_metric.c |
| ZAMO tetrad | Zero Angular Momentum Observer frame | kerr_zamo_tetrad |

## L4 — Fundamental Laws (Complete ✅)

| Law/Theorem | Statement | Verification |
|-------------|-----------|--------------|
| Einstein vacuum eqns | R_μν = 0 (Kerr is exact solution) | kerr_ricci_bl |
| Smarr formula | M = 2Ω_H J + (κ/4π)A | kerr_smarr_formula_residual |
| First law BH thermo | dM = (κ/8π)dA + Ω_H dJ | kerr_smarr_formula_check |
| Area theorem | dA ≥ 0 | kerr_area_increase_ratio |
| Bekenstein-Hawking entropy | S = A/4 | kerr_entropy_geometric |
| Hawking temperature | T = κ/(2π) | kerr_temperature_geometric |
| Lense-Thirring precession | Ω_LT = 2J/r³ | kerr_lense_thirring_frequency |
| Geodesic equation | Carter separation of Hamilton-Jacobi | kerr_geodesic_rk4_step |

## L5 — Computational Methods (Complete ✅)

| Method | Application | Implementation |
|--------|-------------|----------------|
| Brent's method | Root-finding for horizons, turning points | kerr_brent_root |
| Newton-Raphson | Gradient-based root finding | kerr_newton_raphson |
| Cubic solver | Cardano/trigonometric | kerr_cubic_real_roots |
| Quartic solver | Ferrari's method | kerr_quartic_real_roots |
| RK4 fixed-step | Geodesic integration | kerr_rk4_fixed_step |
| RK5(4) adaptive | High-precision ODE integration | kerr_rkf45_step |
| Adams-Bashforth-Moulton | Multi-step predictor-corrector | kerr_adams_bashforth_moulton |
| Chebyshev differentiation | Spectral methods | kerr_chebyshev_diff_matrix |
| Barycentric interpolation | Spectral function evaluation | kerr_chebyshev_interpolate |
| Numerical differentiation | Richardson extrapolation | kerr_numerical_metric_derivative |
| Ray tracing | BH imaging | kerr_raytrace_single |

## L6 — Canonical Systems (Complete ✅)

| System | Properties | Implementation |
|--------|-----------|----------------|
| Schwarzschild limit (a=0) | r_isco=6M, r_ph=3M, kappa=1/(4M) | Tests verify |
| Extremal Kerr (a=M) | r_+=M, kappa=0, eta=29.3% | Tests verify |
| Retrograde extremal | r_isco=9M, r_ph=4M | Tests verify |
| Photon sphere | r_ph=3M (Schw), M<r_ph<4M (Kerr) | kerr_photon_orbit |
| ISCO vs spin | 1M < r_isco < 9M | kerr_isco_radius_bardeen |
| Bound orbits | Radial turning points | kerr_radial_turning_points |

## L7 — Applications (Complete ✅)

| Application | Description | Implementation |
|-------------|-------------|----------------|
| Penrose energy extraction | Up to 29.3% efficiency | kerr_penrose_max_efficiency |
| Superradiant amplification | Scalar wave gain | kerr_superradiance_amplification_scalar |
| Blandford-Znajek jets | AGN jet power ~ 10^45 erg/s | kerr_bz_power_physical |
| Accretion disk (Novikov-Thorne) | Thin disk model | kerr_accretion_disk |
| Iron K-alpha line profile | Relativistic line broadening | kerr_iron_line_profile |
| Black hole shadow | EHT observables | kerr_shadow_observables |
| Gravitational wave ringdown | QNM frequencies (Berti+2006) | kerr_qnm_frequency_real |
| Disk spectrum | Multi-color blackbody | kerr_disk_spectrum |

## L8 — Advanced Topics (Partial ⚠️)

| Topic | Status | Implementation |
|-------|--------|----------------|
| Strong gravitational lensing | Partial | kerr_strong_lensing_images |
| Geodesic deviation | Partial | kerr_geodesic_deviation |
| Tidal tensor | Partial | kerr_tidal_tensor |
| Teukolsky equation (spectral) | Partial | kerr_teukolsky_radial_spectral |
| Collisional Penrose process | Partial | kerr_collisional_cm_energy |
| Superradiant instability | Partial | kerr_superradiant_instability_timescale |
| Closed timelike curves | Partial | kerr_has_ctc_region |
| Spin coefficients | Complete | kerr_spin_coefficients |

## L9 — Research Frontiers (Partial ⚠️)

| Topic | Status | Notes |
|-------|--------|-------|
| Quantum gravity at ring singularity | Documented | K → ∞ at the ring, requires quantum gravity |
| Information loss paradox (Kerr) | Documented | Cauchy horizon instability |
| Black hole spectroscopy | Partial | QNM fitting formulas implemented |
| Kerr in modified gravity | Not implemented | Kerr-Sen, Kerr-Taub-NUT variants |
| Exotic compact objects | Not implemented | Boson stars, gravastars |
| Holographic principle for Kerr | Documented | S ∼ A/4 |
