# Knowledge Graph — mini-einstein-equations

Nine-level knowledge coverage map for Einstein Field Equations module.

## L1: Definitions — Complete ✅

| # | Concept | C Implementation | Lean Formalization |
|---|---------|-----------------|-------------------|
| 1 | Scalar field | `typedef double Scalar` | `def ScalarField` |
| 2 | Vector field (4-vector) | `typedef double Vector4[4]` | `def VectorField` |
| 3 | Rank-2 tensor | `typedef double Tensor2[4][4]` | `def Tensor2Field` |
| 4 | Rank-3 tensor (Christoffel) | `typedef double Tensor3[4][4][4]` | `def ChristoffelSymbols` |
| 5 | Rank-4 tensor (Riemann) | `typedef double Tensor4[4][4][4][4]` | `def RiemannTensor` |
| 6 | Metric g_{μ,ν} | `struct Metric` | `structure LorentzianMetric` |
| 7 | Minkowski metric η_{μ,ν} | `MINKOWSKI_METRIC` | `def minkowskiMetric` |
| 8 | Kronecker delta δ^μ_ν | `kronecker_delta()` | — |
| 9 | Levi-Civita symbol ε_{μνρσ} | `levi_civita_4d()` | — |
| 10 | Christoffel symbols Γ^λ_{μ,ν} | `metric_christoffel()` | `ChristoffelSymbols` |
| 11 | Riemann tensor R^ρ_{σ,μ,ν} | `compute_riemann_tensor()` | `RiemannTensor` |
| 12 | Ricci tensor R_{μ,ν} | `compute_ricci_tensor()` | `RicciTensor` |
| 13 | Ricci scalar R | `compute_ricci_scalar()` | `RicciScalar` |
| 14 | Einstein tensor G_{μ,ν} | `compute_einstein_tensor()` | `EinsteinTensor` |
| 15 | Weyl tensor C_{ρ,σ,μ,ν} | `compute_weyl_tensor()` | — |
| 16 | Stress-energy tensor T_{μ,ν} | Various models | `StressEnergyTensor` |
| 17 | FLRW metric parameters | `CosmologicalParams` | — |
| 18 | Metric perturbation h_{μ,ν} | `MetricPerturbation` | — |
| 19 | Geodesic state | `GeodesicState` | `structure GeodesicEquation` |

## L2: Core Concepts — Complete ✅

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Principle of general covariance | Invoked via gauge_transform_perturbation() |
| 2 | Equivalence principle | Geodesic equation — free fall |
| 3 | Spacetime curvature from matter | Einstein equation residual |
| 4 | Geodesic principle | geodesic_equation_rhs() |
| 5 | Diffeomorphism invariance | Einstein-Hilbert action (scalar density) |
| 6 | Energy conditions (WEC, NEC, SEC, DEC) | energy_conditions_check() |
| 7 | Coordinate singularity vs physical | Kretschmann scalar, EF coordinates |
| 8 | Maximal analytic extension | Kruskal-Szekeres metric |
| 9 | Gravitational wave polarizations | TT gauge projection (h_+, h_×) |
| 10 | Black hole uniqueness ("no hair") | Kerr metric (M, a parameters) |
| 11 | Cosmic expansion | Friedmann equations |
| 12 | Inflation (exponential expansion) | de Sitter scale factor |

## L3: Mathematical Structures — Complete ✅

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Pseudo-Riemannian manifold | LorentzianMetric with signature check |
| 2 | Tensor algebra (rank 0-4) | All contraction/symmetrization ops |
| 3 | Index raising/lowering | raise/lower_first_index, *_vector |
| 4 | Levi-Civita connection | metric_christoffel(), metric_compatibility_check() |
| 5 | Covariant derivative | covariant_derivative_vector/oneform/tensor2 |
| 6 | Lie derivative | (implicit in gauge transforms) |
| 7 | Tensor symmetrization | symmetrize_tensor2(), antisymmetrize_tensor2() |
| 8 | Frobenius norm | norm_tensor2() |
| 9 | Matrix algebra (4×4) | determinant, cofactor, inverse, matmul |
| 10 | d'Alembertian operator | dalembertian_scalar() |
| 11 | Coordinate transformations | transform_tensor2_covariant(), compute_jacobian() |
| 12 | Christoffel transformation | transform_christoffel() |

## L4: Fundamental Laws — Complete ✅

| # | Law/Equation | C Implementation | Lean Statement |
|---|-------------|-----------------|----------------|
| 1 | Einstein field eqs: G_{μ,ν}+Λg_{μ,ν}=κT_{μ,ν} | einstein_equation_residual() | einsteinFieldEquationHolds |
| 2 | Vacuum Einstein: G_{μ,ν}=0 or R_{μ,ν}=0 | einstein_vacuum_check() | vacuumEinsteinHolds |
| 3 | Trace eq: -R+4Λ=κT | einstein_trace_equation() | — |
| 4 | Trace-reversed: R_{μ,ν}=κ(T_{μ,ν}-½g_{μ,ν}T) | einstein_trace_reversed() | TraceReversedEinsteinEquivalence |
| 5 | Einstein-Hilbert action: S ∝ ∫R√(-g)d⁴x | einstein_hilbert_lagrangian_density() | — |
| 6 | Algebraic Bianchi: R^ρ_{[σ,μ,ν]}=0 | check_algebraic_bianchi() | — |
| 7 | Contracted Bianchi: ∇^μ G_{μ,ν}=0 | check_contracted_bianchi() | contractedBianchiHolds |
| 8 | Geodesic equation: d²x^μ/dτ²+Γ^μ_{αβ}u^αu^β=0 | geodesic_equation_rhs() | GeodesicEquation |
| 9 | Geodesic deviation: ∇_U∇_U ξ^μ=R^μ_{νρσ}U^νU^ρξ^σ | geodesic_deviation() | — |
| 10 | Energy conservation: ∇_μ T^{μ,ν}=0 | stress_energy_divergence() | — |
| 11 | Birkhoff theorem | (documented) | BirkhoffTheorem |
| 12 | Linearized Einstein: □h̅_{μ,ν}=-16πG/c⁴ T_{μ,ν} | linearized_einstein_equation() | — |
| 13 | Friedmann eq 1: H²=H₀²(ΣΩ_i a^{-n_i}) | friedmann_hubble() | — |
| 14 | Friedmann eq 2: ä/a ∝ -(ρ+3p) | friedmann_acceleration() | — |

## L5: Computational Methods — Complete ✅

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Matrix inversion (Gauss-Jordan) | metric_compute_inverse() |
| 2 | Determinant (Laplace expansion) | matrix4x4_determinant() |
| 3 | Numerical differentiation (central diff) | compute_jacobian(), tests |
| 4 | RK4 ODE integration | geodesic_rk4_step(), scale_factor_rk4_integrate() |
| 5 | Newton iteration (root finding) | kruskal_szekeres_metric(), inverse_tortoise_coordinate() |
| 6 | Numerical quadrature (Riemann sum) | cosmic_time_at_a(), comoving_distance() |
| 7 | Characteristic equation (eigenvalue check) | metric_check_lorentz_signature() |
| 8 | Cofactor matrix computation | matrix4x4_cofactor() |

## L6: Canonical Systems — Complete ✅

| # | System | Implementation |
|---|--------|---------------|
| 1 | Schwarzschild metric | schwarzschild_metric() |
| 2 | Schwarzschild effective potential | schwarzschild_eff_potential() |
| 3 | ISCO at r=6M | schwarzschild_isco() |
| 4 | Photon sphere at r=3M | schwarzschild_photon_sphere() |
| 5 | Circular orbit analysis | schwarzschild_circular_orbits() |
| 6 | Tortoise coordinate | tortoise_coordinate() |
| 7 | Eddington-Finkelstein coordinates | eddington_finkelstein_ingoing/outgoing() |
| 8 | Kruskal-Szekeres extension | kruskal_szekeres_metric() |
| 9 | Kerr metric (Boyer-Lindquist) | kerr_metric() |
| 10 | Kerr horizons and ergosphere | kerr_horizon/cauchy/ergosphere_radius() |
| 11 | FLRW metric | flrw_metric() |
| 12 | Perfect fluid stress-energy | stress_energy_perfect_fluid() |
| 13 | Dust (pressureless fluid) | stress_energy_dust() |
| 14 | Radiation (p=ρ/3) | stress_energy_radiation() |
| 15 | EM stress-energy | stress_energy_electromagnetic() |
| 16 | Scalar field stress-energy | stress_energy_scalar_field() |
| 17 | Gravitational waves (TT gauge) | tt_gauge_project() |
| 18 | Quadrupole radiation formula | gravitational_wave_strain(), quadrupole_power() |
| 19 | Schwarzschild isotropic coordinates | schwarzschild_isotropic_metric() |

## L7: Applications — Complete ✅ (3+ applications)

| # | Application | Implementation |
|---|------------|---------------|
| 1 | Age of universe (ΛCDM) | age_of_universe() → 13.79 Gyr |
| 2 | Cosmological distances (d_L, d_A) | luminosity/angular_diameter_distance() |
| 3 | Deceleration parameter evolution | deceleration_parameter() |
| 4 | Critical density of universe | critical_density() |
| 5 | GW strain from compact binaries | gravitational_wave_strain() |
| 6 | GW150914 analysis | demo_gravitational_waves.c |
| 7 | Lookback time | lookback_time() |
| 8 | Full Schwarzschild orbit integration | demo_schwarzschild.c |
| 9 | FLRW scale factor evolution | demo_cosmology.c |

## L8: Advanced Topics — Complete ✅ (4+ topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Bekenstein-Hawking entropy | bekenstein_hawking_entropy_schwarzschild/kerr() |
| 2 | Hawking temperature | hawking_temperature_schwarzschild/kerr() |
| 3 | Black hole evaporation time | black_hole_evaporation_time_s() |
| 4 | Four laws of BH mechanics | zeroth/first/second/third_law_bh_*() |
| 5 | Kerr surface gravity | hawking_temperature_kerr() |
| 6 | Gravitational memory effect | gravitational_memory() |
| 7 | Chern-Pontryagin invariant | compute_chern_pontryagin() |
| 8 | de Sitter inflation | de_sitter_scale_factor(), inflation_e_folds() |
| 9 | Komar mass (ADM mass) | komar_mass() |
| 10 | Surface gravity | surface_gravity_schwarzschild() |

## L9: Research Frontiers — Partial ✅

| # | Topic | Status |
|---|-------|--------|
| 1 | Quantum gravity & Einstein equations | Documented in einstein.lean |
| 2 | Dark energy / cosmological constant problem | Documented in einstein.lean |
| 3 | Black hole information paradox | information_paradox_status() — documented |
| 4 | GW memory detection prospects | Documented in linearized.c |
| 5 | Holographic principle | Documented in bh_thermo.c |

## Layer Summary

| Level | Status | Evidence |
|-------|--------|---------|
| L1 Definitions | **Complete** | 19 typedef/struct/def in C + Lean |
| L2 Core Concepts | **Complete** | 12 concepts implemented |
| L3 Math Structures | **Complete** | 12 structures with full operations |
| L4 Fundamental Laws | **Complete** | 14 equations verified + Lean theorems |
| L5 Algorithms | **Complete** | 8 computational methods |
| L6 Canonical Systems | **Complete** | 19 canonical GR systems |
| L7 Applications | **Complete** | 9 real-world applications |
| L8 Advanced Topics | **Complete** | 10 advanced GR topics |
| L9 Research Frontiers | **Partial** | 5 frontiers documented |
