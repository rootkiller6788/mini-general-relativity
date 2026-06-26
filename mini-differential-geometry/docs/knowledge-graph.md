# Knowledge Graph — mini-differential-geometry

## L1: Definitions (Complete)

| # | Concept | C Implementation | Lean/Formal |
|---|---------|-----------------|-------------|
| 1 | Smooth Manifold | `Manifold` struct, `manifold_create()` | — |
| 2 | Coordinate Chart | `Chart` struct, `chart_create()` | — |
| 3 | Atlas | `manifold_add_chart()`, transition maps | — |
| 4 | Tangent Vector | `TangentVector` struct, `tangent_vector_create()` | — |
| 5 | Cotangent Vector (1-form) | `CotangentVector` struct, `cotangent_vector_create()` | — |
| 6 | Tensor (r,s) | `Tensor` struct, `tensor_alloc/create()` | — |
| 7 | Metric Tensor | `Metric` struct, `metric_create()` | — |
| 8 | Metric Signature | `MetricSignature` enum | — |
| 9 | Connection Coefficients | `Connection` struct, `connection_alloc()` | — |
| 10 | Christoffel Symbols | `connection_from_metric()` | — |
| 11 | Riemann Curvature Tensor | `RiemannTensor` struct, `riemann_compute()` | — |
| 12 | Ricci Tensor | `RicciTensor` struct, `ricci_compute()` | — |
| 13 | Ricci Scalar | `ricci_scalar_compute()` | — |
| 14 | Einstein Tensor | `einstein_tensor_compute()` | — |
| 15 | Geodesic | `GeodesicPoint` struct, `geodesic_integrate()` | — |
| 16 | Differential p-form | `DifferentialForm` struct, `form_alloc()` | — |
| 17 | Lie Derivative | `lie_derivative_scalar/vector/covector/tensor02()` | — |
| 18 | Killing Vector Field | `killing_equation_check()` | — |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Coordinate Independence | `tensor_transform()` — tensor transformation law |
| 2 | Metric Compatibility | `connection_verify_metric_compatibility()` |
| 3 | Torsion-Free Condition | `connection_compute_torsion()` |
| 4 | Geodesic Deviation | `geodesic_deviation_rhs()` |
| 5 | Parallel Transport | `parallel_transport()` |
| 6 | Proper Time | `proper_time_total()`, `metric_proper_time_inc()` |
| 7 | Index Raising/Lowering | `tensor_raise_index()`, `tensor_lower_index()` |
| 8 | Tetrad/Vierbein | `tetrad_alloc_from_metric()` (src/dg_tetrad.c) |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Tensor Algebra | `tensor_add/subtract/scale/product()` |
| 2 | Tensor Contraction | `tensor_contract()` |
| 3 | Symmetrization | `tensor_symmetrize_02()`, `antisymmetrize_02()` |
| 4 | Pseudo-Riemannian Structure | `MetricSignature`, `metric_verify_signature()` |
| 5 | Weyl Tensor | `weyl_tensor_compute()` |
| 6 | Sectional Curvature | `sectional_curvature()` |
| 7 | Wedge Product | `form_wedge()` |
| 8 | Exterior Derivative | `form_exterior_derivative()` |
| 9 | Hodge Star | `form_hodge_star()` |
| 10 | Spin Connection | `spin_connection_compute()` (src/dg_tetrad.c) |

## L4: Fundamental Laws (Complete)

| # | Law/Theorem | Implementation |
|---|-------------|---------------|
| 1 | Levi-Civita Connection Formula | `connection_from_metric()` |
| 2 | Riemann Tensor Definition | `riemann_compute()` |
| 3 | First Bianchi Identity | `bianchi_identity_first_verify()` |
| 4 | Second Bianchi Identity | `bianchi_identity_second_verify()` |
| 5 | Geodesic Equation | `geodesic_rhs()`, `geodesic_step_rk4()` |
| 6 | Metric Determinant Formula | `metric_determinant()`, `metric_recompute_inv_det()` |
| 7 | Killing Equation | `killing_equation_check()` |
| 8 | Cartan Structure Equations | `curvature_2form_compute()` (src/dg_tetrad.c) |

## L5: Computational Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | RK4 Geodesic Integration | `geodesic_step_rk4()`, `geodesic_integrate()` |
| 2 | Matrix Inversion (Gaussian Elimination) | `matrix_invert_4x4()` in dg_metric.c |
| 3 | Determinant (Laplace Expansion) | `matrix_det()` in dg_metric.c |
| 4 | Jacobi Eigenvalue Algorithm | `tetrad_from_metric()` in dg_tetrad.c |
| 5 | Finite Difference Derivatives | `connection_from_metric_fd()`, `form_exterior_derivative_fd()` |
| 6 | Transition Map Jacobian (FD) | `manifold_transition_jacobian()` |

## L6: Canonical Systems (Complete)

| # | System | Implementation |
|---|--------|---------------|
| 1 | Minkowski Spacetime | All zero-curvature verification tests |
| 2 | Schwarzschild Geodesics | `geodesic_schwarzschild_init()`, example |
| 3 | FRW Cosmological Metric | `geodesic_frw_null_init()`, example |
| 4 | Spherically Symmetric Curvature | Curvature demo example |
| 5 | Killing Vectors (Schwarzschild) | `killing_schwarzschild()` |

## L7: Applications (Partial+)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | Gravitational Lensing (light bending) | Schwarzschild null geodesic framework |
| 2 | Cosmic Expansion (FRW) | FRW cosmology example |
| 3 | Curvature Singularity Detection | `kretschmann_scalar()` |

## L8: Advanced Topics (Partial+)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Weyl Conformal Tensor | `weyl_tensor_compute()` |
| 2 | Tetrad/Vierbein Formalism | `src/dg_tetrad.c` |
| 3 | de Rham Cohomology (basic) | `form_is_closed()`, `form_is_exact()` |

## L9: Research Frontiers (Partial — documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Numerical Relativity (ADM) | Documented in tetrad formalism |
| 2 | Quantum Gravity Connections | Referenced in docs |
| 3 | Holographic Principle | Referenced in course-tree.md |