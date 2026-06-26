# mini-differential-geometry — 微分几何

> Manifolds · Tensors · Curvature · Geodesics · Differential Forms
> Reference: Wald, Carroll, do Carmo, Nakahara
> Courses: MIT 8.962 · Stanford PHYSICS 230 · Cambridge Part III · Princeton PHY 535
>          Caltech Ph 205 · Oxford CMT · ETH 402-0891

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (18 core definitions)
- **L2 Core Concepts**: Complete (8 concepts)
- **L3 Mathematical Structures**: Complete (10 structures)
- **L4 Fundamental Laws**: Complete (8 theorems/laws)
- **L5 Computational Methods**: Complete (6 methods)
- **L6 Canonical Systems**: Complete (5 classic systems)
- **L7 Applications**: Complete (4 applications: GW/LIGO, GPS, Gravity Probe B, FRW cosmology)
- **L8 Advanced Topics**: Partial+ (3 advanced topics)
- **L9 Research Frontiers**: Partial (documented)

**Score: 16/18 · include/ + src/ = 3470 lines (≥ 3000) · 45 tests passing**

## Core Definitions (L1)

| Concept | Type/Function |
|---------|--------------|
| Smooth Manifold | `Manifold` struct |
| Coordinate Chart | `Chart` struct, `chart_create()` |
| Atlas | `manifold_add_chart()`, transition maps |
| Tangent Vector | `TangentVector` — X = X^μ ∂_μ |
| Cotangent Vector | `CotangentVector` — ω = ω_μ dx^μ |
| Tensor (r,s) | `Tensor` struct with full algebra |
| Metric Tensor | `Metric` — g_{μν} with signature |
| Christoffel Symbols | `Connection` — Γ^ρ_{μν} |
| Riemann Tensor | `RiemannTensor` — R^ρ_{σμν} |
| Ricci Tensor | `RicciTensor` — R_{μν} |
| Einstein Tensor | G_{μν} = R_{μν} - ½Rg_{μν} |
| Geodesic | `GeodesicPoint`, RK4 integrator |
| Differential p-form | `DifferentialForm`, wedge, d, *, interior product |
| Lie Derivative | L_X on scalars, vectors, 1-forms, (0,2) tensors |
| Killing Vector | ξ: L_ξ g = 0, conserved quantities |

## Core Theorems (L4)

| Theorem | Formula | Verification |
|---------|---------|-------------|
| Levi-Civita Connection | Γ^ρ_{μν} = ½g^{ρσ}(∂_μg_{σν} + ∂_νg_{σμ} - ∂_σg_{μν}) | `connection_from_metric()` |
| Riemann Tensor | R^ρ_{σμν} = ∂_μΓ^ρ_{νσ} - ∂_νΓ^ρ_{μσ} + Γ^ρ_{μλ}Γ^λ_{νσ} - Γ^ρ_{νλ}Γ^λ_{μσ} | `riemann_compute()` |
| First Bianchi Identity | R^ρ_{[σμν]} = 0 | `bianchi_identity_first_verify()` |
| Second Bianchi Identity | ∇_{[λ}R^ρ_{|σ|μν]} = 0 | `bianchi_identity_second_verify()` |
| Geodesic Equation | d²x^μ/dτ² + Γ^μ_{αβ}(dx^α/dτ)(dx^β/dτ) = 0 | `geodesic_rhs()` |
| Killing Equation | ∇_μξ_ν + ∇_νξ_μ = 0 | `killing_equation_check()` |
| Cartan Structure Eq. | de^a + ω^a_b ∧ e^b = 0 | `dg_tetrad.c` |
| Einstein Field Eq. | G_{μν} = 8πG T_{μν} | `einstein_tensor_compute()` |

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| RK4 Geodesic Integration | `geodesic_step_rk4()` | O(n_steps × dim³) |
| Matrix Inversion (Gauss) | `matrix_invert_4x4()` in dg_metric.c | O(dim³) |
| Determinant (Laplace) | `matrix_det()` in dg_metric.c | O(dim!) |
| Christoffel Computation | `connection_from_metric()` | O(dim⁴) |
| Riemann Tensor Computation | `riemann_compute()` | O(dim⁵) |
| Jacobi Eigenvalue Algorithm | `tetrad_from_metric()` in dg_tetrad.c | O(dim³) per sweep |

## Canonical Systems (L6)

| System | Example/Test |
|--------|-------------|
| Minkowski Spacetime | Flat spacetime curvature verification |
| Schwarzschild Metric | `geodesic_schwarzschild_init()`, example |
| FRW Cosmology | `geodesic_frw_null_init()`, example |
| Spherically Symmetric Spacetime | Curvature demo example |
| Maximally Symmetric Space | `max_killing_vectors()`, Killing checks |

## Directory Structure

```
mini-differential-geometry/
├── Makefile                          # make test / make examples
├── README.md                         # This file
├── include/                          # 8 header files (552 lines)
│   ├── dg_manifold.h                 # Manifold, Chart, Tangent Space
│   ├── dg_tensor.h                   # Tensor (r,s) algebra
│   ├── dg_metric.h                   # Metric tensor, signature
│   ├── dg_connection.h               # Christoffel, covariant derivative
│   ├── dg_curvature.h                # Riemann, Ricci, Einstein, Weyl
│   ├── dg_geodesic.h                 # Geodesics, parallel transport
│   ├── dg_forms.h                    # Differential forms, Hodge star
│   └── dg_lie.h                      # Lie derivative, Killing vectors
├── src/                              # 10 source files (2918 lines)
│   ├── dg_manifold.c                 # Manifold implementation
│   ├── dg_tensor.c                   # Tensor algebra
│   ├── dg_metric.c                   # Metric operations
│   ├── dg_connection.c               # Connection computation
│   ├── dg_curvature.c                # Curvature tensors
│   ├── dg_geodesic.c                 # Geodesic integration
│   ├── dg_forms.c                    # Exterior calculus
│   ├── dg_lie.c                      # Lie derivatives
│   ├── dg_tetrad.c                   # Tetrad/vierbein formalism
│   └── dg_gravwave.c                 # Gravitational waves, GPS, Lense-Thirring
├── tests/
│   └── test_core.c                   # 45 assertion-based tests
├── examples/
│   ├── schwarzschild_geodesic.c      # Schwarzschild orbit
│   ├── frw_cosmology.c               # FRW expansion
│   └── curvature_demo.c              # Curvature computation
├── docs/
│   ├── knowledge-graph.md            # L1-L9 knowledge map
│   ├── coverage-report.md            # Coverage assessment
│   ├── gap-report.md                 # Gap analysis
│   ├── course-alignment.md           # 9-school curriculum mapping
│   ├── course-tree.md                # Prerequisite dependency tree
│   └── cheatsheet.md                 # Quick reference
├── notebooks/
├── demos/
└── benches/
```

## Quick Start

```bash
make          # compile library
make test     # run 45 tests (all passing)
make examples # build and run examples
make clean    # remove build artifacts
make count    # line count summary
```

## Reference Textbooks

| Topic | Textbook | Author |
|-------|----------|--------|
| General Relativity | *General Relativity* | Wald (1984) |
| GR Introduction | *Spacetime and Geometry* | Carroll (2004) |
| Differential Geometry | *Differential Geometry of Curves and Surfaces* | do Carmo (1976) |
| Geometry & Physics | *Geometry, Topology and Physics* | Nakahara (2003) |
| Numerical Methods | *Computational Physics* | Thijssen (2007) |
| Cosmology | *Modern Cosmology* | Dodelson & Schmidt (2020) |
