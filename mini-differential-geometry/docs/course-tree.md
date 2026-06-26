# Course Tree — mini-differential-geometry

## Prerequisites

```
Linear Algebra ──┐
                 ├──→ Differential Geometry ──→ General Relativity
Real Analysis ───┤                              │
                 │                              ├──→ Cosmology
Topology ────────┘                              ├──→ Gravitational Waves
                                                ├──→ Black Hole Physics
                                                ├──→ Numerical Relativity
                                                └──→ Quantum Gravity (L9)
```

## Internal Dependency Tree

```
dg_manifold.h/c (Manifolds, Charts, Tangent Vectors)
    │
    └── dg_tensor.h/c (Tensor Algebra)
            │
            └── dg_metric.h/c (Metric Tensor)
                    │
                    ├── dg_connection.h/c (Christoffel Symbols, Covariant Derivative)
                    │       │
                    │       └── dg_curvature.h/c (Riemann, Ricci, Einstein, Weyl)
                    │               │
                    │               └── dg_geodesic.h/c (Geodesics, Parallel Transport)
                    │
                    ├── dg_forms.h/c (Differential Forms, Hodge Star)
                    │
                    └── dg_lie.h/c (Lie Derivatives, Killing Vectors)

dg_tetrad.c (depends on dg_metric.h, dg_connection.h, dg_curvature.h)
```

## Knowledge Map

| Module | L1 | L2 | L3 | L4 | L5 | L6 |
|--------|----|----|----|----|----|----|
| dg_manifold | Complete | Partial | — | — | Partial | — |
| dg_tensor | Complete | Complete | Complete | — | — | — |
| dg_metric | Complete | Complete | Complete | Complete | Complete | — |
| dg_connection | Complete | Complete | Complete | Complete | — | — |
| dg_curvature | Complete | Complete | Complete | Complete | — | — |
| dg_geodesic | Complete | Complete | Complete | Complete | Complete | Complete |
| dg_forms | Complete | Complete | Complete | Complete | Partial | — |
| dg_lie | Complete | Complete | Complete | Complete | — | — |
| dg_tetrad | — | Complete | Complete | Complete | Complete | — |

## Research Frontiers (L9)

Connections to active research areas:
- **Numerical Relativity**: Tetrad formalism, geodesic integrator as building blocks
- **Loop Quantum Gravity**: Tetrad variables (Ashtekar connection)
- **String Theory**: Riemannian geometry of Calabi-Yau manifolds
- **AdS/CFT**: Holographic dictionary requires differential geometry
- **Gravitational Wave Astronomy**: Geodesic deviation equation for detector response