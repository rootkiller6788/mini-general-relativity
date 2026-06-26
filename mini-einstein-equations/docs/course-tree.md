# Course Tree — mini-einstein-equations

## Prerequisite Dependency Tree

```
Mathematics Prerequisites
├── Calculus (multivariable)
│   ├── Partial derivatives
│   ├── Chain rule
│   └── Taylor expansion
├── Linear Algebra
│   ├── Matrix operations
│   ├── Determinants
│   ├── Eigenvalues/eigenvectors
│   └── Matrix inverses
├── Differential Equations
│   ├── ODE systems
│   ├── Numerical methods (RK4)
│   └── Boundary value problems
└── Special Relativity (SR)
    ├── Minkowski spacetime
    ├── Lorentz transformations
    ├── 4-vectors
    └── Relativistic mechanics
        └── Stress-energy tensor in SR

Physics Prerequisites
├── Classical Mechanics
│   ├── Lagrangian mechanics
│   ├── Hamiltonian mechanics
│   └── Action principle
├── Newtonian Gravity
│   ├── Poisson equation ∇²φ = 4πGρ
│   ├── Kepler orbits
│   └── Tidal forces
└── Electromagnetism
    ├── Maxwell equations
    ├── Gauge invariance
    └── Wave solutions

↓ (this module builds on the above)

mini-einstein-equations
├── L1: Tensor definitions (tensor.h)
│   └── → L2: Core concepts
├── L2: Core concepts
│   └── → L3: Mathematical structures
├── L3: Mathematical structures (metric.h)
│   ├── → L1: Definitions formalized
│   └── → L4: Fundamental laws
├── L4: Fundamental laws (einstein.h, curvature.h, geodesic.h)
│   ├── → L5: Computational methods
│   └── → L6: Canonical systems
├── L5: Computational methods
│   └── → L6/L7: Applications
├── L6: Canonical systems (coordinate.h, stress_energy.h)
│   └── → L7: Applications
├── L7: Applications (cosmology.h, linearized.h)
│   └── → L8: Advanced topics
├── L8: Advanced topics (bh_thermo.c)
│   └── → L9: Research frontiers
└── L9: Research frontiers (documented)
    └── → Future: mini-quantum-field-theory, mini-cosmology

↓ (modules that depend on this one)

Downstream Modules
├── 14. mini-quantum-field-theory (QFT in curved spacetime)
├── 12. mini-cosmology (FLRW, inflation, CMB)
├── 11. mini-astrophysics (BH astrophysics, GW astronomy)
└── 18. mini-interdisciplinary-physics (GR + thermodynamics)

## Internal File Dependencies

```
include/tensor.h          (no internal deps — base layer)
  ↑
include/metric.h          (depends: tensor.h)
  ↑
include/curvature.h       (depends: metric.h)
include/geodesic.h        (depends: curvature.h)
include/stress_energy.h   (depends: metric.h)
  ↑
include/einstein.h        (depends: curvature.h, stress_energy.h)
  ↑
include/cosmology.h       (depends: einstein.h)
include/linearized.h      (depends: einstein.h)
include/coordinate.h      (depends: metric.h)

src/tensor.c              ← include/tensor.h
src/metric.c              ← include/metric.h (uses tensor.h transitively)
src/curvature.c           ← include/curvature.h
src/geodesic.c            ← include/geodesic.h
src/stress_energy.c       ← include/stress_energy.h
src/einstein.c            ← include/einstein.h
src/cosmology.c           ← include/cosmology.h
src/linearized.c          ← include/linearized.h
src/coordinate.c          ← include/coordinate.h
src/bh_thermo.c           ← include/coordinate.h, include/einstein.h
src/einstein.lean         (standalone — formalizes C implementations)

tests/test_all.c          ← all headers
examples/demo_*.c         ← specific subsets
```

## Learning Path (Recommended Order)

1. **Start**: `tensor.h/c` → data types and basic operations
2. **Then**: `metric.h/c` → metric, Christoffel, covariant derivative
3. **Then**: `curvature.h/c` → Riemann → Ricci → Einstein → Weyl
4. **Core**: `einstein.h/c` → Einstein equations in all forms
5. **Motion**: `geodesic.h/c` → particle motion in curved spacetime
6. **Source**: `stress_energy.h/c` → matter models
7. **Solutions**: `coordinate.h/c` → Schwarzschild, Kerr, FLRW
8. **Observations**: `linearized.h/c`, `cosmology.h/c` → GWs, cosmology
9. **Advanced**: `bh_thermo.c` → quantum aspects of BHs
10. **Formal**: `einstein.lean` → Lean 4 formal verification
