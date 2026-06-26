# Course Tree — mini-numerical-relativity

## Prerequisites

```
mini-classical-mechanics
├── Newtonian gravity
├── Hamiltonian mechanics
└── Canonical formulation

mini-general-relativity
├── Einstein field equations
├── Schwarzschild solution
├── Kerr solution
└── Gravitational waves

mini-computational-physics
├── Finite difference methods
├── Runge-Kutta integration
└── Elliptic PDE solvers
```

## Internal Dependencies

```
nr_grid.h ──────────────────────────────────────────┐
  │ (grid, complex, tensor types)                   │
  ├── nr_utils.h ───────────────────────────────────┤
  │     │ (FD stencils, RK4, SOR, interpolation)     │
  │     ├── nr_adm.h ───────────────────────────────┤
  │     │     │ (3+1 decomposition, constraints)     │
  │     │     ├── nr_bssn.h ────────────────────────┤
  │     │     │     │ (conformal evolution system)   │
  │     │     │     └── nr_gauge.c (slicing, shift)  │
  │     │     └── nr_initial.h (initial data)        │
  │     ├── nr_horizon.h (apparent horizon finder)   │
  │     └── nr_wave.h (GW extraction, strain)        │
  └── [all modules use grid types]                   │
```

## Learning Path

1. **L1-L3**: Grid infrastructure → Mathematical structures (tensors, geometry)
2. **L4**: Einstein constraints → ADM evolution equations
3. **L5**: Finite differences → RK4 → Interpolation → SOR solver
4. **L2**: ADM → BSSN conversion → Conformal decomposition
5. **L6**: Schwarzschild → Brill-Lindquist → Bowen-York → Kerr-Schild
6. **L5**: Horizon finding → Wave extraction → Mode decomposition
7. **L7**: Full GW pipeline → LIGO analysis → Radiated energy
