# Course Tree — mini-black-holes

## Prerequisites (What You Need Before This Module)

```
mini-general-relativity/
├── mini-differential-geometry    ← manifolds, tensors, connection, curvature
├── mini-schwarzschild            ← spherical BH solution, geodesics
├── mini-kerr-metric              ← rotating BH, Boyer-Lindquist coordinates
├── mini-einstein-equations       ← EFEs, stress-energy tensor
└── mini-gravitational-waves      ← linearized GR, quadrupole formula
```

## Key Dependencies Within mini-black-holes

```
black_hole_metrics.h/c           ← foundation: all metric definitions
    ↓
black_hole_thermodynamics.h/c    ← depends on metrics for area, surface gravity
    ↓
black_hole_dynamics.h/c          ← depends on metrics for geodesics, potentials
    ↓
black_hole_waves.h/c             ← depends on metrics for QNM, ringdown
    ↓
black_hole_advanced.h/c          ← depends on thermo for holography, info paradox
```

## Concept Dependency Graph

```
Differential Geometry (manifolds, curvature)
    │
    ├── Schwarzschild Solution ──→ ISCO, photon sphere, gravitational redshift
    │       │
    │       └──→ BH Thermodynamics ──→ T_H, S_BH, four laws, evaporation
    │               │
    │               ├──→ Information Paradox ──→ Page curve, ER=EPR, firewall
    │               │
    │               └──→ Holographic Principle ──→ AdS/CFT, Bekenstein bound
    │
    ├── Kerr Solution ──→ Ergosphere, frame-dragging, Penrose process
    │       │
    │       └──→ Superradiance ──→ Boson clouds, spin-down
    │       │
    │       └──→ ISCO(a) ──→ Bardeen spin evolution
    │
    ├── Reissner-Nordström ──→ Charged BH thermodynamics
    │
    ├── Kerr-Newman ──→ General electrovac BH (no-hair theorem)
    │
    └── BH Perturbation Theory ──→ Teukolsky equation ──→ QNMs
            │
            └──→ Ringdown waveform ──→ GW data analysis, BH spectroscopy
```

## Postrequisites (Modules That Depend on This One)

```
mini-black-holes
    │
    ├──→ mini-cosmological-model    ← primordial BHs, dark matter
    ├──→ mini-numerical-relativity  ← BH merger simulations
    ├──→ mini-quantum-field-theory  ← Hawking radiation from QFT in curved space
    └──→ L9 Research Frontiers      ← quantum gravity, string theory, holography
```

## Learning Path (Suggested Order)

1. **Start**: `mini-differential-geometry` — understand manifolds and curvature
2. **Core**: `mini-schwarzschild` — the simplest BH, all essential concepts
3. **Extend**: `mini-kerr-metric` — rotation adds rich phenomenology
4. **Deepen**: `mini-black-holes` (this module) — thermodynamics, information, GWs
5. **Frontier**: Quantum gravity, AdS/CFT (L8-L9 topics)

## L9 Research Frontiers Dependency

```
String Theory / Quantum Gravity
    │
    ├── BH microstate counting (Strominger-Vafa 1996)
    ├── AdS/CFT correspondence (Maldacena 1997)
    ├── ER=EPR conjecture (Maldacena-Susskind 2013)
    ├── Firewall paradox (AMPS 2012)
    └── Weak gravity conjecture (Arkani-Hamed et al. 2007)
```
