# Course Tree — mini-gravitational-waves

## Prerequisite Dependencies

```
General Relativity (Wald Ch.2-3)
├── Differential Geometry (manifolds, metrics, curvature)
│   ├── Riemann tensor, Ricci tensor, Einstein tensor
│   └── Geodesic equation
├── Linearized Gravity ← THIS MODULE
│   ├── Minkowski perturbation g_{μν} = η_{μν} + h_{μν}
│   ├── Harmonic gauge ∂^μ h̄_{μν} = 0
│   ├── TT gauge (traceless, transverse)
│   ├── Quadrupole formula (Einstein 1918)
│   ├── GW energy-momentum (Isaacson 1968)
│   ├── Binary inspiral + PN expansion
│   ├── GW memory effect
│   ├── Stochastic GW background
│   └── Continuous GW sources
├── GW Detection (matched filtering)
│   ├── Noise-weighted inner product
│   ├── Fisher matrix parameter estimation
│   └── LIGO/LISA/ET noise models
└── GW Data Analysis Applications
    ├── GW150914 waveform modelling
    ├── Pulsar timing arrays (NANOGrav)
    └── Cosmological backgrounds
```

## Knowledge Dependencies

| This Module Requires | From Module |
|---------------------|------------|
| Minkowski metric, SR | mini-classical-mechanics (relativity) |
| Differential geometry | mini-differential-geometry |
| Einstein field equations | mini-einstein-equations |
| Schwarzschild/Kerr metrics | mini-schwarzschild, mini-kerr-metric |
| Numerical integration | mini-computational-physics |

## Downstream Dependencies

| Module | Uses This Module For |
|--------|---------------------|
| mini-black-holes | GW ringdown (QNM analysis) |
| mini-cosmological-model | Inflationary GWs, Ω_GW constraints |
| mini-numerical-relativity | Waveform extraction, NR-GW comparison |
| mini-astrophysics (parent) | Binary population synthesis, merger rates |
