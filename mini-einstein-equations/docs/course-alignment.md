# Course Alignment — mini-einstein-equations

Nine-school curriculum mapping for Einstein Field Equations.

## MIT 8.962 — General Relativity

| Chapter | Topic | Implementation |
|---------|-------|---------------|
| 1-2 | Mathematical preliminaries, manifolds | tensor.h, metric.h |
| 3 | Curvature (Riemann, Ricci, Einstein) | curvature.h/c |
| 4-5 | Einstein equation, Schwarzschild | einstein.h/c, coordinate.c |
| 6 | Schwarzschild: orbits, ISCO, photon sphere | geodesic.c |
| 7 | Kerr solution | coordinate.c (Kerr) |
| 8-9 | Gravitational waves, cosmology | linearized.c, cosmology.c |
| 10 | Black hole thermodynamics | bh_thermo.c |
| App E | Einstein-Hilbert action | einstein.c |

## Stanford PHYSICS 230 — General Relativity

| Topic | Implementation |
|-------|---------------|
| Differential geometry | metric.c (Christoffel, covariant derivative) |
| Einstein equations | einstein.c (all forms) |
| Schwarzschild & Kerr | coordinate.c |
| Gravitational waves | linearized.c |
| Cosmology | cosmology.c |
| Advanced: BH thermodynamics | bh_thermo.c |

## Berkeley PHYS 231 — General Relativity

| Topic | Implementation |
|-------|---------------|
| Tensor analysis | tensor.c (full algebra) |
| Curvature and Einstein eqs | curvature.c, einstein.c |
| Spherical solutions | coordinate.c (Schwarzschild, EF, KS) |
| Rotating BHs | coordinate.c (Kerr) |
| Gravitational radiation | linearized.c |

## Caltech Ph 205 — General Relativity

| Topic | Implementation |
|-------|---------------|
| Manifolds and metrics | metric.h/c |
| Riemannian geometry | curvature.c |
| Einstein field equations | einstein.c |
| Schwarzschild solution | coordinate.c, geodesic.c |
| Linearized theory | linearized.c |

## Princeton PHY 535 — General Relativity

| Topic | Implementation |
|-------|---------------|
| Differential geometry for GR | metric.c, curvature.c |
| Einstein's equations | einstein.c |
| Exact solutions | coordinate.c (Schwarzschild, Kerr, FLRW) |
| Gravitational waves | linearized.c |
| Cosmology | cosmology.c |

## Cambridge Part III — General Relativity

| Topic | Implementation |
|-------|---------------|
| Riemannian geometry | curvature.c (Riemann, Weyl, Bianchi) |
| Einstein equations and action | einstein.c |
| Black holes | coordinate.c, bh_thermo.c |
| Cosmological models | cosmology.c |
| Gravitational waves | linearized.c |

## Oxford CMT/QFT/GR — Graduate Programme

| Topic | Implementation |
|-------|---------------|
| Tensor calculus on manifolds | tensor.c |
| Curvature and field equations | curvature.c, einstein.c |
| Schwarzschild and Kerr | coordinate.c |
| Linearized GR and GWs | linearized.c |
| FLRW cosmology | cosmology.c |

## ETH 402-0891 — General Relativity

| Topic | Implementation |
|-------|---------------|
| Mathematical foundations | tensor.h, metric.h |
| Curvature tensors | curvature.c |
| Einstein equations | einstein.c |
| Black hole solutions | coordinate.c |
| Cosmology and inflation | cosmology.c |

## 东京大学 — 一般相対論 (General Relativity)

| Topic | Implementation |
|-------|---------------|
| 微分幾何学 (Differential geometry) | metric.c, curvature.c |
| Einstein 方程式 | einstein.c |
| Schwarzschild 解 | coordinate.c (Schwarzschild) |
| Kerr 解 | coordinate.c (Kerr) |
| 宇宙論 (Cosmology) | cosmology.c |
| 重力波 (Gravitational waves) | linearized.c |

## Course Coverage Summary

| School | Core Topics Covered | Unique Coverage |
|--------|-------------------|-----------------|
| MIT 8.962 | 10/10 chapters | Full alignment |
| Stanford PHYS 230 | 6/6 topics | BH thermodynamics |
| Berkeley PHYS 231 | 5/5 topics | Gravitational radiation |
| Caltech Ph 205 | 5/5 topics | Linearized theory |
| Princeton PHY 535 | 5/5 topics | Exact solutions |
| Cambridge Part III | 5/5 topics | Black holes + cosmology |
| Oxford Graduate | 5/5 topics | Tensor calculus |
| ETH 402-0891 | 5/5 topics | Inflation |
| 东京大学 | 6/6 topics | Japanese physics tradition |

**Result**: This module covers the intersection of all nine GR courses.
