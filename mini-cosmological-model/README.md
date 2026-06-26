# mini-cosmological-model — 宇宙学模型

> Reference: MIT 8.962 / Wald Ch.5 / Carroll Ch.8 / Dodelson & Schmidt (2020)
> Planck 2018 Cosmological Parameters (A&A 641, A6)

## Module Status: COMPLETE ✅

| Level | Status | Score |
|-------|--------|-------|
| L1 Definitions | Complete | 2 |
| L2 Core Concepts | Complete | 2 |
| L3 Math Structures | Complete | 2 |
| L4 Fundamental Laws | Complete | 2 |
| L5 Computational Methods | Complete | 2 |
| L6 Canonical Systems | Complete | 2 |
| L7 Applications | Complete | 2 |
| L8 Advanced Topics | Partial | 1 |
| L9 Research Frontiers | Partial | 1 |

**Total: 16/18 → COMPLETE**

## Core Definitions (L1)

| Definition | Symbol | Implementation |
|-----------|--------|----------------|
| Hubble parameter | H(z), H₀, h | `friedmann_H_z()`, `CosmoModel.H0` |
| Density parameters | Ωₘ, Ωᵣ, ΩΛ, Ωₖ | `CosmoModel.Omega_*` |
| Redshift | z | `flrw_redshift_from_scale()` |
| Scale factor | a(t) | `CosmoState.a`, `flrw_scale_from_redshift()` |
| Equation of state | w = P/(ρc²) | `friedmann_w_eff()` |
| Deceleration parameter | q | `flrw_deceleration_parameter()` |
| Jerk parameter | j | `flrw_jerk_parameter()` |
| Primordial spectrum | P_R(k) | `ps_primordial_dimensionless()` |
| Critical density | ρ_crit = 3H²/(8πG) | `cosmo_rho_critical()` |
| Proper/comoving distance | d_proper, r_comoving | `flrw_proper_distance()` |

## Core Theorems (L4)

| Theorem | Formula | Implementation | Lean Proof |
|---------|---------|----------------|-------------|
| Friedmann 1st | H² = (8πG/3)ρ − kc²/a² + Λc²/3 | `friedmann_E_z()` | `E_squared` |
| Friedmann 2nd | ä/a = −(4πG/3)(ρ+3P/c²) + Λc²/3 | `friedmann_dH_da()` | — |
| Continuity eq | ρ̇ + 3H(ρ+P/c²) = 0 | `friedmann_rho_a()` | `densityAtScale` |
| Etherington | D_L = (1+z)² D_A | `flrw_luminosity_distance()` | `etherington_reciprocity` |
| Sum rule | Σ Ωᵢ = 1 | `cosmo_model_validate()` | `sum_rule_flat` |
| Poisson eq | k²Φ = −(3/2)ΩₘH₀²δ/a | `transfer_potential_to_density()` | — |

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| RK4 integration | `friedmann_rk4_step()` | O(1) per step |
| Adaptive RK4 | `friedmann_integrate_scale_factor()` | O(n_steps) |
| Simpson integration | `friedmann_cosmic_time()` | O(n) |
| Gauss-Legendre quadrature (5-pt) | `distances_adaptive_integrate()` | O(2^k·5) |
| Spherical Bessel recurrence | `perturbations_cmb_Cl_TT_SW()` | O(l·n_k) |
| Halofit nonlinear correction | `ps_nonlinear_halofit()` | O(n_iter) |

## Canonical Systems (L6)

- Einstein-de Sitter (Ωₘ=1)
- Radiation-dominated (Ωᵣ=1)
- de Sitter (ΩΛ=1)
- Planck 2018 ΛCDM (Ωₘ=0.315, ΩΛ=0.685)
- Matter-radiation equality (z_eq ≈ 3400)
- Acceleration onset (z_acc ≈ 0.67)
- Particle horizon (~46 Glyr)
- Event horizon (~16.5 Glyr)
- Age of universe (13.79 Gyr)
- Distance ladder: D_C, D_M, D_A, D_L, μ

## Nine-School Course Mapping

| School | Course | Topics |
|--------|--------|--------|
| MIT | 8.962 GR / 8.942 Cosmology | FLRW, Friedmann, CMB, LSS |
| Stanford | PHYSICS 230 GR | Conformal time, perturbations, inflation |
| Berkeley | PHYS 231 GR / 242 CM | Friedmann solutions, distance ladder |
| Caltech | Ph 205 GR | de Sitter, quantum fluctuations, CMB pol |
| Princeton | PHY 535 GR | FLRW dynamics, Etherington, Sachs-Wolfe |
| Cambridge | Part III Cosmology | Redshift-distance, growth, BAO |
| Oxford | CMT/QFT/GR | Mathematical cosmology, cosmography |
| ETH | 402-0891 GR | Friedmann-Lemaitre, transfer functions |
| 东京大学 | 宇宙論 | Hubble tension, S8, neutrinos |

## Build & Test

```bash
make          # Build all and run tests
make test     # Run test suites (17 tests)
make examples # Build examples
make clean    # Remove build artifacts

# Run examples
./build/examples/age_of_universe
./build/examples/distance_redshift
./build/examples/cosmic_evolution
```

## Directory Structure

```
mini-cosmological-model/
├── README.md              ← This file (COMPLETE status)
├── Makefile               ← make test passes (17/17)
├── include/               ← 6 headers, 3141 total lines with src/
│   ├── cosmo_model.h          Core types, constants, parameters
│   ├── flrw.h                 FLRW metric, geometry
│   ├── friedmann.h            Friedmann equations, ODE system
│   ├── distances.h            Cosmological distance measures
│   ├── perturbations.h        Density perturbation theory
│   └── power_spectrum.h       Primordial and matter power spectrum
├── src/                   ← 8 C files + 1 Lean file
│   ├── cosmo_model.c          Model init, validation, age
│   ├── flrw.c                 Metric components, geometry
│   ├── friedmann.c            E(z), H(z), RK4, analytic solutions
│   ├── distances.c            All distance measures, horizons
│   ├── perturbations.c        Growth factor, PS mass fn, BAO, CMB
│   ├── power_spectrum.c       P(k), transfer fn, halofit
│   ├── transfer.c             CMB anisotropy, SW, ISW, lensing
│   └── cosmology.lean         Lean 4 formalization (14 theorems)
├── tests/                 ← 2 test files, 17 assertions
│   ├── test_friedmann.c       Friedmann equation tests (10)
│   └── test_distances.c       Distance measure tests (7)
├── examples/              ← 3 end-to-end examples (>30 lines each)
│   ├── age_of_universe.c      Age of universe calculator
│   ├── distance_redshift.c    Distance-redshift relations
│   └── cosmic_evolution.c     Expansion history tracker
└── docs/                  ← 5 knowledge documents
    ├── knowledge-graph.md     L1-L9 full coverage map
    ├── coverage-report.md     Per-layer assessment (16/18)
    ├── gap-report.md          Gap analysis (no critical gaps)
    ├── course-alignment.md    Nine-school curriculum mapping
    └── course-tree.md         Prerequisite dependency tree
```

## Reference Textbooks

| Textbook | Authors | Year |
|----------|---------|------|
| Modern Cosmology | Dodelson & Schmidt | 2020 |
| Spacetime and Geometry | Carroll | 2004 |
| General Relativity | Wald | 1984 |
| Cosmology | Weinberg | 2008 |
| Principles of Physical Cosmology | Peebles | 1993 |
| Galaxy Formation and Evolution | Mo, van den Bosch, White | 2010 |
