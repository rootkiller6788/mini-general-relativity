# mini-numerical-relativity — Numerical Relativity

> Solving Einstein's equations on a computer: 3+1 decomposition, BSSN evolution,
> black hole initial data, horizon finding, and gravitational wave extraction.

**Reference**: MIT 8.962 / Baumgarte & Shapiro (2010) / Alcubierre (2008)

---

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (14 core definitions)
- **L2 Core Concepts**: Complete (9 concepts)
- **L3 Math Structures**: Complete (8 structures)
- **L4 Fundamental Laws**: Complete (8 equations implemented)
- **L5 Computational Methods**: Complete (14 algorithms)
- **L6 Canonical Systems**: Complete (10 systems)
- **L7 Applications**: Complete (6 applications, 3 runnable examples)
- **L8 Advanced Topics**: Partial (4/6 implemented)
- **L9 Research Frontiers**: Partial (documented)
- **Submodule Lines**: 7010 (include/ + src/), ≥ 3000 ✅

---

## Core Definitions (L1)

| Symbol | Name | C Type | Header |
|--------|------|--------|--------|
| α | Lapse function | `nr_gf_t` | `nr_adm.h` |
| β^i | Shift vector | `nr_vector3_t` | `nr_grid.h` |
| γ_{ij} | 3-metric | `nr_sym_tensor3_t` | `nr_adm.h` |
| K_{ij} | Extrinsic curvature | `nr_sym_tensor3_t` | `nr_adm.h` |
| φ | Conformal factor | `nr_gf_t` | `nr_bssn.h` |
| γ̃_{ij} | Conformal metric | `nr_sym_tensor3_t` | `nr_bssn.h` |
| Ã_{ij} | Trace-free curvature | `nr_sym_tensor3_t` | `nr_bssn.h` |
| Γ̃^i | Conformal connection | `nr_vector3_t` | `nr_bssn.h` |
| Ψ_4 | Weyl scalar | `nr_complex_t` | `nr_wave.h` |
| h_+, h_× | GW strain | `nr_strain_series_t` | `nr_wave.h` |

## Core Theorems/Equations (L4)

| Equation | Formula | Implementation |
|----------|---------|---------------|
| ADM 3+1 metric | ds² = −α²dt² + γ_{ij}(dx^i+β^idt)(dx^j+β^jdt) | `nr_adm.h` |
| ADM evolution ∂_t γ_{ij} | −2αK_{ij} + D_iβ_j + D_jβ_i | `nr_adm_rhs_point()` |
| Hamiltonian constraint | H = R + K² − K_{ij}K^{ij} = 0 | `nr_adm_hamiltonian_constraint()` |
| Momentum constraint | M_i = D_jK^j_i − D_iK = 0 | `nr_adm_momentum_constraint()` |
| BSSN ∂_t φ | −αK/6 + β^k∂_kφ + ∂_kβ^k/6 | `nr_bssn_rhs_point()` |
| 1+log slicing | ∂_tα = −2αK + β^i∂_iα | `nr_gauge_1pluslog_point()` |
| Gamma-driver | ∂_tβ^i = 3B^i/4 + β^j∂_jβ^i | `nr_gauge_gamma_driver_point()` |
| Newman-Penrose Ψ_4 | C_{αβγδ}n^αm̄^βn^γm̄^δ | `nr_wave_psi4_point()` |

## Core Algorithms (L5)

| Algorithm | Complexity | Implementation |
|-----------|-----------|---------------|
| FD 2nd/4th/6th/8th order | O(1)/pt | `nr_fd_deriv1()`, `nr_fd_deriv2()` |
| RK4 time integration | O(N³) | `nr_rk4_step()` |
| Trilinear interpolation | O(1) | `nr_interp_trilinear()` |
| SOR elliptic solver | O(N³·iter) | `nr_sor_solve_3d()` |
| Kreiss-Oliger dissipation | O(N³) | `nr_kreiss_oliger_gf()` |
| Flow method AH finder | O(N_iter·N_ang) | `nr_horizon_find_flow()` |
| SWSH mode decomposition | O(N_θ·N_φ·L²) | `nr_wave_decompose_psi4()` |
| Fixed-frequency integration | O(N_t) | `nr_wave_integrate_psi4()` |

## Classic Problems (L6)

1. Schwarzschild in isotropic coords → `example_schwarzschild`
2. Binary BH puncture data → `example_puncture`
3. GW extraction & analysis → `example_wave_extract`

## Nine-School Curriculum Mapping

| School | Course | Coverage |
|--------|--------|----------|
| MIT | 8.962 GR | ADM, Schwarzschild, Kerr, GW |
| Stanford | PHYSICS 230 | Initial value problem |
| Berkeley | PHYS 231 | Numerical methods |
| Caltech | Ph 205 | BH physics, GW extraction |
| Princeton | PHY 535 | 3+1 formalism |
| Cambridge | Part III | BSSN, puncture method |
| Oxford | CMT NR | Horizon finding, wave extraction |
| ETH | 402-0891 | Computational methods |
| Tokyo | Astrophysics | Binary BH simulations |

## Directory Structure

```
mini-numerical-relativity/
├── Makefile              # make test → 11/11 passed
├── README.md             # This file
├── include/              # 6 header files
│   ├── nr_grid.h         # Grid + complex type
│   ├── nr_adm.h          # ADM 3+1 decomposition
│   ├── nr_bssn.h         # BSSN formulation
│   ├── nr_initial.h      # Initial data
│   ├── nr_horizon.h      # Apparent horizon
│   ├── nr_utils.h        # FD, RK4, interpolation
│   └── nr_wave.h         # GW extraction
├── src/                  # 8 source files
│   ├── nr_grid.c         # Grid memory management
│   ├── nr_adm.c          # ADM constraints & evolution
│   ├── nr_bssn.c         # BSSN evolution system
│   ├── nr_initial.c      # Schwarzschild, Kerr, BL, BY
│   ├── nr_horizon.c      # AH finder
│   ├── nr_wave.c         # GW extraction & analysis
│   ├── nr_utils.c        # FD, RK4, SOR, interpolation
│   └── nr_gauge.c        # 1+log, Gamma-driver
├── tests/                # 3 test suites (11 tests)
├── examples/             # 3 runnable examples
├── docs/                 # Knowledge documentation
└── notebooks/            # (reserved)
```

## Quick Start

```bash
make          # Build library + tests + examples
make test     # Run all tests (11/11 pass)
make count    # Count lines: 7010 in include/ + src/
./examples/example_schwarzschild  # Schwarzschild initial data
./examples/example_puncture       # Binary BH puncture data
./examples/example_wave_extract   # GW extraction pipeline
```

## References

- Baumgarte & Shapiro, *Numerical Relativity: Solving Einstein's Equations on the Computer* (Cambridge, 2010)
- Alcubierre, *Introduction to 3+1 Numerical Relativity* (Oxford, 2008)
- Wald, *General Relativity* (Chicago, 1984), Appendix E
- Gourgoulhon, *3+1 Formalism in General Relativity* (Springer, 2012)
- Shibata & Nakamura, PRD 52, 5428 (1995) — BSSN formulation
- Brandt & Brügmann, PRL 78, 3606 (1997) — Puncture method
- Gundlach, PRD 57, 3480 (1998) — Fast flow AH finder
- Baker, Campanelli, Lousto, PRD 65, 044001 (2002) — GW extraction
- Abbott et al. (LIGO/Virgo), PRL 116, 061102 (2016) — GW150914



