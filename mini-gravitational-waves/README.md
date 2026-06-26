# mini-gravitational-waves — Gravitational Waves

> **Mini-Pure-Physics** · General Relativity Sub-module
> Reference: MIT 8.962 · Maggiore Vol.1-2 · Wald Ch.4

## Module Status: COMPLETE ✅

- **L1-L6**: Complete (all core definitions, concepts, math, laws, algorithms, systems)
- **L7**: Complete (4 applications: GW150914, PTA, noise models, detection statistics)
- **L8**: Partial (5/8 advanced topics: higher modes, IMR, eccentric binaries, r-modes, memory stacking)
- **L9**: Partial (5 topics documented: multi-messenger, primordial GWs, cosmic strings, testing GR)

**Score: 16/18** · **include/ + src/ lines: 3492 ≥ 3000** ✅

---

## API Reference

### Core (`gw_core.h`)
| Function | Description |
|----------|-------------|
| `gw_tensor_zero/add/sub/scale` | 3×3 symmetric tensor algebra |
| `gw_tensor_trace` | Trace Tr[T] |
| `gw_tensor_contract` | Double contraction A:B |
| `gw_tensor_tt_project` | TT gauge projection Λ_{ij,kl} |
| `gw_tensor_eigensystem` | 3×3 symmetric eigenvalue decomposition |
| `gw_antenna_pattern` | LIGO/Virgo antenna response F_+, F_× |
| `gw_antenna_pattern_lisa` | LISA antenna response |
| `gw_pol_tensor_plus/cross` | e⁺ and e× polarization basis tensors |
| `gw_pol_rotate` | Polarization rotation by angle ψ |
| `gw_strain_decompose` | h_{ij}^{TT} → h_+, h_× |
| `gw_pn_parameter` | PN expansion parameter x |
| `gw_orbital_separation` | a from Kepler III |
| `gw_isco_frequency` | f_ISCO for Schwarzschild |

### Quadrupole (`gw_quadrupole.h`)
| Function | Description |
|----------|-------------|
| `gw_quadrupole_moment` | Mass quadrupole I_{ij} for N particles |
| `gw_reduced_quadrupole` | Traceless quadrupole Q_{ij} |
| `gw_current_quadrupole` | Current quadrupole J_{ij} |
| `gw_mass_octupole` | Mass octupole (L=3) |
| `gw_quadrupole_strain` | Einstein quadrupole formula h = (2G/c⁴r) Q̈ |
| `gw_quadrupole_strain_from_state` | Strain from particle state (x,v,a) |
| `gw_quadrupole_luminosity` | GW power L_GW = (G/5c⁵)⟨Q_dddot:Q_dddot⟩ |
| `gw_binary_quadrupole_strain` | Circular binary strain h_+, h_× |
| `gw_binary_luminosity_circular` | Circular binary L_GW |
| `gw_energy_density` | Isaacson ρ_GW |

### Binary Systems (`gw_binary.h`)
| Function | Description |
|----------|-------------|
| `gw_chirp_mass` | M_c = (m₁m₂)^{3/5}/(m₁+m₂)^{1/5} |
| `gw_symmetric_mass_ratio` | η = m₁m₂/(m₁+m₂)² |
| `gw_reduced_mass` | μ = m₁m₂/(m₁+m₂) |
| `gw_frequency_derivative` | df/dt (Newtonian chirp) |
| `gw_time_to_coalescence` | τ(f) time to merger |
| `gw_frequency_at_time` | f(Δt) before merger |
| `gw_number_of_cycles` | N_cyc in frequency band |
| `gw_taylor_t2_phase` | TaylorT2 phase Φ(f) |
| `gw_taylor_t4_rhs` | TaylorT4 ODE (df/dt, dΦ/dt) |
| `gw_taylor_f2_phase` | TaylorF2 frequency-domain phase |
| `gw_orbital_evolution_rhs` | Peters eccentric ODE |
| `gw_eccentricity_enhancement_da/de` | f(e), g(e) factors |
| `gw_params_gw150914/bns/nsbh` | Reference system initializers |

### Detection (`gw_detection.h`)
| Function | Description |
|----------|-------------|
| `gw_inner_product` | (a|b) = 4Re∫ a*(f)b(f)/S_n(f) df |
| `gw_snr` | SNR = (d|h)/√(h|h) |
| `gw_optimal_snr` | ρ_opt = √(h|h) |
| `gw_fisher_matrix` | Γ_{ij} = (∂_i h|∂_j h) |
| `gw_fisher_uncertainties` | σ_i = √[(Γ⁻¹)_{ii}] |
| `gw_fisher_covariance` | Cov = Γ⁻¹ |
| `gw_matrix_invert_cholesky` | SPD matrix inversion |
| `gw_cholesky_solve` | Ax=b for SPD A |
| `gw_matched_filter_td` | Time-domain cross-correlation |
| `gw_psd_aligo_design/o3` | LIGO noise PSD |
| `gw_psd_lisa` | LISA noise PSD |
| `gw_psd_et` | Einstein Telescope PSD |
| `gw_horizon_distance` | D_L for SNR=8 |
| `gw_false_alarm_probability/rate` | FAP, FAR |

### Waveforms (`gw_waveforms.h`)
| Function | Description |
|----------|-------------|
| `gw_waveform_inspiral_td` | Time-domain inspiral (TaylorT4) |
| `gw_waveform_fd_spa` | Frequency-domain SPA waveform |
| `gw_waveform_imr_fd` | IMRPhenomD-like FD waveform |
| `gw_waveform_ringdown` | Damped sinusoid ringdown |
| `gw_qnm_220_mode` | Kerr QNM (l=m=2, n=0) |
| `gw_final_mass_spin` | M_f, a_f from m1, m2 |
| `gw_higher_mode_amplitude` | (l,m) harmonic amplitude |
| `gw_waveform_higher_modes` | Multi-mode FD waveform |
| `gw_strain_amplitude_fd` | |h(f)| Newtonian amplitude |
| `gw_characteristic_strain` | h_c(f) = 2f|h(f)| |

### Memory (`gw_memory.h`)
| Function | Description |
|----------|-------------|
| `gw_linear_memory` | Δh from unbound sources |
| `gw_nonlinear_memory` | Christodoulou memory |
| `gw_cbc_memory` | CBC memory amplitude |
| `gw_memory_growth` | Memory vs. frequency |
| `gw_memory_snr_aligo` | Memory SNR in aLIGO |
| `gw_memory_events_for_detection` | N_events for SNR=3 |

### Stochastic Background (`gw_stochastic.h`)
| Function | Description |
|----------|-------------|
| `gw_omega_gw` | Ω_GW(f) power-law model |
| `gw_critical_density` | ρ_c = 3H₀²/8πG |
| `gw_omega_cbc_background` | CBC astrophysical background |
| `gw_omega_inflation` | Inflationary tensor background |
| `gw_omega_phase_transition` | First-order PT spectrum |
| `gw_omega_cosmic_strings` | Nambu-Goto string background |
| `gw_overlap_reduction_hl` | γ(f) for LIGO H-L |
| `gw_stochastic_snr` | Cross-correlation SNR |

### Pulsars (`gw_pulsar.h`)
| Function | Description |
|----------|-------------|
| `gw_cw_strain` | h₀ = (4π²G/c⁴) ε I f²/D |
| `gw_cw_frequency` | f_gw = 2 f_rot |
| `gw_spindown_limit` | h₀^{sd} from f_dot |
| `gw_spindown_limit_age` | Age-based spin-down limit |
| `gw_rmode_strain` | R-mode GW strain |
| `gw_rmode_frequency` | f_gw = 4/3 f_rot |
| `gw_pta_timing_residual` | Earth-term timing residual |
| `gw_hellings_downs` | HD correlation curve C(ξ) |
| `gw_f_statistic` | JKS F-statistic |
| `gw_fstat_fap` | FAP for F-statistic |
| `gw_pulsar_antenna_pattern` | PTA antenna patterns |

---

## Core Theorems

### Einstein Quadrupole Formula (1918)
```
h_{ij}^{TT}(t, r̂) = (2G / c⁴r) Λ_{ij,kl}(r̂) d²Q_{kl}/dt²(t - r/c)
```
where Λ_{ij,kl} = P_{ik}P_{jl} - ½P_{ij}P_{kl} is the TT projector and
P_{ij} = δ_{ij} - n_i n_j is the transverse projector.

### Quadrupole Luminosity
```
L_GW = (G / 5c⁵) ⟨d³Q_{ij}/dt³ · d³Q^{ij}/dt³⟩
```

### Newtonian Chirp
```
df/dt = (96/5) π^{8/3} (G M_c / c³)^{5/3} f^{11/3}
τ(f) = (5/256) (πf)^{-8/3} (G M_c / c³)^{-5/3}
```

### GW150914 Reference Values
```
m₁ = 36 M⊙, m₂ = 29 M⊙ → M_c ≈ 28 M⊙
D_L ≈ 410 Mpc
f_start = 20 Hz → τ ≈ 0.85 s to merger
h_peak ≈ 1 × 10^{-21} at ISCO
```

---

## Physical Constants
| Symbol | Value | Description |
|--------|-------|-------------|
| G | 6.67430×10⁻¹¹ m³/(kg·s²) | Gravitational constant |
| c | 2.99792458×10⁸ m/s | Speed of light |
| M⊙ | 1.98847×10³⁰ kg | Solar mass |
| Mpc | 3.0857×10²² m | Megaparsec |

---

## Directory Structure
```
mini-gravitational-waves/
├── Makefile              # make test, make examples
├── README.md             # This file
├── include/              # 8 header files
│   ├── gw_core.h         #   Constants, tensors, polarization, data types
│   ├── gw_quadrupole.h    #   Quadrupole formula, luminosity
│   ├── gw_binary.h       #   Binary inspiral, chirp, PN phase
│   ├── gw_detection.h    #   Matched filtering, Fisher matrix, PSD
│   ├── gw_waveforms.h    #   Waveform generation (TD/FD)
│   ├── gw_memory.h       #   Linear + nonlinear GW memory
│   ├── gw_stochastic.h   #   SGWB: Omega_GW, overlap reduction
│   └── gw_pulsar.h       #   Continuous waves, PTA, F-statistic
├── src/                  # 9 implementation files
│   ├── gw_core.c         #   Tensor ops, TT projection, antenna patterns
│   ├── gw_quadrupole.c    #   Quadrupole moments, strain, luminosity
│   ├── gw_binary.c       #   Binary dynamics, PN phase, eccentricity
│   ├── gw_detection.c    #   Inner product, SNR, Fisher, Cholesky, PSD
│   ├── gw_waveforms.c    #   SPA, IMR, ringdown, higher modes
│   ├── gw_memory.c       #   Memory effects, detectability
│   ├── gw_stochastic.c   #   SGWB models, overlap reduction, SNR
│   ├── gw_pulsar.c       #   CW strain, spin-down, HD correlation
│   └── gw_formal.lean    #   Lean 4 formalization (theorems)
├── tests/                # 3 test files (core, quadrupole, binary)
├── examples/             # 3 example programs
│   ├── example_gw150914.c
│   ├── example_compact_binary_inspiral.c
│   └── example_pulsar_timing.c
└── docs/                 # 5 knowledge documents
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-alignment.md
    └── course-tree.md
```

---

## Build & Test

```bash
# Build all objects
make

# Run all tests
make test

# Run all examples
make examples

# Count lines of code
make count
```

---

## 九层知识覆盖摘要

| Level | 名称 | 状态 | 条目 |
|-------|------|------|------|
| L1 | 定义 | ✅ Complete | 15 个核心定义 |
| L2 | 核心概念 | ✅ Complete | 15 个核心概念 |
| L3 | 数学结构 | ✅ Complete | 8 个数学结构 |
| L4 | 基本定律 | ✅ Complete | 9 条基本定律 |
| L5 | 计算方法 | ✅ Complete | 10 个算法 |
| L6 | 经典系统 | ✅ Complete | 9 个经典系统 |
| L7 | 应用 | ✅ Complete | 4 个应用 |
| L8 | 进阶主题 | ⚠️ Partial | 5/8 实现 |
| L9 | 研究前沿 | ⚠️ Partial | 5 记录, 0 实现 |

---

## 九校课程映射

| 学校 | 关键课程 | 对应实现 |
|------|---------|---------|
| MIT | 8.962 GR | TT gauge, quadrupole formula |
| Caltech | Ph 205 GR | GW generation, binary inspiral |
| Princeton | PHY 535 GR | GW energy-momentum, GW astronomy |
| Cambridge | Part III GR | Linearized theory, detection |
| Oxford | GR Graduate | Memory effect, stochastic backgrounds |
| Stanford | PHYSICS 230 | Polarizations, PN expansion |
| Berkeley | PHYS 231 | GW generation, detection, data analysis |
| ETH | 402-0891 | GW experiments, parameter estimation |
| 东京大学 | GW Physics | Compact binaries, PTA, GW cosmology |

---

## 参考教材

- Wald, R.M. (1984) *General Relativity*, Ch.4
- Maggiore, M. (2008) *Gravitational Waves Vol.1: Theory and Experiments*
- Maggiore, M. (2018) *Gravitational Waves Vol.2: Astrophysics and Cosmology*
- Blanchet, L. (2014) *Gravitational Radiation from Post-Newtonian Sources* (Living Rev. Relativity)
- Abbott, B.P. et al. (2016) *Observation of Gravitational Waves from a Binary Black Hole Merger* (PRL 116, 061102)
- LIGO Scientific Collaboration & Virgo Collaboration
