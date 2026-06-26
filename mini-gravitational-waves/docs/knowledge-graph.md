# Gravitational Waves — Knowledge Graph

## L1 — Definitions (Complete)

| # | Definition | C Implementation | Lean Formalization |
|---|-----------|-----------------|-------------------|
| 1 | Physical constants (G, c, Msun, Mpc) | `gw_core.h` — `GW_G`, `GW_C`, `GW_MSUN`, `GW_MPC` | — |
| 2 | 3×3 symmetric tensor GwTensor3 | `gw_core.h` — `typedef struct GwTensor3` | `structure GwTensor3` |
| 3 | GW polarization modes (+, ×, scalar, vector) | `gw_core.h` — `enum GwPolarization` | `inductive GwPolarization` |
| 4 | Binary source parameters | `gw_core.h` — `typedef struct GwBinaryParams` | — |
| 5 | Detector response (antenna patterns) | `gw_core.h` — `GwDetectorResponse` | — |
| 6 | GW strain time series | `gw_core.h` — `GwStrainSeries` | — |
| 7 | GW frequency series | `gw_core.h` — `GwFreqSeries` | — |
| 8 | Stochastic background parameters | `gw_core.h` — `GwStochasticParams` | — |
| 9 | Quadrupole moment I_{ij} | `gw_quadrupole.h` — `gw_quadrupole_moment()` | `structure QuadrupoleMoment` |
| 10 | Reduced quadrupole Q_{ij} | `gw_quadrupole.h` — `gw_reduced_quadrupole()` | `reduced_quadrupole()` |
| 11 | Chirp mass M_c | `gw_binary.h` — `gw_chirp_mass()` | `chirp_mass()` |
| 12 | Symmetric mass ratio η | `gw_binary.h` — `gw_symmetric_mass_ratio()` | — |
| 13 | Plus/cross polarization tensors | `gw_core.h` — `gw_pol_tensor_plus/cross()` | `e_plus`, `e_cross` |
| 14 | GW memory tensor Δh_{ij} | `gw_memory.h` — `gw_linear_memory()` | `GwMemory` |
| 15 | Omega_GW stochastic energy density | `gw_stochastic.h` — `gw_omega_gw()` | — |

## L2 — Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Linearized gravity on Minkowski background | Background of all GW theory |
| 2 | TT gauge projection Λ_{ij,kl} | `gw_tensor_tt_project()` |
| 3 | Quadrupole GW strain formula | `gw_quadrupole_strain()` |
| 4 | GW luminosity (quadrupole power) | `gw_quadrupole_luminosity()` |
| 5 | Isaacson stress-energy tensor | `gw_energy_density()` |
| 6 | Antenna pattern functions (ground) | `gw_antenna_pattern()` |
| 7 | Antenna pattern functions (space/LISA) | `gw_antenna_pattern_lisa()` |
| 8 | Frequency evolution df/dt | `gw_frequency_derivative()` |
| 9 | Time to coalescence | `gw_time_to_coalescence()` |
| 10 | Noise-weighted inner product | `gw_inner_product()` |
| 11 | Signal-to-noise ratio (SNR) | `gw_snr()`, `gw_optimal_snr()` |
| 12 | GW memory (linear + nonlinear) | `gw_linear_memory()`, `gw_nonlinear_memory()` |
| 13 | Stochastic GW background energy | `gw_omega_gw()` |
| 14 | Continuous GW strain from NS | `gw_cw_strain()` |
| 15 | Polarization decomposition h→h_+,h_× | `gw_strain_decompose()` |

## L3 — Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|----------|---------------|
| 1 | 3×3 symmetric tensor algebra | `gw_tensor_zero/add/sub/scale/contract()` |
| 2 | TT projection operator | `gw_tensor_tt_project()` |
| 3 | 3×3 eigenvalue decomposition | `gw_tensor_eigensystem()` |
| 4 | Fisher information matrix | `gw_fisher_matrix()`, `gw_fisher_covariance()` |
| 5 | Cholesky decomposition | `gw_matrix_invert_cholesky()`, `gw_cholesky_solve()` |
| 6 | PN expansion coefficients | `gw_pn_phase_coeffs()` (static) |
| 7 | Hellings-Downs correlation function | `gw_hellings_downs()` |
| 8 | PTA timing residual formalism | `gw_pta_timing_residual()` |

## L4 — Fundamental Laws (Complete)

| # | Law / Theorem | Formula | Implementation |
|---|--------------|---------|---------------|
| 1 | Einstein quadrupole formula | h_{ij}^{TT} = (2G/c⁴r) d²Q_{ij}/dt² | `gw_quadrupole_strain()` |
| 2 | Quadrupole luminosity | L_GW = (G/5c⁵) ⟨Q_dddot : Q_dddot⟩ | `gw_quadrupole_luminosity()` |
| 3 | Binary inspiral df/dt (Newtonian) | df/dt = (96/5) π^{8/3} (GMc/c³)^{5/3} f^{11/3} | `gw_frequency_derivative()` |
| 4 | Kepler III + GW balance | a = (GM)^{1/3} / (πf)^{2/3} | `gw_orbital_separation()` |
| 5 | Isaacson GW energy density | ρ_GW = (c²/32πG) ⟨ḣ_{ij} ḣ^{ij}⟩ | `gw_energy_density()` |
| 6 | Spin-down limit | h_0^{sd} ∝ √(I|ḟ| / D²f) | `gw_spindown_limit()` |
| 7 | Memory effect formula | Δh ∼ (G/c⁴r) ΔE_GW | `gw_nonlinear_memory()` |
| 8 | GW strain from ellipticity | h_0 = (4π²G/c⁴) ε I f² / D | `gw_cw_strain()` |
| 9 | TaylorF2 PN phase (up to 3.5PN) | Ψ(f) = 2πf t_c - φ_c - π/4 + ... | `gw_taylor_f2_phase()` |

## L5 — Computational Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | TaylorT2 PN phase evaluation | `gw_taylor_t2_phase()` |
| 2 | TaylorT4 ODE RHS (time domain) | `gw_taylor_t4_rhs()` |
| 3 | TaylorF2 frequency-domain phase | `gw_taylor_f2_phase()` |
| 4 | Stationary Phase Approximation waveform | `gw_waveform_fd_spa()` |
| 5 | Time-domain inspiral integration | `gw_waveform_inspiral_td()` |
| 6 | Peters orbital evolution ODE | `gw_orbital_evolution_rhs()` |
| 7 | Matched filter (time domain) | `gw_matched_filter_td()` |
| 8 | Fisher matrix (numerical diff) | `gw_fisher_matrix()` |
| 9 | Cholesky-based matrix inversion | `gw_matrix_invert_cholesky()` |
| 10 | Eccentricity enhancement functions | `gw_eccentricity_enhancement_da/de()` |

## L6 — Canonical Systems (Complete)

| # | System | Implementation |
|---|--------|---------------|
| 1 | GW150914-like BBH | `gw_params_gw150914()`, `example_gw150914.c` |
| 2 | BNS (binary neutron star) | `gw_params_bns()`, `example_compact_binary_inspiral.c` |
| 3 | NS-BH binary | `gw_params_nsbh()` |
| 4 | Kerr QNM ringdown | `gw_qnm_220_mode()`, `gw_waveform_ringdown()` |
| 5 | Pulsar timing array (PTA) | `gw_pta_timing_residual()`, `gw_hellings_downs()` |
| 6 | CBC stochastic background | `gw_omega_cbc_background()` |
| 7 | Inflationary SGWB | `gw_omega_inflation()` |
| 8 | Phase transition SGWB | `gw_omega_phase_transition()` |
| 9 | Cosmic string SGWB | `gw_omega_cosmic_strings()` |

## L7 — Applications (Complete, 4 applications)

| # | Application | Implementation |
|---|------------|---------------|
| 1 | GW150914 waveform simulation | `example_gw150914.c` — compares with LIGO detection |
| 2 | PTA GW detection (NANOGrav 15yr) | `example_pulsar_timing.c` — HD correlation |
| 3 | LIGO/Virgo noise PSD models | `gw_psd_aligo_design()`, `gw_psd_aligo_o3()`, `gw_psd_lisa()`, `gw_psd_et()` |
| 4 | Detection statistics (FAP, FAR) | `gw_false_alarm_probability()`, `gw_false_alarm_rate()` |

## L8 — Advanced Topics (Partial, 3 topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Higher waveform harmonics (l>2) | `gw_higher_mode_amplitude()`, `gw_waveform_higher_modes()` |
| 2 | Eccentric binary GW emission | `gw_orbital_evolution_rhs()`, eccentric enhancement functions |
| 3 | R-mode GWs from neutron stars | `gw_rmode_strain()`, `gw_rmode_frequency()` |
| 4 | IMR phenomenological waveforms | `gw_waveform_imr_fd()` |
| 5 | Memory detectability and stacking | `gw_memory_snr_aligo()`, `gw_memory_events_for_detection()` |

## L9 — Research Frontiers (Partial, documented)

| # | Topic | Documentation |
|---|-------|--------------|
| 1 | Multi-messenger astronomy (GW+EM) | Standard sirens: `StandardSiren` in Lean, H₀ constraint |
| 2 | Primordial gravitational waves | `gw_omega_inflation()` — tensor-to-scalar ratio r |
| 3 | Testing GR with GW polarization | `GwPolarization` enum includes beyond-GR modes |
| 4 | GW from cosmic strings | `gw_omega_cosmic_strings()` |
| 5 | GW from phase transitions | `gw_omega_phase_transition()` |
