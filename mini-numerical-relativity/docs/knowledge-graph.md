# Knowledge Graph — mini-numerical-relativity

## L1: Definitions (Complete)

| Entry | Implementation | Status |
|-------|---------------|--------|
| 3+1 foliation (Σ_t, α, β^i) | `nr_adm.h` ADM state struct | Complete |
| 3-metric γ_{ij} | `nr_grid.h` nr_sym_tensor3_t | Complete |
| Extrinsic curvature K_{ij} | `nr_adm.h` nr_adm_state_t | Complete |
| Lapse function α | `nr_adm.h` nr_gf_t | Complete |
| Shift vector β^i | `nr_grid.h` nr_vector3_t | Complete |
| BSSN conformal factor φ | `nr_bssn.h` phi field | Complete |
| Conformal metric γ̃_{ij} | `nr_bssn.h` gt tensor | Complete |
| Trace-free curvature Ã_{ij} | `nr_bssn.h` At tensor | Complete |
| Conformal connection Γ̃^i | `nr_bssn.h` Gt vector | Complete |
| Weyl scalar Ψ_4 | `nr_wave.h` psi4_modes_t | Complete |
| Gravitational wave strain h_+, h_× | `nr_wave.h` strain_series_t | Complete |
| Apparent horizon / MOTS | `nr_horizon.h` nr_horizon_t | Complete |
| Grid function / ghost zones | `nr_grid.h` nr_gf_t | Complete |
| Complex number type | `nr_grid.h` nr_complex_t | Complete |

## L2: Core Concepts (Complete)

| Concept | Implementation | Status |
|---------|---------------|--------|
| Spacetime foliation | `nr_adm_rhs()` evolution | Complete |
| Conformal decomposition | `nr_bssn_from_adm()` | Complete |
| Constraint addition / damping | `nr_bssn_enforce_constraints()` | Complete |
| Gauge freedom (slicing + shift) | `nr_gauge.c` 1+log, Gamma-driver | Complete |
| Newman-Penrose formalism | `nr_wave_tetrad_m()`, `nr_wave_psi4_point()` | Complete |
| Wave extraction at finite radius | `nr_wave_psi4_sphere()` | Complete |
| ADM energy (mass at infinity) | `nr_adm_mass()`, `nr_adm_komar_mass()` | Complete |
| Puncture method for BHs | `nr_initial.c` Brill-Lindquist | Complete |
| Trumpet slice (stationary 1+log) | `nr_gauge_set_trumpet()` | Complete |

## L3: Mathematical Structures (Complete)

| Structure | Implementation | Status |
|-----------|---------------|--------|
| 3-Riemannian geometry | `nr_adm_christoffel()`, `nr_adm_ricci()` | Complete |
| Covariant derivative D_i on Σ_t | `nr_adm_rhs_point()` D_i terms | Complete |
| Conformal transformations | `nr_bssn_from_adm()` / `nr_adm_from_bssn()` | Complete |
| Trace-free decomposition | `nr_bssn_enforce_constraints()` tr(Ã)=0 | Complete |
| Spin-weighted spherical harmonics | `nr_wave_spin_weighted_Ylm()` | Complete |
| Legendre polynomials | `nr_legendre_P()` | Complete |
| Spherical harmonics | `nr_spherical_harmonic()` | Complete |
| 3D symmetric tensor algebra | `nr_grid.h` nr_sym_tensor3_t operations | Complete |

## L4: Fundamental Laws (Complete)

| Law / Equation | Implementation | Status |
|----------------|---------------|--------|
| Einstein field equations (3+1) | `nr_adm_rhs()` | Complete |
| Hamiltonian constraint | `nr_adm_hamiltonian_constraint()`, `nr_bssn_hamiltonian_constraint()` | Complete |
| Momentum constraint | `nr_adm_momentum_constraint()`, `nr_bssn_momentum_constraint()` | Complete |
| ADM evolution equations | `nr_adm_rhs_point()` | Complete |
| BSSN evolution system | `nr_bssn_rhs_point()` | Complete |
| 1+log slicing equation | `nr_gauge_1pluslog_point()` | Complete |
| Gamma-driver shift condition | `nr_gauge_gamma_driver_point()` | Complete |
| BSSN Gamma constraint | `nr_bssn_gamma_constraint()` | Complete |

## L5: Computational Methods (Complete)

| Method | Implementation | Status |
|--------|---------------|--------|
| Finite difference (2nd/4th/6th/8th order) | `nr_fd_deriv1()`, `nr_fd_deriv2()` | Complete |
| Mixed derivatives (FD) | `nr_fd_deriv2_mixed()` | Complete |
| Flat-space Laplacian (FD) | `nr_fd_laplacian()` | Complete |
| Upwind advection derivative | `nr_fd_upwind()` | Complete |
| RK4 time integration | `nr_rk4_step()` | Complete |
| RK2 (midpoint) integration | `nr_rk2_step()` | Complete |
| Low-storage RK3 (Williamson) | `nr_rk3_low_storage_step()` | Complete |
| CFL time step condition | `nr_cfl_timestep()` | Complete |
| Trilinear interpolation (3D) | `nr_interp_trilinear()` | Complete |
| Tensor/vector interpolation | `nr_interp_tensor_trilinear()` | Complete |
| SOR elliptic solver (3D) | `nr_sor_iteration_3d()`, `nr_sor_solve_3d()` | Complete |
| Kreiss-Oliger dissipation | `nr_kreiss_oliger_gf()` | Complete |
| Fixed-frequency integration | `nr_wave_integrate_psi4()` | Complete |
| Mode decomposition (SWSH) | `nr_wave_decompose_psi4()` | Complete |

## L6: Canonical Systems (Complete)

| System | Implementation | Status |
|--------|---------------|--------|
| Schwarzschild isotropic | `nr_init_schwarzschild_isotropic()` | Complete |
| Schwarzschild trumpet | `nr_init_schwarzschild_trumpet()` | Complete |
| Kerr-Schild initial data | `nr_init_kerr_schild()`, `nr_init_kerr_schild_bssn()` | Complete |
| Brill-Lindquist (multiple BH) | `nr_init_brill_lindquist()` | Complete |
| Bowen-York (momentum + spin) | `nr_init_bowen_york()` | Complete |
| Apparent horizon (Schwarzschild) | `nr_horizon_find_schwarzschild()` | Complete |
| Horizon analytic formulas | `nr_horizon_schwarzschild_radius()`, `nr_horizon_kerr_radius()` | Complete |
| Harmonic slicing | `nr_gauge_harmonic_point()` | Complete |
| Bona-Masso slicing family | `nr_gauge_bona_masso_f()` | Complete |
| Minkowski flat spacetime | `nr_adm_set_minkowski()`, `nr_bssn_set_minkowski()` | Complete |

## L7: Applications (Complete — 3+ applications)

| Application | Implementation | Status |
|-------------|---------------|--------|
| Schwarzschild BH constraint test | `examples/example_schwarzschild.c` | Complete |
| Binary BH puncture data | `examples/example_puncture.c` | Complete |
| GW extraction pipeline | `examples/example_wave_extract.c` | Complete |
| LIGO/Virgo parameter estimation | `nr_wave_ligo_analysis()` | Complete |
| Radiated energy computation | `nr_wave_radiated_energy()` | Complete |
| Ringdown analysis (QNM) | `nr_wave_final_properties()` | Complete |

## L8: Advanced Topics (Partial)

| Topic | Implementation | Status |
|-------|---------------|--------|
| Flow method horizon finder | `nr_horizon_find_flow()` | Complete |
| Conformal constraint solver | `nr_init_solve_psi()` | Complete |
| Puncture equation solver | `nr_init_solve_puncture_u()` | Complete |
| Horizon area / irreducible mass | `nr_horizon_area()` | Complete |
| CCZ4 formulation | Not implemented | Missing |
| Adaptive mesh refinement | Not implemented | Missing |

## L9: Research Frontiers (Partial — documented)

| Topic | Status |
|-------|--------|
| Binary neutron star mergers | Documented |
| Beyond-GR theories | Documented |
| Exascale NR infrastructure | Documented |
| Machine learning waveform models | Documented |
