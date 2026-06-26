# Knowledge Graph — mini-cosmological-model

## L1: Definitions (Complete)
- **Hubble parameter** H(z), H0, h: `CosmoModel.H0`, `friedmann_H_z()`
- **Density parameters** Ω_m, Ω_r, Ω_Λ, Ω_k: `CosmoModel` struct
- **Redshift z**: `CosmoState.z`, `flrw_redshift_from_scale()`
- **Scale factor a(t)**: `CosmoState.a`, `flrw_scale_from_redshift()`
- **Equation of state w**: `CosmoSpecies` enum, `friedmann_w_eff()`
- **Deceleration parameter q**: `flrw_deceleration_parameter()`
- **Jerk parameter j**: `flrw_jerk_parameter()`
- **Snap parameter s**: `flrw_snap_parameter()`
- **Proper vs comoving distance**: `flrw_proper_distance()`
- **Critical density ρ_crit**: `cosmo_rho_critical()`
- **Primordial power spectrum**: `ps_primordial_dimensionless()`
- **Dimensionless power spectrum**: `ps_matter_dimensionless()`
- **Lookback time**: `cosmo_lookback_time()`
- **Light travel distance**: `distances_light_travel()`
- All defined in `include/cosmo_model.h` + `src/cosmology.lean`

## L2: Core Concepts (Complete)
- **Cosmological principle**: FLRW metric in `flrw.h`
- **Cosmic expansion**: `friedmann_dadt()`
- **Big Bang model**: `cosmo_age_of_universe()`
- **Dark energy / Λ**: `CosmoSpecies.lambda`, `friedmann_z_acceleration_onset()`
- **Dark matter**: `CosmoModel.Omega_cdm`
- **CMB**: `CosmoModel.T_cmb0`, `transfer_Cl_SW()`
- **Inflation**: `ps_primordial_dimensionless()`
- **Structure formation**: `perturbations_growth_factor_D()`

## L3: Mathematical Structures (Complete)
- **FLRW metric**: `flrw_g_tt/rt/thetatheta/phiphi()`
- **Ricci scalar**: `flrw_ricci_scalar()`
- **Friedmann ODE system**: `friedmann_E_z()`, `friedmann_dH_da()`
- **Perturbation theory**: Growth equation
- **Spherical Bessel functions**: `perturbations_cmb_Cl_TT_SW()`
- **Power spectrum**: `ps_matter_power_linear()`
- **Transfer function**: `ps_transfer_BBKS()`, `ps_transfer_EH()`
- **Cosmography**: `flrw_jerk_parameter()`, `flrw_snap_parameter()`

## L4: Fundamental Laws (Complete)
- **Friedmann 1**: H² = 8πGρ/3 − kc²/a² + Λc²/3 → `friedmann_E_z()` + Lean
- **Friedmann 2**: ä/a = −4πG(ρ+3P/c²)/3 + Λc²/3 → `friedmann_dH_da()`
- **Continuity**: ρ̇ + 3H(ρ+P/c²) = 0 → `friedmann_rho_a()` + Lean
- **Etherington**: D_L = (1+z)²D_A → C + Lean proof
- **Sum rule**: ΣΩᵢ = 1 → Lean `sum_rule_flat`
- **Poisson**: k²Φ = −(3/2)ΩₘH₀²δ/a → `transfer_potential_to_density()`

## L5: Computational Methods (Complete)
- RK4 integration (single step + adaptive)
- Simpson integration
- Gauss-Legendre quadrature (5-point adaptive)
- Trapezoidal integration
- Spherical Bessel recurrence
- Finite difference derivatives
- Root finding (iterative)
- Convolution methods
- Lookup table generation

## L6: Canonical Systems (Complete)
- Einstein-de Sitter, Radiation-dominated, de Sitter
- Planck 2018 ΛCDM
- Matter-radiation equality, acceleration onset
- Particle horizon, event horizon
- Age of universe, distance modulus
- All distance measures (D_C, D_M, D_A, D_L)
- Comoving volume, sigma_8
- Linear growth factor, BBKS/EH transfer functions
- Correlation function, linear matter power spectrum

## L7: Applications (Complete - 7)
- Age of universe calculator
- Distance-redshift relation (SNIa Hubble diagram)
- Cosmic expansion history visualization
- CMB temperature spectrum (Sachs-Wolfe)
- BAO feature extraction
- Nonlinear power spectrum (Halofit)
- Growth rate fσ₈ (RSD)

## L8: Advanced Topics (Partial - 8)
- Spherical collapse (δ_c)
- Press-Schechter mass function
- Halo bias (peak-background split)
- Nonlinear Halofit corrections
- CMB lensing potential
- Running spectral index
- Neutrino free-streaming
- Reionization modeling

## L9: Research Frontiers (Documented)
- Hubble tension (H0=67.4 vs 73)
- S8 tension (σ₈√(Ωₘ/0.3) = 0.83 vs 0.76)
- Early dark energy
- Primordial gravitational waves
- Modified gravity (wCDM, w₀wₐ)
