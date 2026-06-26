# Knowledge Graph - mini-schwarzschild

## L1 - Definitions

| Entry | Implementation |
|-------|---------------|
| Schwarzschild metric | schwarzschild_g_tt/g_rr/g_theta_theta/g_phi_phi |
| Schwarzschild radius | schwarzschild_radius() |
| Geometric mass | geometric_mass() |
| Coordinate types | SchwarzschildCoordType enum (standard, EF, KS, PG, isotropic, tortoise) |
| Spacetime point | SchwarzschildPoint struct |
| 4-vector | SchwarzschildVector struct |
| 4x4 symmetric tensor | SchwarzschildSymmetric4x4 struct |
| BH parameters | SchwarzschildBlackHole struct |
| Tetrad basis | schwarzschild_tetrad() |
| Lapse/shift | schwarzschild_lapse(), schwarzschild_shift() |

## L2 - Core Concepts

| Entry | Implementation |
|-------|---------------|
| Event horizon | schwarzschild_is_inside_horizon() |
| Coordinate singularity | Kretschmann finite at horizon |
| Time dilation | schwarzschild_time_dilation_factor() |
| Gravitational redshift | schwarzschild_redshift_at_infinity/finite() |
| Coordinate transformations | tortoise, isotropic, KS, EF, PG |
| Proper distance | schwarzschild_proper_distance() |

## L3 - Mathematical Structures

| Entry | Implementation |
|-------|---------------|
| Christoffel symbols | SchwarzschildChristoffel struct, 13 non-zero components |
| Riemann tensor | SchwarzschildRiemann struct, 6 independent |
| Ricci tensor (=0 in vacuum) | schwarzschild_ricci() |
| Einstein tensor (=0) | schwarzschild_einstein() |
| Weyl tensor | SchwarzschildWeyl struct |
| Kretschmann scalar | schwarzschild_kretschmann() = 12rs^2/r^6 |
| Chern-Pontryagin | schwarzschild_chern_pontryagin() = 0 |
| Tidal tensor | schwarzschild_tidal_tensor() |
| Geodesic deviation | schwarzschild_geodesic_deviation() |
| Bel-Robinson energy | schwarzschild_bel_robinson_super_energy() |

## L4 - Fundamental Laws

| Entry | Implementation |
|-------|---------------|
| Vacuum Einstein eqns | Ricci=0, Einstein=0 verified |
| Birkhoff theorem | theorem birkhoff_theorem (Lean) |
| Geodesic equation | schwarzschild_geodesic_derivs() (8D Hamiltonian) |
| First law BH thermo | schwarzschild_first_law_dM_dA() |
| Second law BH mechanics | schwarzschild_area_increase() |

## L5 - Computational Methods

| Entry | Implementation |
|-------|---------------|
| RK4 integration | schwarzschild_rk4_step_geodesic() |
| Adaptive RK45 (Fehlberg) | schwarzschild_rk45_adaptive_step() |
| Secant method | schwarzschild_find_circular_orbit_secant() |
| Bisection | schwarzschild_find_turning_points_bisection() |
| Shooting method | schwarzschild_shooting_bound_orbit() |
| Ray tracing (null geodesics) | schwarzschild_trace_null_geodesic() |
| Cubic Hermite interpolation | schwarzschild_cubic_hermite_interp() |
| Simpson quadrature | schwarzschild_periastron_advance_numerical() |

## L6 - Canonical Systems

| Entry | Implementation |
|-------|---------------|
| Photon sphere (r_ph=3m) | schwarzschild_photon_sphere_radius() |
| ISCO (r_isco=6m) | schwarzschild_isco_radius() |
| Marginally bound orbit | schwarzschild_marginally_bound_radius() |
| Effective potential | schwarzschild_V_eff_massive/massless() |
| Circular orbits | schwarzschild_L/E_for_circular_orbit() |
| Orbital periods | schwarzschild_orbital_period_keplerian/relativistic() |
| Perihelion precession | schwarzschild_perihelion_precession_per_orbit() |
| Light deflection | schwarzschild_light_deflection_angle() |
| Shapiro delay | schwarzschild_shapiro_time_delay() |
| Einstein ring | schwarzschild_einstein_ring_radius() |
| Gravitational lensing | schwarzschild_lens_magnification() |
| Black hole shadow | schwarzschild_shadow_angular_radius() |
| Accretion disk | schwarzschild_disk_flux/temperature() |
| Orbit classification | schwarzschild_classify_orbit() |
| Radial infall | schwarzschild_proper_time_radial_infall() |

## L7 - Applications (Complete, 6 items)

| Entry | Implementation |
|-------|---------------|
| GPS time dilation | schwarzschild_gps_time_dilation() |
| Mercury precession (43 arcsec/century) | Tested numerically |
| Solar light deflection (1.75 arcsec) | Tested numerically |
| M87* shadow (EHT 2019) | m87_shadow_radius_uas |
| Shapiro delay (Earth-Venus) | schwarzschild_shapiro_earth_venus() |
| Hulse-Taylor binary precession | schwarzschild_hulse_taylor_precession() |

## L8 - Advanced Topics (Complete, 10 items)

| Entry | Implementation |
|-------|---------------|
| Reissner-Nordstrom metric | RN BH with charge, horizons, extremal limit |
| Kerr metric | Kerr BH with spin, ergosphere, ISCO pro/retro |
| Hawking temperature | T_H = hbar*c^3/(8*pi*G*M*k_B) |
| Bekenstein-Hawking entropy | S = k_B*A*c^3/(4*G*hbar) |
| BH evaporation | t_evap ~ M^3, Page (1976) |
| Penrose process | Maximum 29.3% efficiency for extremal Kerr |
| Blandford-Znajek | EM power extraction from Kerr BH |
| Frame dragging | Lense-Thirring precession, kerr_frame_dragging_omega() |
| Quasinormal modes | Ringdown frequency and damping time |
| Kerr ISCO | Bardeen et al. (1972) formula |

## L9 - Research Frontiers (Partial, 4 items documented)

| Entry | Status |
|-------|--------|
| Black hole information paradox | Documented, BH thermodynamics API |
| Naked singularities | RN detection: rn_is_naked_singularity |
| Quantum gravity (Planck scale) | Constants SCHW_M_PL, SCHW_L_PL defined |
| Gravitational wave astronomy | QNM ringdown implemented |
