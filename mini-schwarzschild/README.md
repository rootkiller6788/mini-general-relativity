# mini-schwarzschild - Schwarzschild Solution

> Reference: MIT 8.962 / Wald Ch.6 / Carroll Ch.5 / MTW Ch.23

## Module Status: COMPLETE

| Level | Status |
|-------|--------|
| L1 Definitions | Complete (10 items) |
| L2 Core Concepts | Complete (7 items) |
| L3 Math Structures | Complete (10 items) |
| L4 Fundamental Laws | Complete (5 items) |
| L5 Computational Methods | Complete (8 items) |
| L6 Canonical Systems | Complete (15 items) |
| L7 Applications | Complete (6 items) |
| L8 Advanced Topics | Complete (10 items) |
| L9 Research Frontiers | Partial (4 items documented) |

**Score: 17/18 - COMPLETE**

**Line Count: include/ + src/ = 3829 >= 3000**

## Core Definitions

- Schwarzschild metric: ds^2 = -(1-rs/r)dt^2 + (1-rs/r)^{-1}dr^2 + r^2 dOmega^2
- Schwarzschild radius: rs = 2GM/c^2
- Event horizon at r = rs
- Photon sphere at r = 1.5 rs
- ISCO at r = 3 rs
- Coordinate systems: standard, EF, KS, PG, isotropic, tortoise, Lemaitre

## Core Theorems

- Birkhoff theorem: Schwarzschild is the unique spherically symmetric vacuum solution
- Ricci tensor vanishes: R_{mu nu} = 0 in vacuum
- Kretschmann finite at horizon: coordinate singularity, not curvature
- First law of BH thermodynamics: dM = T_H dS
- Second law of BH mechanics: horizon area never decreases classically
- Hawking temperature: T_H = hbar c^3/(8 pi G M k_B)
- Photon capture: b_crit = 3 sqrt(3) m

## Core Algorithms

- RK4 integration of 8D geodesic system
- Adaptive RK45 (Fehlberg 5(4)) with error control
- Secant method for circular orbit finding
- Bisection for radial turning points
- Shooting method for bound orbit parameters
- Null geodesic ray tracing with capture detection
- Cardano formula for cubic turning point equation
- Simpson quadrature for periastron advance

## Classical Problems

- Mercury perihelion precession (43 arcsec/century)
- Solar light deflection (1.75 arcsec, Eddington 1919)
- Shapiro time delay (Earth-Venus radar)
- Black hole shadow (M87*, EHT 2019)
- Hulse-Taylor binary pulsar precession
- GPS relativistic time dilation correction
- Accretion disk flux and temperature (Shakura-Sunyaev)
- Radial infall proper time to horizon crossing

## Nine-School Curriculum Mapping

| School | Course | Topics |
|--------|--------|--------|
| MIT | 8.962 GR | Schwarzschild, Birkhoff, Christoffel, Riemann, geodesics |
| Stanford | PHYSICS 230 | BH thermodynamics, Penrose process, causal structure |
| Cambridge | Part III GR | Exact solutions (Schwarzschild, RN, Kerr) |
| Caltech | Ph 205 GR | Experimental tests, black hole shadow |
| ETH | 402-0891 GR | Orbits, rotating BHs, frame dragging |
| Princeton | PHY 535 GR | Gravitational lensing, accretion disks |
| Berkeley | PHYS 231 GR | Geodesic motion, effective potential |
| Oxford | CMT/GR | Weyl tensor, tidal forces, QNMs |
| Tokyo | GR/Cosmology | Schwarzschild, Kerr, collapse |

## API Reference

| Function | Description |
|----------|------------|
| schwarzschild_metric() | Schwarzschild metric components |
| schwarzschild_radius() | rs = 2GM/c^2 |
| schwarzschild_christoffel() | Christoffel symbols |
| schwarzschild_riemann() | Riemann curvature tensor |
| schwarzschild_kretschmann() | Kretschmann scalar invariant |
| schwarzschild_geodesic_derivs() | 8D geodesic system derivative |
| schwarzschild_photon_sphere_radius() | r_ph = 3m |
| schwarzschild_isco_radius() | r_isco = 6m |
| schwarzschild_time_dilation_factor() | sqrt(1-rs/r) |
| schwarzschild_redshift_at_infinity() | Gravitational redshift z |
| schwarzschild_light_deflection_angle() | alpha = 4GM/(bc^2) |
| schwarzschild_perihelion_precession_per_orbit() | Delta phi per orbit |
| schwarzschild_shapiro_time_delay() | Radar delay |
| schwarzschild_hawking_temperature() | T_H |
| schwarzschild_bekenstein_hawking_entropy() | S_BH |
| kerr_metric_covariant() | Kerr metric |
| kerr_penrose_max_efficiency() | Penrose process efficiency |
| schwarzschild_qnm_fundamental_l2_n0() | QNM ringdown |

## Build and Test



## References

- Wald, R.M. (1984) General Relativity, University of Chicago Press
- Carroll, S. (2004) Spacetime and Geometry, Addison-Wesley
- Misner, Thorne, Wheeler (1973) Gravitation, W.H. Freeman
- Chandrasekhar, S. (1983) The Mathematical Theory of Black Holes, Oxford
- Poisson, E. (2004) A Relativist's Toolkit, Cambridge
- Will, C.M. (2018) Theory and Experiment in Gravitational Physics, Cambridge
