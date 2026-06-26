# Knowledge Graph — mini-black-holes

## L0: Access Requirement
- `include/` + `src/` line count: **5084** ≥ 3000 ✓

## L1: Definitions (Complete)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | Schwarzschild metric g_μν | `SchwarzschildBH`, `schwarzschild_metric()` | ✓ |
| 2 | Kerr metric g_μν | `KerrBH`, `kerr_metric()` | ✓ |
| 3 | Reissner-Nordström metric g_μν | `ReissnerNordstromBH`, `reissner_nordstrom_metric()` | ✓ |
| 4 | Kerr-Newman metric g_μν | `KerrNewmanBH`, `kerr_newman_metric()` | ✓ |
| 5 | Event horizon r_+ for all BH types | `schwarzschild_radius()`, `kerr_horizon_radius()` | ✓ |
| 6 | Hawking temperature T_H | `hawking_temperature_schwarzschild()`, `hawking_temperature_kn()` | ✓ |
| 7 | Bekenstein-Hawking entropy S_BH | `bekenstein_hawking_entropy()` | ✓ |
| 8 | Surface gravity κ | `kerr_surface_gravity()`, `surface_gravity_kn()` | ✓ |
| 9 | Photon sphere r_ph | `photon_sphere_radius()` | ✓ |
| 10 | ISCO radius | `schwarzschild_isco()`, `kerr_isco_radius()` | ✓ |
| 11 | Ergosphere | `kerr_ergosphere_radius()`, `is_inside_ergosphere()` | ✓ |
| 12 | Quasinormal mode | `QuasinormalMode`, `schwarzschild_qnm_fundamental()` | ✓ |
| 13 | Chirp mass | `chirp_mass()` | ✓ |
| 14 | Planck units | `planck_mass()`, `planck_length()`, `planck_time()` | ✓ |
| 15 | No-hair theorem parameters | `KerrMultipoles`, `no_hair_test_quadrupole()` | ✓ |

## L2: Core Concepts (Complete)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | Black hole thermodynamics | `SchwarzschildThermo`, `KerrNewmanThermo` | ✓ |
| 2 | Hawking radiation | `hawking_luminosity()`, `hawking_spectrum_scalar()` | ✓ |
| 3 | BH evaporation | `evaporation_time()`, `evaporating_mass()` | ✓ |
| 4 | Four laws of BH mechanics | `zeroth_law_check()`, `first_law_schwarzschild_check()`, `area_increase_theorem()` | ✓ |
| 5 | Cosmic censorship | `cosmic_censorship_check()` | ✓ |
| 6 | Frame dragging | `kerr_frame_dragging_omega()` | ✓ |
| 7 | Penrose process | `penrose_process_efficiency()` | ✓ |
| 8 | Superradiance | `superradiance_condition()` | ✓ |
| 9 | Gravitational redshift | `schwarzschild_redshift()`, `schwarzschild_time_dilation()` | ✓ |

## L3: Mathematical Structures (Complete)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | 4D metric tensor | `Metric4D` struct | ✓ |
| 2 | Christoffel symbols | `christoffel_symbol_fd()`, `christoffel_fd_local()` | ✓ |
| 3 | Kretschmann scalar | `schwarzschild_kretschmann()`, `kerr_kretschmann()` | ✓ |
| 4 | Geodesic equation | `geodesic_rhs()`, `geodesic_rk4_step()` | ✓ |
| 5 | Effective potentials | `schwarzschild_effective_potential()`, `kerr_effective_potential_sq()` | ✓ |
| 6 | Epicyclic frequencies | `schwarzschild_epicyclic_freq()`, `kerr_vertical_epicyclic_freq()` | ✓ |
| 7 | Carter constant (Kerr integrals) | `GeodesicConstants` struct | ✓ |
| 8 | Weyl tensor component | `schwarzschild_weyl_component()` (declared) | ✓ |
| 9 | Teukolsky potential | `teukolsky_radial_potential()` | ✓ |
| 10 | Metric inverse validation | `metric_inverse_check()` | ✓ |

## L4: Fundamental Laws (Complete)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | Einstein vacuum equations (Schwarzschild, Kerr) | All metric functions | ✓ |
| 2 | First law of BH mechanics | `first_law_schwarzschild_check()`, `first_law_kn_check()` | ✓ |
| 3 | Area theorem (Hawking 1971) | `area_increase_theorem()`, `area_theorem_merger_check()` | ✓ |
| 4 | Zeroth law (κ uniform) | `zeroth_law_check()` | ✓ |
| 5 | Third law (κ→0 unreachable) | `surface_gravity_kn()` | ✓ |
| 6 | Hawking area bound | `hawking_area_bound_energy()` | ✓ |
| 7 | Cosmic censorship | `cosmic_censorship_check()` (Lean + C) | ✓ |
| 8 | Area theorem formalization | `areaTheoremNonDecreasing` (Lean) | ✓ |

## L5: Computational Methods (Complete)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | RK4 geodesic integration | `geodesic_rk4_step()`, `geodesic_integrate()` | ✓ |
| 2 | Newton-Raphson orbit finding | `schwarzschild_find_circular_orbit()` | ✓ |
| 3 | Finite-difference Christoffel symbols | `christoffel_symbol_fd()`, `christoffel_fd_local()` | ✓ |
| 4 | QNM frequency computation | `schwarzschild_qnm_fundamental()`, `kerr_qnm_fundamental()` | ✓ |
| 5 | Matched filter SNR | `matched_filter_snr()` | ✓ |
| 6 | Numerical proper distance integration | `schwarzschild_proper_distance()` | ✓ |

## L6: Canonical Systems (Complete)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | Schwarzschild BH (full solution) | `SchwarzschildBH` + all metric functions | ✓ |
| 2 | Kerr BH (full solution) | `KerrBH` + all metric functions | ✓ |
| 3 | Reissner-Nordström BH | `ReissnerNordstromBH` | ✓ |
| 4 | Kerr-Newman BH | `KerrNewmanBH` | ✓ |
| 5 | Schwarzschild QNM spectrum | `schwarzschild_qnm_fundamental()`, overtones | ✓ |
| 6 | Kerr QNM spectrum | `kerr_qnm_fundamental()` | ✓ |
| 7 | Binary BH merger remnant | `final_mass_after_merger()`, `final_spin_after_merger()` | ✓ |

## L7: Applications (Complete — 7 applications)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | GW150914 parameter estimation | `chirp_mass()`, `merger_peak_frequency()` | ✓ |
| 2 | GW150914 area theorem test | `area_theorem_merger_check()`, test verified | ✓ |
| 3 | Ringdown waveform generation | `ringdown_waveform()`, `ringdown_multimode()` | ✓ |
| 4 | BH shadow (EHT) | `kerr_shadow_radius()`, `schwarzschild_critical_impact()` | ✓ |
| 5 | Gravitational wave luminosity | `gw_luminosity_binary()`, `radiated_energy_fraction()` | ✓ |
| 6 | Standard siren cosmology | `luminosity_distance_from_gw()` | ✓ |
| 7 | BH spin evolution (Bardeen) | `bardeen_spin_evolution()` | ✓ |

## L8: Advanced Topics (Complete — 7 topics)
| # | Knowledge Item | Code Artifact | Status |
|---|---------------|---------------|--------|
| 1 | Holographic principle | `holographic_bound_bits()`, `spherical_holographic_bound()` | ✓ |
| 2 | Information paradox / Page curve | `page_time()`, `page_curve_entropy()`, `page_information()` | ✓ |
| 3 | ER = EPR conjecture | `er_bridge_throat_radius()`, `er_bridge_length()` | ✓ |
| 4 | Firewall paradox (AMPS 2012) | `firewall_required()` | ✓ |
| 5 | No-hair theorem multipole test | `kerr_multipoles()`, `no_hair_test_quadrupole()` | ✓ |
| 6 | AdS/CFT and BH thermodynamics | `ads_schwarzschild_temperature()`, `hawking_page_temperature()` | ✓ |
| 7 | Bekenstein bound | `bekenstein_bound()` | ✓ |

## L9: Research Frontiers (Partial — 4 topics documented)
| # | Knowledge Item | Code / Doc Artifact | Status |
|---|---------------|---------------------|--------|
| 1 | Quantum gravity Planck scale | `planck_mass()`, `planck_length()`, `planck_temperature()` | ✓ |
| 2 | BH remnants as dark matter | `remnant_mass()`, `semiclassical_mass_limit()` | ✓ |
| 3 | PBH constraints | `critical_mass_cmb()`, `evaporation_time()` | ✓ |
| 4 | Teukolsky-Starobinsky identities | `teukolsky_starobinsky_check()` | ✓ |
