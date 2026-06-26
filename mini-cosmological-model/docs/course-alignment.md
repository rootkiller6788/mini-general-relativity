# Course Alignment — mini-cosmological-model

## MIT — 8.962 GR / 8.942 Cosmology
| Chapter | Topic | Implementation |
|---------|-------|----------------|
| Wald Ch.5 | FLRW metric, Friedmann eqs | flrw.c, friedmann.c |
| Carroll Ch.8 | Cosmological models | cosmo_model.c |
| Dodelson Ch.2 | Background cosmology | friedmann_E_z() |
| Dodelson Ch.7-8 | Structure formation | perturbations.c |
| Dodelson Ch.9-10 | CMB anisotropies | transfer.c |

## Stanford — PHYSICS 230 GR / 450 Adv QM
| Topic | Implementation |
|-------|----------------|
| FLRW geometry, conformal time | flrw_conformal_time() |
| Perturbation theory | perturbations_growth_factor_D() |
| Inflationary spectrum | ps_primordial_dimensionless() |

## Berkeley — PHYS 231 GR / 242 CM
| Topic | Implementation |
|-------|----------------|
| Friedmann solutions | friedmann_scale_*() |
| Distance ladder | distances_luminosity(), distances_distance_modulus() |
| LSS | ps_correlation_function() |

## Caltech — Ph 205 GR / 230 QFT
| Topic | Implementation |
|-------|----------------|
| de Sitter space | friedmann_scale_de_sitter(), distances_event_horizon() |
| Quantum fluctuations | ps_primordial_dimensionless() |
| CMB polarization | transfer_Cl_EE() |

## Princeton — PHY 535 GR / 545 QFT
| Topic | Implementation |
|-------|----------------|
| FLRW dynamics | friedmann_rk4_step() |
| Etherington reciprocity | test_etherington_reciprocity() |
| Sachs-Wolfe effect | transfer_Cl_SW() |

## Cambridge — Part III Cosmology
| Topic | Implementation |
|-------|----------------|
| Redshift-distance relation | distances_comoving_line_of_sight() |
| Growth of perturbations | perturbations_growth_factor_D() |
| BAO physics | perturbations_bao_oscillation() |

## Oxford — CMT/QFT/GR Graduate Courses
| Topic | Implementation |
|-------|----------------|
| Mathematical cosmology | flrw_ricci_scalar() |
| Cosmography | flrw_jerk_parameter(), flrw_snap_parameter() |
| Spherical collapse | perturbations_delta_c() |

## ETH — 402-0891 GR
| Topic | Implementation |
|-------|----------------|
| Friedmann-Lemaitre models | friedmann_H_z() |
| Transfer functions | ps_transfer_BBKS() |
| Nonlinear structure | ps_nonlinear_halofit() |

## 东京大学 — 宇宙論
| Topic | Implementation |
|-------|----------------|
| Hubble tension | examples/age_of_universe.c |
| S8 tension | transfer_S8() |
| Neutrino cosmology | transfer_neutrino_free_streaming() |

## Core Textbooks
| Textbook | Chapters |
|----------|---------|
| Dodelson & Schmidt (2020) | Ch. 2, 4, 7-10, 13 |
| Carroll (2004) | Ch. 8 |
| Wald (1984) | Ch. 5 |
| Weinberg (2008) | Ch. 1-3, 6-8 |
| Peebles (1993) | Ch. 5, 13-14 |
| Mo, van den Bosch, White (2010) | Ch. 4, 7-8 |
