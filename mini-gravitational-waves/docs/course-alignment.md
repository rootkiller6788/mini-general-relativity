# Course Alignment — mini-gravitational-waves

## Curriculum Mapping

### MIT — 8.962 General Relativity
- Ch.4: Linearized Gravity → `gw_core.h` (tensors, TT gauge)
- Ch.4.4: Quadrupole formula → `gw_quadrupole.h`
- Gravitational wave solutions → `gw_waveforms.h`

### Caltech — Ph 205 General Relativity  
- Linearized field equations → `gw_core.c`
- GW generation and detection → `gw_quadrupole.c`, `gw_detection.c`
- Binary inspiral → `gw_binary.c`

### Princeton — PHY 535 General Relativity
- Perturbation theory → `gw_core.h` (TT projection)
- GW energy-momentum → Isaacson tensor in `gw_quadrupole.c`
- GW astronomy → `gw_detection.c`, examples

### Cambridge — Part III General Relativity
- Linearized theory, TT gauge → `gw_core.h`
- Quadrupole formula and applications → `gw_quadrupole.c`
- GW detection principles → `gw_detection.c`

### Oxford — CMT/QFT/GR Graduate Courses
- GW memory effect → `gw_memory.c`
- Stochastic backgrounds → `gw_stochastic.c`
- Continuous wave sources → `gw_pulsar.c`

### Stanford — PHYSICS 230 General Relativity
- Linearized Einstein equations → `gw_core.c` (TT projection)
- GW polarizations → `gw_core.c` (plus/cross tensors)
- PN expansion → `gw_binary.c` (TaylorT2/4/F2)

### Berkeley — PHYS 231 General Relativity
- GW generation → `gw_quadrupole.c`
- Binary coalescence → `gw_binary.c`
- Detection and data analysis → `gw_detection.c`

### ETH — 402-0891 General Relativity
- Linearized theory → all core files
- GW experiments (LIGO/LISA) → noise models in `gw_detection.c`
- Parameter estimation → Fisher matrix in `gw_detection.c`

### Tokyo — Gravitational Wave Physics
- GW from compact binaries → `gw_binary.c`, `gw_waveforms.c`
- Pulsar timing → `gw_pulsar.c`
- GW cosmology → `gw_stochastic.c` (standard sirens in Lean)

## Reference Textbooks

| Textbook | Author | Chapters Covered |
|----------|--------|-----------------|
| General Relativity | Wald (1984) | Ch.4 — Linearized Gravity |
| Gravitational Waves Vol.1 | Maggiore (2008) | Ch.1-7 — Theory |
| Gravitational Waves Vol.2 | Maggiore (2018) | Ch.23 — Stochastic Backgrounds |
| Living Reviews in Relativity | Blanchet (2014) | PN formalism |
