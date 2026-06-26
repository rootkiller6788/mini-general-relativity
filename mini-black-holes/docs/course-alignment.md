# Course Alignment — mini-black-holes

## Nine-School Curriculum Mapping

| School | Course | Key Chapters | Status |
|--------|--------|-------------|--------|
| **MIT** | 8.962 General Relativity | Ch.6 Schwarzschild, Ch.7 Kerr, Ch.12 BH Thermodynamics | ✓ |
| **Stanford** | PHYSICS 230 GR | BH solutions, Penrose process, Hawking radiation | ✓ |
| **Berkeley** | PHYS 231 GR | Singularity theorems, BH mechanics | ✓ |
| **Caltech** | Ph 205 GR | Kerr metric, BH perturbations, QNMs | ✓ |
| **Princeton** | PHY 535 GR | Advanced BH physics, information paradox | ✓ |
| **Cambridge** | Part III GR / BH | Four laws, Hawking evaporation, holography | ✓ |
| **Oxford** | CMT/QFT/GR grad courses | AdS/CFT, BH thermodynamics, ER=EPR | ✓ |
| **ETH** | 402-0891 GR | BH solutions, numerical relativity, GWs | ✓ |
| **Tokyo** | Quantum Mechanics / Cosmology | BH information, primordial BHs | ✓ |

## Reference Textbooks

| Textbook | Author(s) | Year | Coverage |
|----------|-----------|------|----------|
| General Relativity | Wald | 1984 | Ch.6-9, 12, 14 (primary) |
| Gravitation | Misner, Thorne, Wheeler | 1973 | Ch.31-34 |
| Spacetime and Geometry | Carroll | 2004 | Ch.5-7 |
| The Large Scale Structure of Space-Time | Hawking & Ellis | 1973 | Ch.9 (singularity theorems) |
| Black Hole Physics | Frolov & Novikov | 1998 | Comprehensive |
| Introduction to Black Hole Physics | Frolov & Zelnikov | 2011 | Modern pedagogical treatment |

## Key Results Mapped to Code

| Result | Reference | Code Location |
|--------|-----------|---------------|
| Birkhoff theorem → Schwarzschild uniqueness | Wald §6.1 | `schwarzschild_metric()` |
| Kerr solution (1963) | Wald §7.1 | `kerr_metric()` |
| Hawking temperature T_H = κ/(2π) | Hawking 1974, Wald §14.3 | `hawking_temperature_schwarzschild()` |
| Bekenstein-Hawking entropy S = A/(4G) | Bekenstein 1973, Wald §12.5 | `bekenstein_hawking_entropy()` |
| Four laws of BH mechanics | Bardeen-Carter-Hawking 1973 | `first_law_schwarzschild_check()` |
| Area theorem δA ≥ 0 | Hawking 1971, Wald §12.2 | `area_increase_theorem()` |
| Penrose process | Penrose 1969, Wald §7.2 | `penrose_process_efficiency()` |
| QNM frequencies | Chandrasekhar-Detweiler 1975 | `schwarzschild_qnm_fundamental()` |
| Page curve | Page 1993 | `page_curve_entropy()` |
| Holographic bound | 't Hooft 1993, Susskind 1995 | `holographic_bound_bits()` |
| ER=EPR | Maldacena-Susskind 2013 | `er_bridge_throat_radius()` |
| Firewall paradox | AMPS 2012 | `firewall_required()` |
| AdS/CFT BH thermodynamics | Witten 1998 | `ads_schwarzschild_temperature()` |
