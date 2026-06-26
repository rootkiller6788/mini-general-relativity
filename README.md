# Mini General Relativity

A collection of **from-scratch, zero-dependency C implementations** of core topics in Einstein's general relativity — from differential geometry foundations to black hole thermodynamics, gravitational waves, and numerical relativity. Each sub-module maps to MIT 8.962 (General Relativity) and other top-tier physics curricula, translating textbook equations into runnable, readable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|-----------|--------|-------------|
| [mini-differential-geometry](mini-differential-geometry/) | Affine connection, covariant derivative, Riemann/Ricci/Einstein/Weyl tensors, differential forms, exterior calculus, Hodge duality, Lie derivatives, geodesics | MIT 8.962, Cambridge Part III, ETH 402-0891 |
| [mini-schwarzschild](mini-schwarzschild/) | Schwarzschild metric, event horizon, geodesic motion (ISCO, photon sphere), gravitational time dilation, redshift, light deflection, perihelion precession, Shapiro delay | MIT 8.962, Wald Ch.6, Carroll Ch.5 |
| [mini-kerr-metric](mini-kerr-metric/) | Kerr metric in Boyer-Lindquist & Kerr-Schild coordinates, horizons, ergosphere, ring singularity, Penrose process, geodesics, superradiance | MIT 8.962, Wald Ch.7, Kerr (1963) |
| [mini-einstein-equations](mini-einstein-equations/) | Einstein field equations, Einstein-Hilbert action, stress-energy tensor, curvature, geodesics, coordinate systems, FLRW cosmology | MIT 8.962, Weinberg |
| [mini-gravitational-waves](mini-gravitational-waves/) | Linearized gravity, TT gauge, quadrupole formula, compact binary inspiral/merger, matched filtering, SNR, pulsar timing, gravitational wave memory | MIT 8.962, Maggiore (2008) |
| [mini-black-holes](mini-black-holes/) | BH metrics (Schwarzschild/Kerr/RN/Kerr-Newman), thermodynamics, Hawking radiation, Bekenstein-Hawking entropy, four laws of BH mechanics, information paradox, Penrose process, superradiance | MIT 8.962, Hawking & Ellis (1973) |
| [mini-cosmological-model](mini-cosmological-model/) | FLRW metric, Friedmann equations, scale factor evolution, cosmological distances, power spectrum, cosmological perturbations | MIT 8.962, Dodelson, Carroll Ch.8 |
| [mini-numerical-relativity](mini-numerical-relativity/) | ADM 3+1 decomposition, BSSN formulation, 3D Cartesian grids, extrinsic curvature, horizon finding, initial data | MIT 8.962, Baumgarte & Shapiro (2010), Gourgoulhon (2012) |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes
- **Physics-first design** — tensor operations, geodesic integration, and metric calculations directly mirror the mathematical formalism of GR textbooks (Wald, Carroll, MTW)

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-differential-geometry
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-general-relativity/
├── mini-differential-geometry/  # Affine connection, curvature tensors, differential forms, Lie derivatives
├── mini-schwarzschild/          # Schwarzschild metric, geodesics, time dilation, light deflection
├── mini-kerr-metric/            # Kerr metric, horizons, ergosphere, Penrose process, superradiance
├── mini-einstein-equations/     # Einstein field equations, Einstein-Hilbert action, stress-energy
├── mini-gravitational-waves/    # Linearized gravity, TT gauge, binary inspiral, matched filtering
├── mini-black-holes/            # BH thermodynamics, Hawking radiation, entropy, information paradox
├── mini-cosmological-model/     # FLRW metric, Friedmann equations, scale factor, perturbations
└── mini-numerical-relativity/   # ADM decomposition, BSSN formulation, 3D grids, horizon finding

```

## License

MIT
