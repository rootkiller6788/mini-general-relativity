# Course Tree — mini-cosmological-model

## Prerequisites

```
mini-cosmological-model
├── Mathematics
│   ├── Calculus (ODE, PDE, integration)
│   ├── Differential Geometry (metric, curvature)
│   ├── Probability & Statistics
│   └── Numerical Methods (RK4, Gauss-Legendre)
├── Physics
│   ├── General Relativity (Einstein eqs → Friedmann eqs)
│   ├── Thermodynamics (CMB, Boltzmann eq)
│   ├── Quantum Field Theory (inflation, vacuum energy)
│   └── Fluid Dynamics (perfect fluid EM tensor)
└── Astronomy
    ├── Stellar physics (SNIa standard candles)
    ├── Galaxy surveys (power spectrum)
    ├── CMB observations (Planck, WMAP)
    └── Weak lensing (S8 parameter)
```

## Internal Dependencies

```
cosmo_model.h/c     (types, constants)
    ↓
friedmann.h/c       (E(z), H(z), ODE, solutions)
    ↓        ↓
flrw.h/c    distances.h/c     (metric, horizons)
    ↓
perturbations.h/c  (growth, PS mass fn, BAO)
    ↓
power_spectrum.h/c (P(k), transfer fn, halofit)
    ↓
transfer.h         (CMB, SW, ISW, lensing)
```

## L9 Research Frontiers (Documented)
- **Hubble tension**: 5σ discrepancy H0=67.4 vs 73
- **S8 tension**: 2-3σ σ₈√(Ωₘ/0.3) discrepancy
- **Early dark energy**: Proposed Hubble tension solution
- **Primordial GW**: r < 0.036 (Planck/BICEP)
- **Modified gravity**: f(R), DGP, Horndeski alternatives
