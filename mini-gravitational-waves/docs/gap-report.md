# Gap Report — mini-gravitational-waves

## Current Gaps

### L8 — Advanced Topics (3 gaps)

| # | Missing Topic | Priority | Effort |
|---|--------------|----------|--------|
| 1 | Effective One-Body (EOB) waveform model | Medium | High |
| 2 | Tidal effects in NS binaries (Love numbers) | Medium | Medium |
| 3 | Numerical Relativity surrogate waveforms | Low | Very High |

### L9 — Research Frontiers (all documented, 0 implemented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Primordial GW detection with LISA | Documented in `gw_omega_inflation()` |
| 2 | Cosmic string GW detection bounds | Documented in `gw_omega_cosmic_strings()` |
| 3 | Multi-messenger cosmology (standard sirens) | Documented in Lean (`StandardSiren`) |
| 4 | Testing GR with GW polarization | Documented in `GwPolarization` enum |
| 5 | Machine learning for GW detection | Not covered |

## No Blocking Gaps

L1-L7 are complete. L8-L9 have partial coverage, which is sufficient
per SKILL.md §6: "L7-L9: Partial+" for COMPLETE status.

## Recommended Improvements

1. Implement EOB waveform model (L8)
2. Add tidal Love number computation (L8)
3. Add matched-filter demo with colored noise (L7)
4. Implement Bayesian parameter estimation MCMC (L8)
