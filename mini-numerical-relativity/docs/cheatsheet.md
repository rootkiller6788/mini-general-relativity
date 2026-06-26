# 数值相对论

> MIT 8.962 / Baumgarte & Shapiro / Alcubierre

```python

def adm_formulation_3plus1(gamma_ij, K_ij, alpha_lapse, beta_shift):
    """ADM 3+1 decomposition: ds^2 = -alpha^2*dt^2 + gamma_ij*(dx^i+beta^i)(dx^j+beta^j)"""
    return gamma_ij, K_ij

def bssn_formulation_variables(phi, gamma_tilde, K, A_tilde, Gamma_tilde):
    """BSSNOK formulation for numerical stability."""
    return phi, gamma_tilde, K, A_tilde, Gamma_tilde

def puncture_initial_data(M1, M2, x1, x2, P1, P2, S1, S2):
    """Bowen-York puncture initial data for binary BHs."""
    return M1+M2

def apparent_horizon_finder(expansion_Theta):
    """Theta = 0 on apparent horizon (marginally outer trapped surface)."""
    return abs(expansions_Theta) < 1e-10

def gravitational_wave_extraction(psi4, r, l=2, m=2):
    """Extract h from Weyl scalar psi4 at large r."""
    return None

def excision_method(r_excise):
    """Avoid singularity by excising interior region."""
    return r_excise

def gauge_conditions_1pluslog(alpha, phi):
    """1+log slicing: d(alpha)/dt = -2*alpha*K"""
    return -2.0*alpha

def gamma_driver_shift(beta_i, Gamma_i, eta=1.0):
    """Gamma-driver shift condition."""
    return None

```
