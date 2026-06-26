# 黑洞物理

> MIT 8.962 / Wald Ch.8 / Hawking & Ellis

```python

G_N = 6.67430e-11; C_LIGHT = 299792458.0; HBAR = 1.054571817e-34
K_B = 1.380649e-23

def hawking_temperature(M, hbar=HBAR, c=C_LIGHT, G=G_N, k_B=K_B):
    """T_H = hbar*c^3/(8*pi*G*M*k_B)"""
    return hbar*c**3/(8.0*np.pi*G*M*k_B)

def bekenstein_hawking_entropy(A, k_B=K_B, l_P=1.616255e-35):
    """S_BH = k_B*A/(4*l_P^2)"""
    return k_B*A/(4.0*l_P**2)

def black_hole_evaporation_time(M, G=G_N, hbar=HBAR, c=C_LIGHT):
    """tau ~ 5120*pi*G^2*M^3/(hbar*c^4)"""
    return 5120.0*np.pi*G**2*M**3/(hbar*c**4)

def no_hair_theorem():
    """BH characterized by M, J, Q only."""
    return ['M', 'J', 'Q']

def area_increase_theorem(A_initial, A_final):
    """delta_A >= 0 for any classical process"""
    return A_final >= A_initial

def black_hole_thermodynamics_laws(T_H, dM, Omega_H, dJ, Phi_H, dQ):
    """dM = T_H*dS + Omega_H*dJ + Phi_H*dQ (First Law)"""
    return None

def information_loss_paradox():
    """Pure -> mixed state evolution in BH evaporation."""
    return "Unresolved"

def holographic_principle_bound(A, l_P=1.616255e-35):
    """N_bits <= A/(4*l_P^2*ln(2))"""
    return A/(4.0*l_P**2*np.log(2.0))

```
