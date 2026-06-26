# Einstein 方程

> MIT 8.962 / Wald Ch.4 / Carroll Ch.4

```python

G_N = 6.67430e-11; C_LIGHT = 299792458.0
KAPPA = 8.0*np.pi*G_N/C_LIGHT**4

def einstein_field_equations(G_mu_nu, T_mu_nu, kappa=KAPPA):
    """G_{mu,nu} = kappa*T_{mu,nu}"""
    return G_mu_nu - kappa*T_mu_nu

def einstein_hilbert_action(R, g_det):
    """S = (c^4/(16*pi*G))*integral(R*sqrt(-g)*d^4x)"""
    return R*np.sqrt(-g_det)

def stress_energy_perfect_fluid(rho, P, u_mu, g_mu_nu):
    """T_{mu,nu} = (rho+P)*u_mu*u_nu + P*g_{mu,nu}"""
    return (rho+P)*np.outer(u_mu, u_mu) + P*g_mu_nu

def stress_energy_electromagnetic(F_mu_nu, g_mu_nu):
    """T_{mu,nu} = F_{mu,alpha}*F_nu^alpha - g_{mu,nu}*F_{alpha,beta}*F^{alpha,beta}/4"""
    return None

def energy_conditions_check(rho, P):
    """Weak: rho>=0, rho+P>=0. Strong: rho+3P>=0. Dominant: rho>=|P|."""
    return {'WEC': rho>=0 and rho+P>=0, 'SEC': rho+3*P>=0, 'DEC': rho>=abs(P)}

def cosmological_constant_term(g_mu_nu, Lambda=1.1e-52):
    """G_{mu,nu} + Lambda*g_{mu,nu} = kappa*T_{mu,nu}"""
    return Lambda*g_mu_nu

def birkhoff_theorem_statement():
    """Spherically symmetric vacuum solution is static and unique (Schwarzschild)."""
    return True

def linearized_gravity(h_mu_nu, box_operator):
    """Box(hbar_{mu,nu}) = -16*pi*G/c^4*T_{mu,nu} in Lorentz gauge"""
    return None

```
