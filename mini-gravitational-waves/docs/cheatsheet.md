# 引力波

> MIT 8.962 / Maggiore / Wald Ch.4

```python

G_N = 6.67430e-11; C_LIGHT = 299792458.0

def quadrupole_formula_h_ij(I_ddot_ij, r, G=G_N, c=C_LIGHT):
    """h_ij = (2*G/(c^4*r))*d^2(I_ij)/dt^2 (TT gauge)"""
    return 2.0*G/(c**4*r)*I_ddot_ij

def quadrupole_power(I_ddot_ddot_ij, G=G_N, c=C_LIGHT):
    """P = (G/(5*c^5))*<I_ddot_ddot_ij>"""
    return G/(5.0*c**5)*I_ddot_ddot_ij

def two_body_chirp_mass(m1, m2):
    """M_c = (m1*m2)^(3/5)/(m1+m2)^(1/5)"""
    return (m1*m2)**0.6/(m1+m2)**0.2

def gw_strain_from_binary(M_c, D_L, f, G=G_N, c=C_LIGHT):
    """h ~ (G*M_c/c^2/D_L)*(G*M_c/c^3*pi*f)^(2/3)"""
    return G*M_c/c**2/D_L*(G*M_c/c**3*np.pi*f)**(2.0/3.0)

def gw_frequency_evolution(f, M_c, G=G_N, c=C_LIGHT):
    """df/dt = 96/5*pi^(8/3)*(G*M_c/c^3)^(5/3)*f^(11/3)"""
    return 96.0/5.0*np.pi**(8.0/3.0)*(G*M_c/c**3)**(5.0/3.0)*f**(11.0/3.0)

def memory_effect_Delta_h(A, B=1.0):
    """Gravitational wave memory: permanent displacement after GW passage."""
    return A

def continuous_wave_from_pulsar(epsilon, I, f_rot, D):
    """h ~ (4*pi^2*G/(c^4))*epsilon*I*f_rot^2/D"""
    return 4.0*np.pi**2*G_N/C_LIGHT**4*epsilon*I*f_rot**2/D

def stochastic_background_energy_density(Omega_GW):
    """Omega_GW = (1/rho_c)*d(rho_GW)/d(ln f)"""
    return Omega_GW

```
