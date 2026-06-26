# 宇宙学模型

> MIT 8.962 / Wald Ch.5 / Carroll Ch.8

```python

C_LIGHT = 299792458.0; G_N = 6.67430e-11

def flrw_metric(a_t, k, r):
    """ds^2 = -c^2*dt^2 + a(t)^2*(dr^2/(1-k*r^2) + r^2*dOmega^2)"""
    return a_t, k

def friedmann_equation_first(a_dot, a, rho, k, G=G_N, c=C_LIGHT):
    """(a_dot/a)^2 = 8*pi*G*rho/3 - k*c^2/a^2 + Lambda*c^2/3"""
    return a_dot**2/a**2 - 8.0*np.pi*G*rho/3.0

def friedmann_equation_second(a_ddot, a, rho, P, G=G_N, c=C_LIGHT):
    """a_ddot/a = -4*pi*G*(rho+3*P/c^2)/3 + Lambda*c^2/3"""
    return a_ddot/a + 4.0*np.pi*G*(rho+3.0*P/c**2)/3.0

def continuity_equation(rho_dot, rho, P, a_dot, a):
    """rho_dot + 3*a_dot/a*(rho+P/c^2) = 0"""
    return rho_dot + 3.0*a_dot/a*(rho+P/C_LIGHT**2)

def scale_factor_evolution(t, Omega_m, Omega_L, H0=70.0):
    """a(t) ~ (t/t_0)^(2/3) for matter-dominated"""
    return t**(2.0/3.0)

def particle_horizon(a_t, c=C_LIGHT):
    """d_H = c*a*integral_0^t(dt/a)"""
    return c*a_t

def event_horizon(a_t, c=C_LIGHT):
    """d_E = c*a*integral_t^inf(dt/a)"""
    return c*a_t

def curvature_density_parameter(Omega_k, H0, a):
    """Omega_k = -k*c^2/(a^2*H0^2)"""
    return Omega_k

```
