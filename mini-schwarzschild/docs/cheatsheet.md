# Schwarzschild 解

> MIT 8.962 / Wald Ch.6 / Carroll Ch.5

```python

G_N = 6.67430e-11; C_LIGHT = 299792458.0

def schwarzschild_metric(r, M, G=G_N, c=C_LIGHT):
    """ds^2 = -(1-R_s/r)*c^2*dt^2 + dr^2/(1-R_s/r) + r^2*dOmega^2"""
    R_s = 2.0*G*M/c**2
    g_tt = -(1.0 - R_s/r)
    g_rr = 1.0/(1.0 - R_s/r)
    return g_tt, g_rr

def schwarzschild_radius(M, G=G_N, c=C_LIGHT):
    """R_s = 2*G*M/c^2"""
    return 2.0*G*M/c**2

def gravitational_time_dilation(r, M, G=G_N, c=C_LIGHT):
    """dt/dtau = 1/sqrt(1-R_s/r)"""
    R_s = 2.0*G*M/c**2
    return 1.0/np.sqrt(1.0 - R_s/r)

def gravitational_redshift(z, r_emitter, M, G=G_N, c=C_LIGHT):
    """1+z = 1/sqrt(1-R_s/r)"""
    R_s = 2.0*G*M/c**2
    return 1.0/np.sqrt(1.0 - R_s/r_emitter)

def photon_sphere_radius(M, G=G_N, c=C_LIGHT):
    """r_ph = 3*G*M/c^2 = 1.5*R_s"""
    return 3.0*G*M/c**2

def innermost_stable_circular_orbit(M, G=G_N, c=C_LIGHT):
    """r_ISCO = 6*G*M/c^2 = 3*R_s (Schwarzschild)"""
    return 6.0*G*M/c**2

def orbital_period_keplerian(r, M, G=G_N):
    """T = 2*pi*sqrt(r^3/(G*M))"""
    return 2.0*np.pi*np.sqrt(r**3/(G*M))

def perihelion_precession_per_orbit(M, a, e, G=G_N, c=C_LIGHT):
    """Delta_phi = 6*pi*G*M/(a*c^2*(1-e^2)) [rad/orbit]"""
    return 6.0*np.pi*G*M/(a*c**2*(1.0-e**2))

def light_deflection_angle(b, M, G=G_N, c=C_LIGHT):
    """alpha = 4*G*M/(b*c^2) [rad]"""
    return 4.0*G*M/(b*c**2)

def shapiro_time_delay(Delta_t, r1, r2, M, G=G_N, c=C_LIGHT):
    """Delta_t = 2*G*M/c^3*ln(4*r1*r2/b^2)"""
    R_s = 2.0*G*M/c**2
    return R_s/c*np.log(4.0*r1*r2/schwarzschild_radius(M)**2)

```
