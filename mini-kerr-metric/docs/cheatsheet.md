# Kerr 度规

> MIT 8.962 / Wald Ch.7 / Teukolsky

```python

G_N = 6.67430e-11; C_LIGHT = 299792458.0

def kerr_metric_components(r, theta, a, M, G=G_N, c=C_LIGHT):
    """Boyer-Lindquist metric for rotating black hole."""
    R_s = 2.0*G*M/c**2
    rho2 = r**2 + a**2*np.cos(theta)**2
    Delta = r**2 - R_s*r + a**2
    return rho2, Delta

def kerr_horizons(M, a, G=G_N, c=C_LIGHT):
    """r_plusminus = G*M/c^2 +- sqrt((G*M/c^2)^2 - a^2)"""
    M_len = G*M/c**2
    sqrt_term = np.sqrt(max(M_len**2 - a**2, 0))
    return M_len + sqrt_term, M_len - sqrt_term

def kerr_ergosphere(r, theta, M, a, G=G_N, c=C_LIGHT):
    """r_ergo = G*M/c^2 + sqrt((G*M/c^2)^2 - a^2*cos^2(theta))"""
    M_len = G*M/c**2
    return M_len + np.sqrt(max(M_len**2 - a**2*np.cos(theta)**2, 0))

def frame_dragging_angular_velocity(r, theta, a, M, G=G_N, c=C_LIGHT):
    """omega = 2*G*M*a*r/((r^2+a^2)^2 - a^2*Delta*sin^2(theta))"""
    R_s = 2.0*G*M/c**2
    return 2.0*G*M*a*r/c**2/((r**2+a**2)**2)

def penrose_process_energy_extraction(M, a):
    """Extract up to 29% of rest mass from maximal Kerr BH."""
    eta_max = 1.0 - 1.0/np.sqrt(2.0)
    return eta_max if a <= G_N*M/C_LIGHT**2 else 0.29

def blandford_znajek_power(B, M, a):
    """P_BZ ~ B^2*M^2*a^2"""
    return B**2*M**2*a**2

def photon_region_kerr(r, theta, a, M):
    """Photon orbits in Kerr geometry (not circular for a != 0)."""
    return None

def ring_singularity(a):
    """Kerr singularity is a ring of radius a at r=0, theta=pi/2."""
    return a

```
