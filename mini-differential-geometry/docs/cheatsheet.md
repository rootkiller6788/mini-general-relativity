# 微分几何

> MIT 8.962 / Wald / Carroll GR

```python

def metric_tensor_components(g_mu_nu):
    """4x4 metric tensor g_{mu,nu}"""
    return np.asarray(g_mu_nu)

def christoffel_symbols(g, dg):
    """Gamma^rho_{mu,nu} = (1/2)*g^{rho,sigma}*(dg_{sigma mu,nu} + dg_{sigma nu,mu} - dg_{mu nu,sigma})"""
    pass

def riemann_tensor(Gamma, dGamma):
    """R^rho_{sigma,mu,nu} = d_mu(Gamma^rho_{nu,sigma}) - d_nu(Gamma^rho_{mu,sigma}) + ..."""
    pass

def ricci_tensor(R_riemann):
    """R_{mu,nu} = R^rho_{mu,rho,nu}"""
    return None

def ricci_scalar(g_inv, R_mu_nu):
    """R = g^{mu,nu}*R_{mu,nu}"""
    return None

def einstein_tensor(R_mu_nu, R, g_mu_nu):
    """G_{mu,nu} = R_{mu,nu} - R*g_{mu,nu}/2"""
    return R_mu_nu - 0.5*R*g_mu_nu

def covariant_derivative_scalar(V, dV, Gamma):
    """nabla_mu(V^nu) = d_mu(V^nu) + Gamma^nu_{mu,rho}*V^rho"""
    return None

def geodesic_equation(x_mu, dx_mu_dtau, Gamma):
    """d^2(x^mu)/d(tau)^2 + Gamma^mu_{alpha,beta}*dx^alpha/dtau*dx^beta/dtau = 0"""
    pass

def parallel_transport(v_mu, dx_dtau, Gamma):
    """dV^mu/dtau + Gamma^mu_{alpha,beta}*V^alpha*dx^beta/dtau = 0"""
    pass

def proper_time(g_mu_nu, dx_dlambda):
    """d(tau)^2 = g_{mu,nu}*dx^mu*dx^nu"""
    return np.sqrt(abs(np.dot(np.dot(dx_dlambda, g_mu_nu), dx_dlambda)))

```
