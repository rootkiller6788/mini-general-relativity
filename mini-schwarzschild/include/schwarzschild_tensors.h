/**
 * schwarzschild_tensors.h — Curvature tensors in Schwarzschild spacetime
 *
 * Reference: Wald (1984) Ch.3, Ch.6; Carroll (2004) Ch.3, Ch.5
 * MIT 8.962 — General Relativity
 *
 * Covers:
 *   L3 — Mathematical Structures: Christoffel symbols, Riemann/Ricci/Einstein/Weyl tensors
 *   L4 — Fundamental Laws: Vacuum Einstein equations R_{mu nu} = 0
 */

#ifndef SCHWARZSCHILD_TENSORS_H
#define SCHWARZSCHILD_TENSORS_H

#include "schwarzschild_defs.h"

typedef struct {
    double Gamma_ttr;
    double Gamma_trt;
    double Gamma_rtt;
    double Gamma_rrr;
    double Gamma_rthetatheta;
    double Gamma_rphiphi;
    double Gamma_thetartheta;
    double Gamma_thetathetar;
    double Gamma_thetaphphi;
    double Gamma_phirphi;
    double Gamma_phiphir;
    double Gamma_phithetaphi;
    double Gamma_phiphitheta;
} SchwarzschildChristoffel;

void schwarzschild_christoffel(const SchwarzschildPoint *p, double rs,
                                SchwarzschildChristoffel *chris);
double schwarzschild_christoffel_component(int mu, int alpha, int beta,
                                            double r, double theta, double rs);
void schwarzschild_covariant_derivative_vector(
    const SchwarzschildVector *V, const double dV[4][4],
    const SchwarzschildPoint *p, double rs, double cov_deriv[4][4]);
void schwarzschild_covariant_derivative_oneform(
    const SchwarzschildVector *omega, const double domega[4][4],
    const SchwarzschildPoint *p, double rs, double cov_deriv[4][4]);

typedef struct {
    double R_trtr;
    double R_tthetattheta;
    double R_tphitphi;
    double R_rthetartheta;
    double R_rphirphi;
    double R_thetaphithetaphi;
} SchwarzschildRiemann;

void schwarzschild_riemann(const SchwarzschildPoint *p, double rs,
                            SchwarzschildRiemann *riemann);
double schwarzschild_riemann_component(int mu, int nu, int alpha, int beta,
                                        double r, double theta, double rs);
void schwarzschild_ricci(const SchwarzschildPoint *p, double rs,
                          SchwarzschildSymmetric4x4 *ricci);
double schwarzschild_ricci_scalar(const SchwarzschildPoint *p, double rs);
void schwarzschild_einstein(const SchwarzschildPoint *p, double rs,
                             SchwarzschildSymmetric4x4 *einstein);

typedef struct {
    double C_trtr;
    double C_tthetattheta;
    double C_tphitphi;
    double C_rthetartheta;
    double C_rphirphi;
    double C_thetaphithetaphi;
} SchwarzschildWeyl;

void schwarzschild_weyl(const SchwarzschildPoint *p, double rs,
                         SchwarzschildWeyl *weyl);
double schwarzschild_kretschmann(double r, double rs);
double schwarzschild_chern_pontryagin(const SchwarzschildPoint *p, double rs);
double schwarzschild_euler_invariant(double r, double rs);
void schwarzschild_tidal_tensor(const SchwarzschildPoint *p, double rs,
                                 double E[3][3]);
void schwarzschild_geodesic_deviation(const double deviation[3],
                                       const SchwarzschildPoint *p, double rs,
                                       double relative_accel[3]);
double schwarzschild_bel_robinson_super_energy(const SchwarzschildPoint *p, double rs);

#endif /* SCHWARZSCHILD_TENSORS_H */