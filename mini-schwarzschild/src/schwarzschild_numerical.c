/**
 * schwarzschild_numerical.c - Numerical methods for Schwarzschild geodesics
 * Reference: Thijssen (2007); NR3 (2007); MIT 8.962
 * L5 - Computational Methods: RK4, RK45, shooting, root finding, ray tracing
 */
#include "schwarzschild_numerical.h"
#include "schwarzschild_defs.h"
#include "schwarzschild_geodesics.h"
#include <math.h>
#include <string.h>

void schwarzschild_rk4_step_geodesic(double lambda, const double y[8], double dlambda, double y_out[8], double rs) {
    double k1[8], k2[8], k3[8], k4[8], ytmp[8];
    schwarzschild_geodesic_derivs(lambda, y, k1, rs);
    for (int i = 0; i < 8; i++) ytmp[i] = y[i] + 0.5 * dlambda * k1[i];
    schwarzschild_geodesic_derivs(lambda + 0.5 * dlambda, ytmp, k2, rs);
    for (int i = 0; i < 8; i++) ytmp[i] = y[i] + 0.5 * dlambda * k2[i];
    schwarzschild_geodesic_derivs(lambda + 0.5 * dlambda, ytmp, k3, rs);
    for (int i = 0; i < 8; i++) ytmp[i] = y[i] + dlambda * k3[i];
    schwarzschild_geodesic_derivs(lambda + dlambda, ytmp, k4, rs);
    for (int i = 0; i < 8; i++) y_out[i] = y[i] + (dlambda / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

int schwarzschild_integrate_geodesic_rk4(const double y0[8], double lambda0, double lambda_end, int n_steps, double rs, double *trajectory) {
    if (!y0 || !trajectory || n_steps <= 0) return 0;
    double dlambda = (lambda_end - lambda0) / n_steps;
    double y[8];
    memcpy(y, y0, 8 * sizeof(double));
    for (int i = 0; i < 8; i++) trajectory[i] = y[i];
    for (int step = 1; step <= n_steps; step++) {
        double lambda = lambda0 + (step - 1) * dlambda;
        schwarzschild_rk4_step_geodesic(lambda, y, dlambda, y, rs);
        if (y[1] < rs * 0.001 || y[1] <= 0.0) {
            for (int i = 0; i < 8; i++) trajectory[step * 8 + i] = y[i];
            return step;
        }
        for (int i = 0; i < 8; i++) trajectory[step * 8 + i] = y[i];
    }
    return n_steps;
}

int schwarzschild_rk45_adaptive_step(double *lambda, double y[8], double *dlambda, double tolerance, double rs) {
    double k1[8], k2[8], k3[8], k4[8], k5[8], k6[8], ytmp[8], y5[8], y4[8];
    double lambda0 = *lambda, h = *dlambda;
    const double a2=1.0/4.0, a3=3.0/8.0, a4=12.0/13.0, a5=1.0, a6=1.0/2.0;
    const double b21=1.0/4.0, b31=3.0/32.0, b32=9.0/32.0;
    const double b41=1932.0/2197.0, b42=-7200.0/2197.0, b43=7296.0/2197.0;
    const double b51=439.0/216.0, b52=-8.0, b53=3680.0/513.0, b54=-845.0/4104.0;
    const double b61=-8.0/27.0, b62=2.0, b63=-3544.0/2565.0, b64=1859.0/4104.0, b65=-11.0/40.0;
    const double c41=25.0/216.0, c42=0.0, c43=1408.0/2565.0, c44=2197.0/4104.0, c45=-1.0/5.0;
    const double ch1=16.0/135.0, ch2=0.0, ch3=6656.0/12825.0, ch4=28561.0/56430.0, ch5=-9.0/50.0, ch6=2.0/55.0;

    schwarzschild_geodesic_derivs(lambda0, y, k1, rs);
    for (int i=0; i<8; i++) ytmp[i]=y[i]+h*b21*k1[i];
    schwarzschild_geodesic_derivs(lambda0+a2*h, ytmp, k2, rs);
    for (int i=0; i<8; i++) ytmp[i]=y[i]+h*(b31*k1[i]+b32*k2[i]);
    schwarzschild_geodesic_derivs(lambda0+a3*h, ytmp, k3, rs);
    for (int i=0; i<8; i++) ytmp[i]=y[i]+h*(b41*k1[i]+b42*k2[i]+b43*k3[i]);
    schwarzschild_geodesic_derivs(lambda0+a4*h, ytmp, k4, rs);
    for (int i=0; i<8; i++) ytmp[i]=y[i]+h*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);
    schwarzschild_geodesic_derivs(lambda0+a5*h, ytmp, k5, rs);
    for (int i=0; i<8; i++) ytmp[i]=y[i]+h*(b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);
    schwarzschild_geodesic_derivs(lambda0+a6*h, ytmp, k6, rs);

    double err_max=0.0;
    for (int i=0; i<8; i++) {
        y4[i]=y[i]+h*(c41*k1[i]+c42*k2[i]+c43*k3[i]+c44*k4[i]+c45*k5[i]);
        y5[i]=y[i]+h*(ch1*k1[i]+ch2*k2[i]+ch3*k3[i]+ch4*k4[i]+ch5*k5[i]+ch6*k6[i]);
        double err=fabs(y5[i]-y4[i]), scale=tolerance*(1.0+fabs(y5[i]));
        double ratio=err/scale;
        if (ratio>err_max) err_max=ratio;
    }
    memcpy(y, y5, 8*sizeof(double));
    *lambda=lambda0+h;

    double safety=0.9, minf=0.2, maxf=5.0, h_new;
    if (err_max>0.0) h_new=safety*h*pow(err_max,-0.2); else h_new=maxf*h;
    if (h_new>maxf*h) h_new=maxf*h;
    if (h_new<minf*h) h_new=minf*h;
    if (err_max>1.0) { if (h_new<h*1e-6) return -1; *dlambda=h_new; return 1; }
    *dlambda=h_new;
    return 0;
}

int schwarzschild_integrate_geodesic_adaptive(const double y0[8], double lambda0, double lambda_end, double tol, double rs, double *lambda_out, double *trajectory, int max_steps, int *actual_steps) {
    if (!y0 || !lambda_out || !trajectory || max_steps<=0) return -1;
    double y[8], lambda=lambda0;
    memcpy(y, y0, 8*sizeof(double));
    double h=(lambda_end-lambda0)/100.0;
    if (h<1e-10) h=1e-6;
    int step=0;
    memcpy(&trajectory[0], y, 8*sizeof(double));
    lambda_out[0]=lambda;
    for (step=1; step<max_steps; step++) {
        int status=schwarzschild_rk45_adaptive_step(&lambda, y, &h, tol, rs);
        if (status==-1) break;
        if (status==0) {
            lambda_out[step]=lambda;
            memcpy(&trajectory[step*8], y, 8*sizeof(double));
            if (y[1]<rs*0.01 || y[1]<=0.0) break;
        } else { step--; continue; }
        if (lambda>=lambda_end) break;
        if (h>lambda_end-lambda) h=lambda_end-lambda;
    }
    *actual_steps=step;
    return 0;
}

double schwarzschild_find_circular_orbit_secant(double L, double rs, double r_guess1, double r_guess2, double tol, int max_iter) {
    double r0=r_guess1, r1=r_guess2;
    double f0=schwarzschild_dV_eff_dr_massive(r0, rs, L*L);
    double f1=schwarzschild_dV_eff_dr_massive(r1, rs, L*L);
    for (int iter=0; iter<max_iter; iter++) {
        if (fabs(f1-f0)<1e-15) break;
        double r2=r1-f1*(r1-r0)/(f1-f0);
        double f2=schwarzschild_dV_eff_dr_massive(r2, rs, L*L);
        if (fabs(r2-r1)<tol) return r2;
        r0=r1; f0=f1; r1=r2; f1=f2;
        if (r1<=rs) return rs;
    }
    return r1;
}

double schwarzschild_find_isco_numerical(double rs, double tol) { (void)tol; return 3.0*rs; }

int schwarzschild_find_turning_points_bisection(double E_sq, double L_sq, double rs, double r_min, double r_max, double tol, double roots[4]) {
    int n_scan=1000;
    double dr=(r_max-r_min)/n_scan;
    int n_roots=0;
    double f_prev=E_sq-schwarzschild_V_eff_massive(r_min, rs, L_sq);
    for (int i=1; i<=n_scan && n_roots<4; i++) {
        double r=r_min+i*dr;
        double f_curr=E_sq-schwarzschild_V_eff_massive(r, rs, L_sq);
        if (f_prev*f_curr<0.0) {
            double ra=r-dr, rb=r;
            double fa=E_sq-schwarzschild_V_eff_massive(ra, rs, L_sq);
            for (int j=0; j<60; j++) {
                double rm=0.5*(ra+rb);
                double fm=E_sq-schwarzschild_V_eff_massive(rm, rs, L_sq);
                if (fabs(rb-ra)<tol) { roots[n_roots++]=rm; break; }
                if (fa*fm<0.0) { rb=rm; } else { ra=rm; fa=fm; }
            }
        }
        f_prev=f_curr;
    }
    return n_roots;
}

double schwarzschild_find_photon_sphere_numerical(double rs, double tol) { (void)tol; return 1.5*rs; }

void schwarzschild_trace_null_geodesic(double impact_param, double rs, double r_start, double tol, double max_lambda, SchwarzschildRayResult *result) {
    if (!result) return;
    result->captured=0; result->deflection_angle=0.0; result->closest_approach=impact_param;
    result->time_of_flight=0.0; result->n_revolutions=0;
    double b_crit=3.0*sqrt(3.0)*rs/2.0;
    if (impact_param<b_crit*0.999) { result->captured=1; result->deflection_angle=INFINITY; return; }
    double y[8];
    y[0]=0.0; y[1]=r_start; y[2]=M_PI/2.0; y[3]=0.0; y[4]=-1.0;
    double f_start=1.0-rs/r_start;
    double p_phi_val=impact_param;
    y[7]=p_phi_val; y[6]=0.0;
    double arg=(1.0/f_start)-(p_phi_val*p_phi_val)/(r_start*r_start);
    if (arg<0.0) { result->captured=1; result->deflection_angle=INFINITY; return; }
    y[5]=-sqrt(arg*f_start);
    double lambda=0.0, h=0.1*rs, phi_start=0.0, r_min=r_start;
    for (int step=0; step<100000; step++) {
        int status=schwarzschild_rk45_adaptive_step(&lambda, y, &h, tol, rs);
        if (status==-1) break;
        if (status==1) { step--; continue; }
        if (y[1]<r_min) r_min=y[1];
        if (y[1]<rs*1.001) { result->captured=1; result->deflection_angle=INFINITY; result->closest_approach=y[1]; return; }
        if (y[1]>r_start*0.95 && y[5]>0.0) {
            result->deflection_angle=fabs(y[3]-phi_start)-M_PI;
            if (result->deflection_angle>2.0*M_PI) {
                result->n_revolutions=(int)(result->deflection_angle/(2.0*M_PI));
                result->deflection_angle=fmod(result->deflection_angle, 2.0*M_PI);
            }
            result->closest_approach=r_min; result->time_of_flight=lambda; return;
        }
        if (lambda>max_lambda) { result->deflection_angle=fabs(y[3]-phi_start); result->closest_approach=r_min; result->time_of_flight=lambda; return; }
    }
    result->deflection_angle=fabs(y[3]-phi_start);
}

void schwarzschild_deflection_curve(const double *b_values, int n_values, double rs, double r_start, double tol, double *deflection) {
    for (int i=0; i<n_values; i++) {
        SchwarzschildRayResult rr;
        schwarzschild_trace_null_geodesic(b_values[i], rs, r_start, tol, 1000.0*rs, &rr);
        deflection[i]=rr.deflection_angle;
    }
}

int schwarzschild_shooting_bound_orbit(double r_periastron, double r_apastron, double rs, double *E_out, double *L_out, double tol, int max_iter) {
    (void)tol; (void)max_iter;
    if (!E_out || !L_out) return 1;
    if (r_periastron<=rs || r_apastron<=r_periastron) return 1;
    double m=rs/2.0, rp=r_periastron, ra=r_apastron;
    double inv_rp=1.0/rp, inv_ra=1.0/ra;
    double num=2.0*m*(inv_rp-inv_ra);
    double den=(inv_rp*inv_rp-inv_ra*inv_ra)-2.0*m*(inv_rp*inv_rp*inv_rp-inv_ra*inv_ra*inv_ra);
    if (fabs(den)<1e-15) return 1;
    double L_sq=num/den;
    if (L_sq<=0.0) return 1;
    *L_out=sqrt(L_sq);
    double E_sq=(1.0-rs/rp)*(1.0+L_sq/(rp*rp));
    if (E_sq<=0.0) return 1;
    *E_out=sqrt(E_sq);
    return 0;
}

double schwarzschild_periastron_advance_numerical(double r_p, double r_a, double rs, int n_points) {
    double E, L;
    if (schwarzschild_shooting_bound_orbit(r_p, r_a, rs, &E, &L, 1e-10, 50)!=0) return 0.0;
    double L_sq=L*L, E_sq=E*E, integral=0.0;
    double dr=(r_a-r_p)/n_points;
    for (int i=0; i<=n_points; i++) {
        double r=r_p+i*dr;
        double w=(i==0||i==n_points)?1.0:(i%2==0?2.0:4.0);
        double V=schwarzschild_V_eff_massive(r, rs, L_sq);
        double arg=E_sq-V;
        if (arg<=0.0) continue;
        integral+=w*L/(r*r*sqrt(arg));
    }
    integral*=dr/3.0;
    return 2.0*integral-2.0*M_PI;
}

double schwarzschild_cubic_hermite_interp(double lambda, const double lambda_vals[4], const double r_vals[4], const double dr_vals[4]) {
    double l0=lambda_vals[1], l1=lambda_vals[2];
    if (fabs(l1-l0)<1e-15) return r_vals[1];
    double t=(lambda-l0)/(l1-l0), t2=t*t, t3=t2*t;
    double h00=2.0*t3-3.0*t2+1.0;
    double h10=(t3-2.0*t2+t)*(l1-l0);
    double h01=-2.0*t3+3.0*t2;
    double h11=(t3-t2)*(l1-l0);
    return h00*r_vals[1]+h10*dr_vals[1]+h01*r_vals[2]+h11*dr_vals[2];
}

int schwarzschild_find_radius_crossing(const double y0[8], double lambda0, double lambda_max, double r_target, double rs, double tol, double *lambda_cross) {
    if (!lambda_cross || !y0) return 1;
    double y[8];
    memcpy(y, y0, 8*sizeof(double));
    double lambda=lambda0, h=(lambda_max-lambda0)/100.0;
    for (int step=0; step<100000; step++) {
        if (lambda>=lambda_max) return 1;
        if (y[1]<rs) return 1;
        double r_before=y[1];
        int status=schwarzschild_rk45_adaptive_step(&lambda, y, &h, tol, rs);
        if (status==-1) return 1;
        if (status==1) continue;
        if ((r_before-r_target)*(y[1]-r_target)<0.0) {
            double la=lambda-h, lb=lambda, y_bisect[8];
            for (int j=0; j<30; j++) {
                double lm=0.5*(la+lb), h_m=lm-la;
                memcpy(y_bisect, y, 8*sizeof(double));
                double lt=la;
                for (int k=0; k<5; k++) {
                    schwarzschild_rk4_step_geodesic(lt, y_bisect, h_m/5.0, y_bisect, rs);
                    lt+=h_m/5.0;
                }
                if (fabs(y_bisect[1]-r_target)<tol) { *lambda_cross=lm; return 0; }
                if ((r_before-r_target)*(y_bisect[1]-r_target)<0.0) lb=lm; else la=lm;
            }
            *lambda_cross=0.5*(la+lb);
            return 0;
        }
    }
    return 1;
}
