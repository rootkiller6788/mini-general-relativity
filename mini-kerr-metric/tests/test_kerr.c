#include "kerr_metric.h"
#include "kerr_horizons.h"
#include "kerr_geodesics.h"
#include "kerr_physics.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int run=0, passed=0;
#define T(n) do{run++;printf("  %s ... ",n);}while(0)
#define P() do{printf("PASS\n");passed++;}while(0)
#define F(m) do{printf("FAIL: %s\n",m);return;}while(0)
#define A(c,m) do{if(!(c)){F(m);}}while(0)
#define N(a,b,t,m) do{if(fabs((a)-(b))>(t)){F(m);}}while(0)

static void t0(void){T("Sigma");N(kerr_sigma(2,M_PI/2,1),4,1e-10,"");P();}
static void t1(void){T("Delta Schw");N(kerr_delta(4,2,0),0,1e-10,"");P();}
static void t2(void){T("Delta Ex");N(kerr_delta(1,1,1),0,1e-10,"");P();}
static void t3(void){
  T("Metric BL Schw");
  KerrBlackHole bh={2,0,1}; KerrBLPoint pt={0,10,M_PI/2,0}; KerrMetricBL m;
  kerr_metric_bl(&bh,&pt,&m);
  N(m.g_tt,-0.6,1e-10,"g_tt"); N(m.g_rr,10.0/6.0,1e-10,"g_rr");
  N(m.g_phph,100,1e-10,"g_phph"); N(m.g_tph,0,1e-10,"g_tphi"); P();}
static void t4(void){
  T("Metric sym"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,5,1,0.5};
  KerrMetricBL m; kerr_metric_bl(&bh,&pt,&m);
  N(m.g_tph,m.g_pht,1e-15,"sym"); P();}
static void t5(void){
  T("Determinant"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,5,1,0.5};
  A(kerr_metric_determinant_bl(&bh,&pt)<0,"det<0"); P();}
static void t6(void){
  T("Frame drag"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,5,M_PI/2,0};
  A(kerr_frame_dragging_omega(&bh,&pt)>0,"omega>0"); P();}
static void t7(void){
  T("Lapse>0"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,10,M_PI/2,0};
  A(kerr_lapse_bl(&bh,&pt)>0,"lapse>0"); P();}
static void t8(void){
  T("Christoffel sym"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,5,1,0.5};
  double ch[64]; kerr_christoffel_bl(&bh,&pt,ch);
  double mx=0;
  for(int mu=0;mu<4;mu++)for(int nu=0;nu<4;nu++)for(int rho=0;rho<4;rho++){
    double d=fabs(ch[16*mu+4*nu+rho]-ch[16*mu+4*rho+nu]);
    if(d>mx)mx=d;}
  N(mx,0,1e-8,"Christoffel asym"); P();}
static void t9(void){
  T("Ricci vanish"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,5,1,0.5};
  double ricci[16]; kerr_ricci_bl(&bh,&pt,ricci);
  double mx=0; for(int i=0;i<16;i++){if(fabs(ricci[i])>mx)mx=fabs(ricci[i]);}
  A(mx<1,"Ricci small"); P();}
static void t10(void){
  T("g*g^-1=I"); KerrBlackHole bh={1,0.5,1}; KerrBLPoint pt={0,5,1,0.5};
  KerrMetricBL m; KerrMetricInvBL inv;
  kerr_metric_bl(&bh,&pt,&m); kerr_inverse_metric_bl(&bh,&pt,&inv);
  double g[16]={m.g_tt,m.g_tr,m.g_tth,m.g_tph,m.g_rt,m.g_rr,m.g_rth,m.g_rph,
    m.g_tht,m.g_thr,m.g_thth,m.g_thph,m.g_pht,m.g_phr,m.g_phth,m.g_phph};
  double gi[16]={inv.inv_tt,inv.inv_tr,inv.inv_tth,inv.inv_tph,
    inv.inv_rt,inv.inv_rr,inv.inv_rth,inv.inv_rph,
    inv.inv_tht,inv.inv_thr,inv.inv_thth,inv.inv_thph,
    inv.inv_pht,inv.inv_phr,inv.inv_phth,inv.inv_phph};
  double mx=0;
  for(int mu=0;mu<4;mu++)for(int rho=0;rho<4;rho++){
    double s=0;for(int nu=0;nu<4;nu++)s+=gi[mu*4+nu]*g[nu*4+rho];
    double dev=fabs(s-(mu==rho?1:0));if(dev>mx)mx=dev;}
  N(mx,0,1e-10,"identity"); P();}
static void t11(void){
  T("Horizon Schw"); KerrBlackHole bh={2,0,1}; double rp,rm;
  kerr_horizon_radii(&bh,&rp,&rm);
  N(rp,4,1e-10,"r_plus"); N(rm,0,1e-10,"r_minus"); P();}
static void t12(void){
  T("Horizon ex"); KerrBlackHole bh={1,1,1}; double rp,rm;
  kerr_horizon_radii(&bh,&rp,&rm);
  N(rp,1,1e-10,"r_plus=1"); N(rm,1,1e-10,"r_minus=1"); P();}
static void t13(void){
  T("Surface grav Schw=1/4M"); KerrBlackHole bh={1,0,1};
  N(kerr_surface_gravity(&bh),0.25,1e-10,"kappa=0.25"); P();}
static void t14(void){
  T("Kappa extremal=0"); KerrBlackHole bh={1,1,1};
  N(kerr_surface_gravity(&bh),0,1e-10,"kappa=0"); P();}
static void t15(void){
  T("Area Schw=16pi M^2"); KerrBlackHole bh={1,0,1};
  N(kerr_horizon_area(&bh),16*M_PI,1e-10,"A=16pi"); P();}
static void t16(void){
  T("Entropy bound"); KerrBlackHole bh={1,0.8,1};
  A(kerr_entropy_bound_ratio(&bh)<=1,"S_Kerr<=S_Schw"); P();}
static void t17(void){
  T("Smarr formula"); KerrBlackHole bh={1,0.5,1};
  N(kerr_smarr_formula_residual(&bh),0,1e-10,"Smarr"); P();}
static void t18(void){
  T("ISCO Schw=6M"); KerrBlackHole bh={1,0,1}; KerrISCOData isco;
  kerr_isco(&bh,1,&isco); N(isco.r_isco,6,1e-10,"r_isco=6"); P();}
static void t19(void){
  T("ISCO prograde ex=M"); KerrBlackHole bh={1,1,1}; KerrISCOData isco;
  kerr_isco(&bh,1,&isco); N(isco.r_isco,1,1e-10,"r_isco=1"); P();}
static void t20(void){
  T("ISCO retrograde ex=9M"); KerrBlackHole bh={1,1,1}; KerrISCOData isco;
  kerr_isco(&bh,0,&isco); N(isco.r_isco,9,1e-10,"r_isco=9"); P();}
static void t21(void){
  T("Photon orbit=3M"); KerrBlackHole bh={1,0,1}; KerrPhotonOrbit orb;
  kerr_photon_orbit(&bh,1,&orb); N(orb.r_photon,3,1e-10,"r_ph=3M"); P();}
static void t22(void){
  T("Penrose max eff"); KerrBlackHole bh={1,1,1};
  N(kerr_penrose_max_efficiency(&bh),1-1/sqrt(2),1e-10,"eta~0.2929"); P();}
static void t23(void){
  T("Penrose Schw=0"); KerrBlackHole bh={1,0,1};
  N(kerr_penrose_max_efficiency(&bh),0,1e-10,"eta=0"); P();}
static void t24(void){
  T("Superradiance ON"); KerrBlackHole bh={1,0.9,1}; SuperradianceCondition c;
  kerr_superradiance_condition(&bh,0.1,1,&c);
  A(c.is_superradiant,"omega=0.1 < Omega_H"); P();}
static void t25(void){
  T("Superradiance OFF"); KerrBlackHole bh={1,0.5,1}; SuperradianceCondition c;
  kerr_superradiance_condition(&bh,10,1,&c);
  A(!c.is_superradiant,"omega=10 >> Omega_H"); P();}
static void t26(void){
  T("BZ power>0"); KerrBlackHole bh={1,0.9,1}; BlandfordZnajek bz;
  kerr_blandford_znajek(&bh,1,&bz); A(bz.power>0,"BZ>0"); P();}
static void t27(void){
  T("Geodesic step"); KerrBlackHole bh={1,0.5,1};
  KerrGeodesicConstants gc={1,3,0,-1,1,0.5};
  KerrGeodesicState st={{0,10,M_PI/2,0}};
  st.p_t=-1; st.p_r=-0.1; st.p_theta=0; st.p_phi=3; st.affine=0;
  A(kerr_geodesic_rk4_step(&bh,&st,0.1,&gc)==0,"step ok"); P();}
static void t28(void){
  T("BL-KS roundtrip"); KerrBLPoint bi={0,5,1,2}, bo; KerrKSPoint ks;
  kerr_bl_to_ks(&bi,0.5,&ks);
  A(kerr_ks_to_bl(&ks,0.5,&bo)==0,"KS->BL ok");
  N(bo.r,bi.r,1e-6,"r preserve"); P();}

int main(void){
  printf("\n=== mini-kerr-metric Test Suite ===\n\n");
  t0();t1();t2();t3();t4();t5();t6();t7();t8();t9();
  t10();t11();t12();t13();t14();t15();t16();t17();t18();t19();
  t20();t21();t22();t23();t24();t25();t26();t27();t28();
  printf("\n=== %d/%d tests passed ===\n",passed,run);
  return (passed==run)?0:1;
}
