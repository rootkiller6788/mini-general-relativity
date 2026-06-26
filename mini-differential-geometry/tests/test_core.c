#include "dg_manifold.h"
#include "dg_tensor.h"
#include "dg_metric.h"
#include "dg_connection.h"
#include "dg_curvature.h"
#include "dg_geodesic.h"
#include "dg_forms.h"
#include "dg_lie.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define EPS 1e-10
static int passed = 0, failed = 0;
#define T(n) do{printf("  TEST: " n " ... ");}while(0)
#define P() do{printf("PASS\n");passed++;}while(0)
#define F(m) do{printf("FAIL: " m "\n");failed++;}while(0)
#define C(c,m) do{if(c)P();else F(m);}while(0)
#define CE(a,b,m) C(fabs((a)-(b))<EPS,m)
#define CT(c,m) C(c,m)

static void t_manifold(void) {
    T("manifold_create");
    Manifold m = manifold_create(4, "Minkowski");
    C(m.dim==4 && m.num_charts==0, "init");
}
static void t_chart(void) {
    T("chart");
    const char n[][DG_MAX_COORD_NAME]={"t","r","q","p"};
    double lo[4]={0,0,0,0}, hi[4]={1,2,3,4};
    Chart c=chart_create(4,n,lo,hi);
    C(c.dim==4, "dim");
    double pt[4]={0.5,1,1.5,2};
    C(chart_contains_point(&c,pt), "in domain");
    pt[0]=99;
    C(!chart_contains_point(&c,pt), "out domain");
}
static void t_tan_cotan(void) {
    T("tangent/cotangent");
    double vc[4]={1,2,3,4}, wc[4]={5,6,7,8};
    TangentVector v=tangent_vector_create(4,vc,0);
    CotangentVector w=cotangent_vector_create(4,wc,0);
    CE(cotangent_act_on_tangent(&w,&v), 70.0, "omega(X)=70");
}
static void t_tensor_alloc(void) {
    T("tensor_alloc");
    Tensor *t=tensor_alloc(1,1,4);
    C(t!=NULL, "alloc");
    int idx[2]={0,0};
    tensor_set(t,idx,3.14);
    CE(tensor_get(t,idx),3.14,"set/get");
    tensor_free(t);
}
static void t_tensor_add(void) {
    T("tensor_add");
    Tensor *a=tensor_alloc(0,2,4),*b=tensor_alloc(0,2,4);
    int i1[2]={0,1},i2[2]={1,0};
    tensor_set(a,i1,2);tensor_set(a,i2,2);
    tensor_set(b,i1,3);tensor_set(b,i2,3);
    Tensor *c=tensor_add(a,b);
    C(c!=NULL,"add"); CE(tensor_get(c,i1),5.0,"2+3=5");
    tensor_free(a);tensor_free(b);tensor_free(c);
}
static void t_tensor_product(void) {
    T("tensor_product");
    Tensor *v=tensor_alloc(1,0,4),*w=tensor_alloc(0,1,4);
    int i[1]={2},j[1]={3},k[2]={2,3};
    tensor_set(v,i,5);tensor_set(w,j,7);
    Tensor *p=tensor_product(v,w);
    C(p&&p->rank_r==1&&p->rank_s==1,"type(1,1)");
    CE(tensor_get(p,k),35.0,"5*7=35");
    tensor_free(v);tensor_free(w);tensor_free(p);
}
static void t_tensor_contract(void) {
    T("tensor_contract");
    Tensor *d=tensor_alloc(1,1,4);
    for(int mu=0;mu<4;mu++){int i[2]={mu,mu};tensor_set(d,i,1.0);}
    Tensor *tr=tensor_contract(d,0,0);
    C(tr->rank_r==0&&tr->rank_s==0,"scalar");
    CE(tr->components[0],4.0,"trace=4");
    tensor_free(d);tensor_free(tr);
}
static void t_tensor_sym(void) {
    T("tensor_symmetrize");
    double d[16]={0};d[1*4+0]=3;d[0*4+1]=7;
    Tensor *t=tensor_create(0,2,4,d,0);
    C(!tensor_is_symmetric_02(t),"asymmetric");
    Tensor *s=tensor_symmetrize_02(t);
    C(tensor_is_symmetric_02(s),"symmetrized");
    int i1[2]={0,1},i2[2]={1,0};
    CE(tensor_get(s,i1),5.0,"(3+7)/2");
    tensor_free(t);tensor_free(s);
}
static void t_metric(void) {
    T("metric_minkowski");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    CE(metric_get(g,0,0),-1.0,"g00=-1");
    CE(metric_get(g,1,1),1.0,"g11=1");
    CE(metric_determinant(g),-1.0,"det=-1");
    C(metric_verify_signature(g),"signature");
    double dx[4]={2,1,0,0};
    CE(metric_line_element(g,dx),-3.0,"ds2=-3");
    C(metric_classify_dx(g,dx)==-1,"timelike");
    CE(metric_proper_time_inc(g,dx),sqrt(3.0),"dtau=sqrt(3)");
    metric_free(g);
}
static void t_metric_inv(void) {
    T("metric_inverse");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    int ok=1;
    for(int mu=0;mu<4;mu++)for(int rho=0;rho<4;rho++){
        double s=0;
        for(int nu=0;nu<4;nu++)s+=metric_get(g,mu,nu)*metric_inv_get(g,nu,rho);
        if(fabs(s-(mu==rho?1.0:0.0))>EPS)ok=0;
    }
    C(ok,"g*g_inv=I");
    metric_free(g);
}
static void t_conn(void) {
    T("connection_minkowski");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0};
    Connection *c=connection_from_metric(g,dg);
    double mx=0;
    for(int r=0;r<4;r++)for(int m=0;m<4;m++)for(int n=0;n<4;n++)
        if(fabs(connection_get(c,r,m,n))>mx)mx=fabs(connection_get(c,r,m,n));
    C(mx<EPS,"Christoffels=0");
    connection_free(c);metric_free(g);
}
static void t_covd(void) {
    T("covariant_derivative_scalar");
    double dp[4]={1,2,3,4};
    CE(covariant_derivative_scalar(dp,4,0),1.0,"nabla_0=1");
    CE(covariant_derivative_scalar(dp,4,3),4.0,"nabla_3=4");
}
static void t_torsion(void) {
    T("torsion");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0};
    Connection *c=connection_from_metric(g,dg);
    CE(connection_compute_torsion(c),0.0,"torsion-free");
    connection_free(c);metric_free(g);
}
static void t_riemann(void) {
    T("riemann_minkowski");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0},dG[256]={0};
    Connection *c=connection_from_metric(g,dg);
    RiemannTensor *R=riemann_compute(c,dG);
    double mx=0;
    for(int r=0;r<4;r++)for(int s=0;s<4;s++)for(int m=0;m<4;m++)for(int n=0;n<4;n++)
        if(fabs(riemann_get(R,r,s,m,n))>mx)mx=fabs(riemann_get(R,r,s,m,n));
    C(mx<EPS,"Riemann=0");
    riemann_free(R);connection_free(c);metric_free(g);
}
static void t_einstein(void) {
    T("einstein_tensor");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0},dG[256]={0};
    Connection *c=connection_from_metric(g,dg);
    RiemannTensor *R=riemann_compute(c,dG);
    RicciTensor *Ric=ricci_compute(R);
    double Rs=ricci_scalar_compute(Ric,g);
    double *G=einstein_tensor_compute(Ric,g,Rs);
    double mx=0;
    for(int i=0;i<16;i++)if(fabs(G[i])>mx)mx=fabs(G[i]);
    C(mx<EPS,"G_munu=0");
    free(G);ricci_free(Ric);riemann_free(R);connection_free(c);metric_free(g);
}
static void t_bianchi(void) {
    T("bianchi_identity");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0},dG[256]={0};
    Connection *c=connection_from_metric(g,dg);
    RiemannTensor *R=riemann_compute(c,dG);
    CE(bianchi_identity_first_verify(R),0.0,"Bianchi=0");
    riemann_free(R);connection_free(c);metric_free(g);
}
static void t_geo_rhs(void) {
    T("geodesic_rhs");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0};
    Connection *c=connection_from_metric(g,dg);
    double y[8]={0,0,0,0,1,0,0,0},dy[8]={0};
    geodesic_rhs(0,y,dy,(void*)c);
    CE(dy[0],1.0,"dt/dtau=1"); CE(dy[4],0.0,"a^t=0");
    connection_free(c);metric_free(g);
}
static void t_geo_step(void) {
    T("geodesic_step");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    double dg[64]={0};
    Connection *c=connection_from_metric(g,dg);
    GeodesicPoint p; p.dim=4;p.chart_index=0;p.tau=0;
    p.pos[0]=0;p.pos[1]=0;p.pos[2]=0;p.pos[3]=0;
    p.vel[0]=1;p.vel[1]=0;p.vel[2]=0;p.vel[3]=0;
    C(geodesic_step_rk4(&p,c,1.0)==0,"step ok");
    CE(p.pos[0],1.0,"t inc by 1");
    connection_free(c);metric_free(g);
}
static void t_vel_norm(void) {
    T("velocity_norm");
    double e[16]={-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    Metric *g=metric_create(4,SIG_LORENTZIAN_MPPP,e,0);
    GeodesicPoint p; p.dim=4;p.chart_index=0;
    p.vel[0]=1;p.vel[1]=0;p.vel[2]=0;p.vel[3]=0;
    CE(geodesic_velocity_norm(&p,g),-1.0,"norm=-1");
    metric_free(g);
}
static void t_schw(void) {
    T("schwarzschild_init");
    GeodesicPoint p;
    geodesic_schwarzschild_init(10.0,5.0,1.0,&p);
    C(p.dim==4,"4D"); CE(p.pos[1],10.0,"r0=10"); CE(p.vel[3],0.05,"uphi");
}
static void t_lie_scalar(void) {
    T("lie_derivative_scalar");
    double X[4]={1,0,0,0},df[4]={3,2,1,0};
    CE(lie_derivative_scalar(X,df,4),3.0,"L_Xf=3");
}
static void t_killing(void) {
    T("killing_schwarzschild");
    int nk;
    double *k=killing_schwarzschild(4,&nk);
    C(k&&nk==4,"4 Killing"); CE(k[0],1.0,"xi_t^0=1");
    free(k);
}
static void t_maxk(void) {
    T("max_killing");
    C(max_killing_vectors(4)==10,"n=4->10");
    C(max_killing_vectors(3)==6,"n=3->6");
}
static void t_form_create(void) {
    T("form_create");
    DifferentialForm *f=form_alloc(1,4,0);
    C(f&&f->degree==1,"1-form");
    form_free(f);
}

int main(void) {
    printf("=== mini-diff-geom Core Tests ===\n\n");
    printf("--- L1: Manifolds ---\n");
    t_manifold();t_chart();t_tan_cotan();
    printf("\n--- L1-L3: Tensors ---\n");
    t_tensor_alloc();t_tensor_add();t_tensor_product();t_tensor_contract();t_tensor_sym();
    printf("\n--- L1-L4: Metric ---\n");
    t_metric();t_metric_inv();
    printf("\n--- L1-L4: Connection ---\n");
    t_conn();t_covd();t_torsion();
    printf("\n--- L1-L4: Curvature ---\n");
    t_riemann();t_einstein();t_bianchi();
    printf("\n--- L5-L6: Geodesics ---\n");
    t_geo_rhs();t_geo_step();t_vel_norm();t_schw();
    printf("\n--- L1-L2: Lie ---\n");
    t_lie_scalar();t_killing();t_maxk();
    printf("\n--- L2: Forms ---\n");
    t_form_create();
    printf("\n=== %d passed, %d failed ===\n",passed,failed);
    return failed>0?1:0;
}