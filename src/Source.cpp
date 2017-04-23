#include "Source.h"
#include "Impedance.h"
#include "TriangleQuadrature.h"
#include "MathUtils.h"
#include "RWGTree.h"
#include "GreensFunction.h"
#include "Options.h"
#include "Matrix.h"
#include "DenseMatrix.h"
#include "HMatrix.h"
#include "MatrixOperations.h"


DipolePointSource :: DipolePointSource(GreensFunction *g_, Point &p_, Vector3 &j_) : g(g_), p(p_), j(j_) {}

void EHNonSmoothIntegral(HalfRWG &sh, const Vector3 &rt, Real *ISn3h, Real *ISn1, Real *IS1, Vector3 *ILn1, Vector3 *IL1, Vector3 *IL3)
{
  Real tol = 1e-6 * sh.e->length;

  Vector3 p0 = rt - *sh.p0;
  Vector3 p1 = rt - *sh.p1;
  Vector3 p2 = rt - *sh.p2;

  Real R02 = norm(p0);
  Real R0 = sqrt(R02);
  Real R12 = norm(p1);
  Real R1 = sqrt(R12);
  Real R22 = norm(p2);
  Real R2 = sqrt(R22);
  
  Vector3 a0 = 1. / R0 * p0;
  Vector3 a1 = 1. / R1 * p1;
  Vector3 a2 = 1. / R2 * p2;

  Vector3 &n = sh.n;

  Real x = 1.0 + dot(a0,a1+a2) + dot(a1,a2);
  Real y = sqrt(norm(dot(a0,a1*a2)));

  Real h = dot(p0,n);

  *ISn3h = (h<0.?-2.:2.) * atan2(y,x);

  Real num, den;
  Real s010 = sh.normv01i * dot(sh.v01,p0);
  Real s011 = sh.normv01i * dot(sh.v01,p1);
  num = (R1 - s011);
  den = (R0 - s010);
  Real ILn1_01 = fabs(num) < tol || fabs(den) < tol?0.:log(num / den);

  Real s121 = sh.normv12i * dot(sh.v12,p1);
  Real s122 = sh.normv12i * dot(sh.v12,p2);
  num = (R2 - s122);
  den = (R1 - s121);
  Real ILn1_12 = fabs(num) < tol || fabs(den) < tol?0.:log(num / den);

  Real s202 = sh.normv20i * dot(sh.v20,p2);
  Real s200 = sh.normv20i * dot(sh.v20,p0);
  num = (R0 - s200);
  den = (R2 - s202);
  Real ILn1_20 = fabs(num) < tol || fabs(den) < tol?0.:log(num / den);

  Real t01 = dot(sh.m01,p0);
  Real t12 = dot(sh.m12,p1);
  Real t20 = dot(sh.m20,p2);
  
  *ISn1 = -h * *ISn3h - (t01 * ILn1_01 + t12 * ILn1_12 + t20 * ILn1_20);

  Real R02_01 = (R02 - s010 * s010);
  Real IL1_01 = 0.5 * (R02_01 * ILn1_01 - (R1 * s011 - R0 * s010));

  Real R02_12 = (R12 - s121 * s121);
  Real IL1_12 = 0.5 * (R02_12 * ILn1_12 - (R2 * s122 - R1 * s121));

  Real R02_20 = (R22 - s202 * s202);
  Real IL1_20 = 0.5 * (R02_20 * ILn1_20 - (R0 * s200 - R2 * s202));
  
  Real IL3_01 = 0.75 * R02_01 * IL1_01 - 0.25 * (R12 * R1 * s011 - R02 * R0 * s010);
  Real IL3_12 = 0.75 * R02_12 * IL1_12 - 0.25 * (R22 * R2 * s122 - R12 * R1 * s121);
  Real IL3_20 = 0.75 * R02_20 * IL1_20 - 0.25 * (R02 * R0 * s200 - R22 * R2 * s202);

  *IS1 = ONETHIRD * (h * h * *ISn1 - (t01 * IL1_01 + t12 * IL1_12 + t20 * IL1_20));
  *ILn1 = ILn1_01 * sh.m01 + ILn1_12 * sh.m12 + ILn1_20 * sh.m20;
  *IL1 = IL1_01 * sh.m01 + IL1_12 * sh.m12 + IL1_20 * sh.m20;
  *IL3 = IL3_01 * sh.m01 + IL3_12 * sh.m12 + IL3_20 * sh.m20;
}

class EHSmoothIntegrand : public TriangleIntegrand {
public:
  ComplexVector3 DG;
  ComplexVector3 Gphi[3];
  ComplexVector3 DGxphi[3];
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  EHSmoothIntegrand(GreensFunction *g, vector<RWGIndexPair> &s) : g(g), s(s) {}
  
  void setTestPoint(Vector3 &v) {
    rt = v;
    memset(&DG,0,3*sizeof(Complex));
    memset(Gphi,0,9*sizeof(Complex));
    memset(DGxphi,0,9*sizeof(Complex));
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs;
    ComplexVector3 gg;
    g->scalar_gradient_smooth(&rt,&v,&gs,&gg);
    gs *= w;
    gg *= w;
    DG += gg; 
    for(unsigned int i=0; i<s.size(); i++) {
      Vector3 phi = 0.5 * (v - *(s[i].first)->p0);
      Gphi[i] += gs * phi; 
      DGxphi[i] += gg * phi;
    }
  }
};

void DipolePointSource :: EHTriangle(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options)
{
  EHSmoothIntegrand smooth(g,t);
  TriangleQuadratureGauss21 innerquad;
  smooth.setTestPoint(p);
  HalfRWG &th = *t[0].first;
  innerquad.integrate(smooth,*th.p0,th.v01,th.v02);
  
  Real ISn3h;
  Real ISn1;
  Real IS1;
  Vector3 ILn1;
  Vector3 IL1;
  Vector3 IL3;
  EHNonSmoothIntegral(th,p,&ISn3h,&ISn1,&IS1,&ILn1,&IL1,&IL3);
  Complex &k2 = g->k2;
  Real Ai = th.Ai;
  Real Ai2 = 0.5 * Ai;
  ComplexVector3 IL13 = Ai2 * (IL1 + (ONETHIRD * k2) * IL3);
  Complex ISn11 = Ai2 * (ISn1 + k2 * IS1);
  Vector3 &n = th.n;
  Real h = dot(p - *th.p0,n);
  Vector3 rho = p - h * n;
  ComplexVector3 K31 = Ai2 * ((ILn1 + ISn3h * n) + k2 * (IL1 - h * ISn1 * n));
  ComplexVector3 E1 = smooth.DG + Ai * ((ISn3h - k2 * h * ISn1) * n + (ILn1 + k2 * IL1));

  bool bExternal = g->isExternal();

  for(unsigned int i=0; i<t.size(); i++) {
    int ii = t[i].second<<1;
    Real c = ONEOVER4PI * t[i].first->e->length;
    if(bExternal != t[i].first->bPositive) {
      c = -c;
    }
    Complex Ei = c * g->Zeps * dot(j,E1 - g->ksquared * (smooth.Gphi[i] + IL13 + ISn11 * (rho - *t[i].first->p0)));
    Complex Hi = c * dot(j,-(smooth.DGxphi[i] + K31 * (p - *t[i].first->p0)));
    
    B(ii,0) += Ei;
    B(ii+1,0) += Hi;
  }
}

ComplexDenseMatrix *Bdense(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options)
{
  GreensFunction *g1 = t->g;
  GreensFunction *g2 = gExt;

  ComplexDenseMatrix *Bp = new ComplexDenseMatrix(2*t->rwgs,1,0);
  ComplexDenseMatrix &B = *Bp;

  TriangleRWGMap tt;

  for(int i=0; i<t->rwgs; i++) {
    RWG *ti = t->rwg[i];
    tt[ti->pos.t].push_back(RWGIndexPair(&ti->pos,i));
    tt[ti->neg.t].push_back(RWGIndexPair(&ti->neg,i));
  }

  for(SourceSet::iterator j=sources.begin(); j != sources.end(); ++j) {
    Source *s = *j;
    GreensFunction *g;
    if(s->isInRegion(g1)) {
      g = g1;
    } else if(s->isInRegion(g2)) {
      g = g2;
    } else {
      g = NULL;
    }
    if(g) {
      for(TriangleRWGMap::iterator i = tt.begin(); i != tt.end(); ++i) {
        s->EHTriangle(g,i->first,i->second,B,options);
      }
    }
  }
  // progress.finishedRHS(B);
  return Bp;
}

ComplexMatrix * BH(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options)
{
  if(!t->sub) {
    return Bdense(t,gExt,sources,progress,options);
  }
  ComplexHMatrix *Hp = new ComplexHMatrix(t->subs,1);
  ComplexHMatrix &H = *Hp;
  for(int i=0; i<t->subs; i++) {
    H.r[i] = 2*t->sub[i]->rwgs;
  }
  H.c[0] = 1;
  for(int i=0; i<t->subs; i++) {
    ComplexMatrix *M = BH(t->sub[i],gExt,sources,progress,options);
    H.set(i,0,M);
  }
  return Hp;
}

bool DipolePointSource :: isInRegion(GreensFunction *g)
{
  return (g == this->g);
}
