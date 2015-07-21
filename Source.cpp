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


DipolePointSource :: DipolePointSource(GreensFunction *g_, const Point &p_, Vector3 &j_) : g(g_), p(p_), j(j_) {}

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
  TriangleQuadratureGauss91 innerquad;
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



class EHSmoothIntegrandLongitudinal : public TriangleIntegrand {
public:
  ComplexVector3 DG;
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  EHSmoothIntegrandLongitudinal(GreensFunction *g, vector<RWGIndexPair> &s) : g(g), s(s) {}
  
  void setTestPoint(Vector3 &v) {
    rt = v;
    memset(&DG,0,3*sizeof(Complex));
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs;
    ComplexVector3 gg;
    g->scalar_gradient_smooth(&rt,&v,&gs,&gg);
    gg *= w;
    DG += gg; 
  }
};

void DipolePointSource :: EHTriangleLongitudinal(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options)
{
  EHSmoothIntegrandLongitudinal smooth(g,t);
  TriangleQuadratureGauss91 innerquad;
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
  Vector3 &n = th.n;
  Real h = dot(p - *th.p0,n);
  ComplexVector3 E1 = smooth.DG + Ai * ((ISn3h - k2 * h * ISn1) * n + (ILn1 + k2 * IL1));

  bool bExternal = g->isExternal();

  for(unsigned int i=0; i<t.size(); i++) {
    int ii = t[i].second<<1;
    Real c = ONEOVER4PI * t[i].first->e->length;
    if(bExternal != t[i].first->bPositive) {
      c = -c;
    }
    Complex Ei = c * g->Zeps * dot(j,E1);
    //B(ii,0) -= Ei;    
  }
}

ComplexVector3Pair getField(const SourceSet &sources, RWGTree *t, const ComplexDenseMatrix &JM, GreensFunction *g, const Point &p, bool bLong, Options &options)
{
  ComplexVector3Pair EH(ComplexVector3(0,0,0),ComplexVector3(0,0,0));
  if((!bLong) && g->isExternal()) {
    for(SourceSet::const_iterator j=sources.begin(); j != sources.end(); ++j) {
      Source *source = *j;
      EH += source->getField(p);
    }
  }

  Vector3 x(1,0,0);
  Vector3 y(0,1,0);
  Vector3 z(0,0,1);
  DipolePointSource sx(g,p,x);
  DipolePointSource sy(g,p,y);
  DipolePointSource sz(g,p,z);

  ComplexDenseMatrix Bx(2*t->rwgs,1,0);
  ComplexDenseMatrix By(2*t->rwgs,1,0);
  ComplexDenseMatrix Bz(2*t->rwgs,1,0);
  TriangleRWGMap tt;

  for(int i=0; i<t->rwgs; i++) {
    RWG *ti = t->rwg[i];
    tt[ti->pos.t].push_back(RWGIndexPair(&ti->pos,i));
    tt[ti->neg.t].push_back(RWGIndexPair(&ti->neg,i));
  }
  
  for(TriangleRWGMap::iterator i = tt.begin(); i != tt.end(); ++i) {

    if(bLong) {
      if(g->gLong) {
        sx.EHTriangleLongitudinal(g->gLong,i->first,i->second,Bx,options);
        sy.EHTriangleLongitudinal(g->gLong,i->first,i->second,By,options);
        sz.EHTriangleLongitudinal(g->gLong,i->first,i->second,Bz,options);
      }
    } else {
      sx.EHTriangle(g,i->first,i->second,Bx,options);
      sy.EHTriangle(g,i->first,i->second,By,options);
      sz.EHTriangle(g,i->first,i->second,Bz,options);
    }
  }
  
  ComplexVector3 E(0,0,0);
  ComplexVector3 H(0,0,0);

  Complex Z2 = g->eps / g->mu;
  for(int i=0; i<t->rwgs; i++) {
    int ii = (i<<1);
    E.x += JM(ii,0) * Bx(ii,0) - JM(ii+1,0) * Bx(ii+1,0);
    E.y += JM(ii,0) * By(ii,0) - JM(ii+1,0) * By(ii+1,0);
    E.z += JM(ii,0) * Bz(ii,0) - JM(ii+1,0) * Bz(ii+1,0);

    H.x += JM(ii,0) * Bx(ii+1,0) + Z2 * JM(ii+1,0) * Bx(ii,0);
    H.y += JM(ii,0) * By(ii+1,0) + Z2 * JM(ii+1,0) * By(ii,0);
    H.z += JM(ii,0) * Bz(ii+1,0) + Z2 * JM(ii+1,0) * Bz(ii,0);
  }
  EH += ComplexVector3Pair(E,H);
  return EH;
}

PlaneWaveSource :: PlaneWaveSource(GreensFunction *g, const Vector3 &n, const ComplexVector3 &_E0) : g(g), k(g->k), n(n), E0(_E0) 
{
  E0 *= 1.0 / sqrt(norm(E0));
  H0 = k * (n * E0);
  H0 *= 1.0 / (g->omega * g->mu);
}


class PlaneWaveIntegrand : public TriangleIntegrand {
public:
  Complex E[3];
  Complex H[3];
  PlaneWaveSource *src;
  vector<RWGIndexPair> &s;

  PlaneWaveIntegrand(PlaneWaveSource *src, vector<RWGIndexPair> &s) : src(src), s(s) {
    memset(E,0,3*sizeof(Complex));    
    memset(H,0,3*sizeof(Complex));    
  }
  
  void weighPoint(Real w, Vector3 &v) { 

    Complex I(0,1);
    Complex ph = exp(I * src->k * dot(v,src->n));
    ComplexVector3 E1 = src->E0 * ph;
    ComplexVector3 H1 = src->H0 * ph;    

    for(unsigned int i=0; i<s.size(); i++) {
      Vector3 phi = 0.5 * (v - *(s[i].first)->p0);
      E[i] += w * dot(phi, E1); 
      H[i] += w * dot(phi, H1);
    }
  }
};

void PlaneWaveSource :: EHTriangle(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options)
{
  PlaneWaveIntegrand integrand(this,t);
  HalfRWG &th = *t[0].first;
  TriangleQuadratureGauss91 quad;
  quad.integrate(integrand,*th.p0,th.v01,th.v02);

  bool bExternal = g->isExternal();

  for(unsigned int i=0; i<t.size(); i++) {
    int ii = t[i].second<<1;
    Real c = t[i].first->e->length;
    if(bExternal != t[i].first->bPositive) {
      c = -c;
    }
    Complex Ei = c * integrand.E[i];
    Complex Hi = c * integrand.H[i];
    
    B(ii,0) += Ei;
    B(ii+1,0) += Hi;
  }
}

ComplexVector3Pair DipolePointSource :: getField(const Point &p)
{
  abort();  
  ComplexVector3 E;
  ComplexVector3 H;

  return ComplexVector3Pair(E,H);
}

ComplexVector3Pair PlaneWaveSource :: getField(const Point &p)
{
  Complex I(0,1);
  Complex ph = exp(I * g->k * dot(n, p));
  ComplexVector3 E = ph * E0;
  ComplexVector3 H = ph * H0;

  return ComplexVector3Pair(E,H);
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

bool PlaneWaveSource :: isInRegion(GreensFunction *g)
{
  return (g == this->g);
}

Real extinctionCrossSection(RWGTree *t, const ComplexDenseMatrix &JM, GreensFunction *g, PlaneWaveSource *plane)
{
  Complex I(0,1);
  ComplexVector3 J(0,0,0);
  ComplexVector3 M(0,0,0);

  for(int i=0; i<t->rwgs; i++) {
    RWG *rwg = t->rwg[i];
    int ii = i<<1;
    Complex ph1 = exp(-I * g->k * dot(plane->n, center(rwg->pos.t)));
    Complex ph2 = exp(-I * g->k * dot(plane->n, center(rwg->neg.t)));
    ComplexVector3 f = ONESIXTH * rwg->e->length * ( ph1 * (rwg->pos.v01 + rwg->pos.v02) - ph2 * (rwg->neg.v01 + rwg->neg.v02));
    J += JM(ii,0) * f;
    M += JM(ii+1,0) * f;
  }
  // 1 / 4pi cancels in optical theorem
  ComplexVector3 E = I * g->omega * g->mu * (plane->n * J) * plane->n - I * g->k * (plane->n * M);

  Complex f0 = dot(E, plane->E0) / dot(plane->E0,plane->E0);
  Real sigma = imag(f0) / g->normk;

  return sigma;
}
