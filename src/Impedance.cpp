#include "Impedance.h"
#include "GreensFunction.h"
#include "FunctionArgs.h"
#include "TriangleQuadrature.h"
#include "Triangle.h"
#include "Options.h"
#include "MatrixOperations.h"
#include "ComplexVector3.h"
#include "Source.h"

//#define DEBUG

void LKEdgeIntegral(HalfRWG &th, const Vector3 &rs, Real *ISn1, Vector3 *IL1)
{
  Vector3 p0 = rs - *th.p0;
  Vector3 p1 = rs - *th.p1;
  Vector3 p2 = rs - *th.p2;

  Real R02 = norm(p0);
  Real R0 = sqrt(R02);
  Real R12 = norm(p1);
  Real R1 = sqrt(R12);
  Real R22 = norm(p2);
  Real R2 = sqrt(R22);
  Real tol = 1e-8 * th.e->length;

  Vector3 a0 = 1. / R0 * p0;
  Vector3 a1 = 1. / R1 * p1;
  Vector3 a2 = 1. / R2 * p2;

  Vector3 &n = th.n;
  Real h = dot(p0,n);

  Real x = 1.0 + dot(a0,a1+a2) + dot(a1,a2);
  Real y = sqrt(norm2(dot(a0,a1*a2)));
  Real ISn3h = (h<0?-2.:2.) * atan2(y,x);

  Real t01 = dot(th.m01,p0);
  Real t12 = dot(th.m12,p1);
  Real t20 = dot(th.m20,p2);

  Real s010 = th.normv01i * dot(th.v01,p0);
  Real s011 = th.normv01i * dot(th.v01,p1);
  Real num, den;
  num = (R1 - s011);
  den = (R0 - s010);
  Real ILn1_01 = fabs(num) < tol || fabs(den) < tol?0.:log(num / den);

  Real s121 = th.normv12i * dot(th.v12,p1);
  Real s122 = th.normv12i * dot(th.v12,p2);
  num = (R2 - s122);
  den = (R1 - s121);
  Real ILn1_12 = fabs(num) < tol || fabs(den) < tol?0.:log(num / den);

  Real s202 = th.normv20i * dot(th.v20,p2);
  Real s200 = th.normv20i * dot(th.v20,p0);
  num = (R0 - s200);
  den = (R2 - s202);
  Real ILn1_20 = fabs(num) < tol || fabs(den) < tol?0.:log(num / den);
  
  *ISn1 = -h * ISn3h - (t01 * ILn1_01 + t12 * ILn1_12 + t20 * ILn1_20);

  Real R02_01 = (R02 - s010 * s010);
  Real IL1_01 = 0.5 * (R02_01 * ILn1_01 - (R1 * s011 - R0 * s010));

  Real R02_12 = (R12 - s121 * s121);
  Real IL1_12 = 0.5 * (R02_12 * ILn1_12 - (R2 * s122 - R1 * s121));

  Real R02_20 = (R22 - s202 * s202);
  Real IL1_20 = 0.5 * (R02_20 * ILn1_20 - (R0 * s200 - R2 * s202));

  *IL1 = (IL1_01 * th.m01 + IL1_12 * th.m12 + IL1_20 * th.m20);
}

class LKSmoothInnerIntegrand : public TriangleIntegrand {
public:
  Complex G;
  ComplexVector3 Gphi[3];
  ComplexVector3 DGxphi[3];
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  LKSmoothInnerIntegrand(GreensFunction *g, vector<RWGIndexPair> &s) : g(g), s(s) {}
  
  void setTestPoint(Vector3 &v) {
    rt = v;
    G = 0.0;
    memset(Gphi,0,9*sizeof(Complex));
    memset(DGxphi,0,9*sizeof(Complex));
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs;
    ComplexVector3 gg;
    g->scalar_gradient_smooth(&rt,&v,&gs,&gg);
    gs *= w;
    gg *= w;
    G += gs;
    for(unsigned int i=0; i<s.size(); i++) {
      Vector3 phi = 0.5 * (v - *(s[i].first)->p0);
      Gphi[i] += gs * phi; 
      DGxphi[i] += gg * phi;
    }
  }
};

class LKInnerIntegrand : public TriangleIntegrand {
public:
  Complex G;
  ComplexVector3 Gphi[3];
  ComplexVector3 DGxphi[3];
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  LKInnerIntegrand(GreensFunction *g, vector<RWGIndexPair> &s) : g(g), s(s) {}
  
  void setTestPoint(Vector3 &v) {
    rt = v;
    G = 0.0;
    memset(Gphi,0,9*sizeof(Complex));
    memset(DGxphi,0,9*sizeof(Complex));
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs;
    ComplexVector3 gg;
    g->scalar_gradient(&rt,&v,&gs,&gg);
    gs *= w;
    gg *= w;
    G += gs;
    for(unsigned int i=0; i<s.size(); i++) {
      Vector3 phi = 0.5 * (v - *(s[i].first)->p0);
      Gphi[i] += gs * phi; 
      DGxphi[i] += gg * phi;
    }
  }
};

class LKFarOuterIntegrand : public TriangleIntegrand {
public:
  GreensFunction *g;
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  Complex *L;
  Complex *K;
  LKInnerIntegrand inner;
  TriangleQuadrature &innerquad;

  LKFarOuterIntegrand(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, Complex *L, Complex *K, TriangleQuadrature &innerquad) 
    : g(g),t(t),s(s),L(L),K(K),inner(g,s),innerquad(innerquad) {}

  void weighPoint(Real w, Vector3 &v) {
    inner.setTestPoint(v);
    HalfRWG &sh = *s[0].first;
    innerquad.integrate(inner,*sh.p0,sh.v01,sh.v02);
    int ij = 0;

    for(unsigned int i=0; i<t.size(); i++) {
      Vector3 phi = 0.5 * (v - *(t[i].first->p0));
      for(unsigned int j=0; j<s.size(); j++) {
        L[ij] += w * (g->ksquared * dot(phi,inner.Gphi[j]) - inner.G);
        K[ij] -= w * dot(phi,inner.DGxphi[j]);
        ij++;
      }
    }
  }  
};

class LKCloseOuterIntegrand : public TriangleIntegrand {
public:
  GreensFunction *g;
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  Complex *L;
  Complex *K;
  LKSmoothInnerIntegrand inner_smooth;
  TriangleQuadratureGauss91 innerquad;
  Vector3 vts[9];
  bool bInPlane;

  LKCloseOuterIntegrand(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, Complex *L, Complex *K) 
    : g(g),t(t),s(s),L(L),K(K),inner_smooth(g,s) {
    int ij = 0;
    for(unsigned int i=0; i<t.size(); i++) {
      for(unsigned int j=0; j<s.size(); j++) {
        vts[ij] = *s[j].first->p0 - *t[i].first->p0;
        ij++;
      }
    }
  }
  
  void weighPoint(Real w, Vector3 &v) {
    inner_smooth.setTestPoint(v);
    HalfRWG &sh = *s[0].first;
    innerquad.integrate(inner_smooth,*sh.p0,sh.v01,sh.v02);
    Real ISn3h;
    Real ISn1;
    Real IS1;
    Vector3 ILn1;
    Vector3 IL1;
    Vector3 IL3;
    EHNonSmoothIntegral(sh,v,&ISn3h,&ISn1,&IS1,&ILn1,&IL1,&IL3);
    Complex &k2 = g->k2;
    Real Ai = sh.Ai;
    Real Ai2 = 0.5 * Ai;
    ComplexVector3 IL13 = Ai2 * (IL1 + (ONETHIRD * k2) * IL3);
    Complex ISn11 = Ai2 * (ISn1 + k2 * IS1);
    Vector3 &n = sh.n;
    Real h = dot(v - *sh.p0,n);
    if(fabs(h) < 1e-8 * sh.e->length) {
      ComplexVector3 K21[3];
      for(unsigned int j=0; j<s.size(); j++) {
        K21[j] = IL13 + ISn11 * (v - *s[j].first->p0);
      }
      Complex L2 = 2.0 * ISn11 + inner_smooth.G;
      int ij = 0;
      for(unsigned int i=0; i<t.size(); i++) {
        Vector3 phi = 0.5 * (v - *t[i].first->p0);
        for(unsigned int j=0; j<s.size(); j++) {
          L[ij] += w * (g->ksquared * dot(phi,inner_smooth.Gphi[j] + K21[j]) - L2);
          if(t[i].first->t == s[j].first->t) {
            //K[ij] += w * Ai2 * (g->isExternal()?TWOPI:-TWOPI) * dot(phi,n * (v - *s[j].first->p0));
          }
          ij++;
        }
      }
    } else {
      Vector3 ISn3hn = Ai2 * ISn3h * n;
      ComplexVector3 K21[3];
      ComplexVector3 K41[3];
      Vector3 rho = v - h * n;
      ComplexVector3 K31 = Ai2 * k2 * (IL1 - h * ISn1 * n);
      for(unsigned int j=0; j<s.size(); j++) {
        K21[j] = IL13 + ISn11 * (rho - *s[j].first->p0);
        K41[j] = K31 * (v - *s[j].first->p0);
      }
      Complex L2 = 2.0 * ISn11 + inner_smooth.G;
      int ij = 0;
      for(unsigned int i=0; i<t.size(); i++) {
        Vector3 phi = 0.5 * (v - *(t[i].first)->p0);
        for(unsigned int j=0; j<s.size(); j++) {
          L[ij] += w * (g->ksquared * dot(phi,inner_smooth.Gphi[j] + K21[j]) - L2);
          K[ij] -= w * dot(phi,inner_smooth.DGxphi[j] + vts[ij] * ISn3hn + K41[j]);
          ij++;
        }
      }
    }
  }
};

class LKNeighborEdgeIntegrand : public Integrand {
public:
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  Vector3 *K;

  LKNeighborEdgeIntegrand(vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, Vector3 *K) : t(t), s(s), K(K) {}

  void weighPoint(Real w, Vector3 &v) {
    Real ISn1;
    Vector3 IL1;
    HalfRWG &th = *t[0].first;
    LKEdgeIntegral(th,v,&ISn1,&IL1);
    Vector3 &n = th.n;
    Real h = dot(v - *th.p0,n);
    Vector3 rho = v - h * n;
    for(unsigned int i=0; i<t.size(); i++) {
      Vector3 rho_p0 = (rho - *(t[i].first)->p0);
      K[i] += w * (IL1 + ISn1 * rho_p0);
    }
  }
};

void LKTriangle(GreensFunction *g, Triangle *trit, Triangle *tris, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, ComplexDenseMatrix &Z, Options &options)
{
  Complex L[9];
  Complex K[9];
  memset(L,0,9*sizeof(Complex));
  memset(K,0,9*sizeof(Complex));

  HalfRWG &th = *(t[0].first);

  if(trit == tris) {
#ifdef DEBUG
    printf("Same\n");
#endif
    TriangleQuadratureGauss105 outerquad;
    LKCloseOuterIntegrand outerintegrand(g,t,s,L,K);
    outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
  } else {
    int shared = sharedpoints(trit,tris);
    if(shared > 0) {
#ifdef DEBUG
      printf("Neighbors\n");
#endif
      TriangleQuadratureGauss105 outerquad;
      LKCloseOuterIntegrand outerintegrand(g,t,s,L,K);
      outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      QuadratureGauss15 quad;
      Vector3 K01[3];
      Vector3 K12[3];
      Vector3 K20[3];
      memset(K01,0,9*sizeof(Real));
      memset(K12,0,9*sizeof(Real));
      memset(K20,0,9*sizeof(Real));
      HalfRWG &sh = *(s[0].first);
      LKNeighborEdgeIntegrand integrand01(t,s,K01);
      quad.integrate(integrand01,*sh.p0,sh.v01);
      LKNeighborEdgeIntegrand integrand12(t,s,K12);
      quad.integrate(integrand12,*sh.p1,sh.v12);
      LKNeighborEdgeIntegrand integrand20(t,s,K20);
      quad.integrate(integrand20,*sh.p2,sh.v20);

      int ij = 0;
      Real c = 0.125 * t[0].first->Ai * s[0].first->Ai;
      for(unsigned int i=0; i<t.size(); i++) {
        for(unsigned int j=0; j<s.size(); j++) {
          Vector3 v = *t[i].first->p0 - *s[j].first->p0;
          K[ij] += c * (dot(v * sh.m01,K01[i]) / sh.normv01i + dot(v * sh.m12,K12[i]) / sh.normv12i + dot(v * sh.m20,K20[i]) / sh.normv20i);
          ij++;
        }
      }
    } else {
#ifdef DEBUG
      printf("Far\n");
#endif
      Point pt = center(trit);
      Point ps = center(tris);
      Vector3 v = pt - ps;
      Real d = norm(v);
      Real temax = max(max(trit->e12->length,trit->e23->length),trit->e31->length);
      Real semax = max(max(tris->e12->length,tris->e23->length),tris->e31->length);
      Real emax2 = Square(max(temax,semax));
      if(d * options.farThresh6 > emax2) {
        TriangleQuadratureGauss6 innerquad;
        TriangleQuadratureGauss6 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh10 > emax2) {
        TriangleQuadratureGauss10 innerquad;
        TriangleQuadratureGauss10 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh15 > emax2) {
        TriangleQuadratureGauss15 innerquad;
        TriangleQuadratureGauss15 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh21 > emax2) {
        TriangleQuadratureGauss21 innerquad;
        TriangleQuadratureGauss21 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else {
        TriangleQuadratureGauss36 innerquad;
        TriangleQuadratureGauss36 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      }
    }
  }

  int ij = 0;
  for(unsigned int i=0; i<t.size(); i++) {
    HalfRWG &th = *(t[i].first);
    int ii = t[i].second<<1;
    for(unsigned int j=0; j<s.size(); j++) {
      HalfRWG &sh = *(s[j].first);
      int jj = s[j].second<<1;
      Real c = ONEOVER4PI * th.e->length * sh.e->length;
      if(th.bPositive != sh.bPositive) {
        c = -c;
      }
      Complex Lij = c * L[ij];
      Complex Kij = c * K[ij];
      Z(ii,jj) += g->Zeps * Lij;
      Z(ii,jj+1) += Kij;
      Z(ii+1,jj) -= Kij;
      Z(ii+1,jj+1) += g->Zmu * Lij;
      ij++;
    }
  }
}


ComplexDenseMatrix *Zts(RWGTree *t, 
                        RWGTree *s,
                        GreensFunction *g,
                        Options &options)
{
  TriangleRWGMap tt;
  TriangleRWGMap ss;
  ComplexDenseMatrix *dense = new ComplexDenseMatrix(2*t->rwgs,2*s->rwgs);
  dense->zero();
  for(int i=0; i<t->rwgs; i++) {
    RWG *ti = t->rwg[i];
    tt[ti->pos.t].push_back(RWGIndexPair(&ti->pos,i));
    tt[ti->neg.t].push_back(RWGIndexPair(&ti->neg,i));
  }
  for(int i=0; i<s->rwgs; i++) {
    RWG *si = s->rwg[i];
    ss[si->pos.t].push_back(RWGIndexPair(&si->pos,i));
    ss[si->neg.t].push_back(RWGIndexPair(&si->neg,i));
  }
  for(TriangleRWGMap::iterator i = tt.begin(); i != tt.end(); ++i) {
    for(TriangleRWGMap::iterator j = ss.begin(); j != ss.end(); ++j) {
      LKTriangle(g,i->first,j->first,i->second,j->second,*dense,options);
    }
  }
  return dense;
}

void Zdense(RWGTree *t, 
            RWGTree *s,
            ComplexDenseMatrix *dense,
            GreensFunction *gExt,
            Options &options)
{
  TriangleRWGMap tt;
  TriangleRWGMap ss;
  dense->resize(2*t->rwgs,2*s->rwgs);
  dense->zero();

  for(int i=0; i<t->rwgs; i++) {
    RWG *ti = t->rwg[i];
    tt[ti->pos.t].push_back(RWGIndexPair(&ti->pos,i));
    tt[ti->neg.t].push_back(RWGIndexPair(&ti->neg,i));
  }
  for(int i=0; i<s->rwgs; i++) {
    RWG *si = s->rwg[i];
    ss[si->pos.t].push_back(RWGIndexPair(&si->pos,i));
    ss[si->neg.t].push_back(RWGIndexPair(&si->neg,i));
  }
  bool bConnected = (t->g == s->g);
  for(TriangleRWGMap::iterator i = tt.begin(); i != tt.end(); ++i) {
    for(TriangleRWGMap::iterator j = ss.begin(); j != ss.end(); ++j) {
      LKTriangle(gExt,i->first,j->first,i->second,j->second,*dense,options);
      if(bConnected) {
        LKTriangle(t->g,i->first,j->first,i->second,j->second,*dense,options);
      }
    }
  }
}
