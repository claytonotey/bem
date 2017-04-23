#include "ImpedanceNew.h"
#include "GreensFunction.h"
#include "FunctionArgs.h"
#include "TriangleQuadrature.h"
#include "Triangle.h"
#include "Options.h"
#include "MatrixOperations.h"
#include "ComplexVector3.h"
#include "Source.h"
#include "Nonlocal.h"

//#define DEBUG

/* n x RWG functions cawnnot be used as both basis and testing functions, as the results are
   divergent.  For the present purposes, the testing functions are always RWG, and the basis
   functions may be RWG or n x RWG.  For ease of computation, the <RWG | LK | n x RWG> elements
   are calculated using 
  <RWG | L | n x RWG> = <n x RWG | L | RWG> 
  <RWG | K | n x RWG> = -<n x RWG | K | RWG> 
  So the same inner integrals may be used, K1, K2 K3 and the testing is done by
  explicitly taking the dot product with <n x RWG| instead of <RWG| for the smooth terms,
  or by converting the <n x RWG| integration into a line integral,
  or in the case of K4, by using the line integral from Taskinen.
*/

void LKEdgeIntegral(HalfRWG &th, const Vector3 &rs, Real *ISn1, Vector3 *IL1, Real *IS1 = NULL)
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
  Real tol = 1e-6 * th.e->length;

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

  if(IS1) *IS1 = ONETHIRD * (h * h * *ISn1 - (t01 * IL1_01 + t12 * IL1_12 + t20 * IL1_20));
}

class LKSmoothInnerIntegrand : public TriangleIntegrand {
public:
  Complex G;
  ComplexVector3 Gphi[3];
  ComplexVector3 DGxphi[3];
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  LKSmoothInnerIntegrand(GreensFunction *g, vector<RWGIndexPair> &s, int basisType) : g(g), s(s) {}
  
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
  Complex DGdotphi[3];
  ComplexVector3 DGxphi[3];
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;
  int basis;

  // basisTypeT is the test basis, not the source basis
  // if the test basis is n x RWG, then DGdotphi is required
  // otherwise if RWG is the test basis, then DG|RWG> reduced to G D dot |RWG>
  // and G is used
  LKInnerIntegrand(GreensFunction *g, vector<RWGIndexPair> &s, int basisTypeT) : g(g), s(s), basis(basisTypeT)  {}
  
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
      if(basis == nCrossRWGBasis) {
        phi = (s[i].first)->n * phi;
        DGdotphi[i] += dot(gg,phi);
      }
      Gphi[i] += gs * phi; 
      DGxphi[i] += gg * phi;
    }
  }
};

class LK1Integrand : public TriangleIntegrand {
public:
  Vector3 &rt;
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  ComplexVector3 Gphi[3];
  int basis;

  LK1Integrand(Vector3 &rt, GreensFunction *g, vector<RWGIndexPair> &s, int basisType) : rt(rt), g(g), s(s), basis(basisType) {
    memset(Gphi,0,9*sizeof(Complex));
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs = g->scalar(&rt,&v);
    gs *= w;
    for(unsigned int i=0; i<s.size(); i++) {
      Vector3 phi = (v - *(s[i].first)->p0);
      if(basis == nCrossRWGBasis) {
        phi = (s[i].first)->n * phi;
      }
      Gphi[i] += gs * phi;
    }
  }
};

class LK2Integrand : public TriangleIntegrand {
public:
  Vector3 &rt;
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  ComplexVector3 L[3];
  ComplexVector3 K[3];
  int basis;

  LK2Integrand(Vector3 &rt, GreensFunction *g, vector<RWGIndexPair> &s, int basisType = RWGBasis) : rt(rt), g(g), s(s), basis(basisType) {
    memset(L,0,9*sizeof(Complex));
    memset(K,0,9*sizeof(Complex));
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs;
    ComplexVector3 gg;
    g->scalar_gradient(&rt,&v,&gs,&gg);
    gs *= w;
    gg *= w;
    for(unsigned int i=0; i<s.size(); i++) {
      Vector3 phi = (v - *(s[i].first)->p0);
      if(basis == nCrossRWGBasis) {
        phi = (s[i].first)->n * phi;
      }
      L[i] += g->ksquared * gs * phi - 2. * gg;
      K[i] -= gg * phi;
    }
  }
};

// far integrand should work with n x f.  the DD term can be treatedXS with K5

class LKFarOuterIntegrand : public TriangleIntegrand {
public:
  GreensFunction *g;
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  int basis;
  Complex *L;
  Complex *K;
  LKInnerIntegrand inner;
  TriangleQuadrature &innerquad;

  LKFarOuterIntegrand(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, int basisType, Complex *L, Complex *K, TriangleQuadrature &innerquad) 
    : g(g),t(t),s(s),basis(basisType),L(L),K(K),inner(g,s,basisType),innerquad(innerquad) {}

  void weighPoint(Real w, Vector3 &v) {
    inner.setTestPoint(v);
    HalfRWG &sh = *s[0].first;
    innerquad.integrate(inner,*sh.p0,sh.v01,sh.v02);
    int ij = 0;

    for(unsigned int i=0; i<t.size(); i++) {
      Vector3 phi = 0.5 * (v - *(t[i].first->p0));
      if(basis == nCrossRWGBasis) {
        phi = (s[i].first)->n * phi;
      }
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
  int basis;
  Complex *L;
  Complex *K;
  LKSmoothInnerIntegrand inner_smooth;
  TriangleQuadratureGauss91 innerquad;
  
  LKCloseOuterIntegrand(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, int basisType, Complex *L, Complex *K) 
    : g(g),t(t),s(s),basis(basisType),L(L),K(K),inner_smooth(g,s,basisType) {
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

    if(basis == nCrossRWGBasis) {
      HalfRWG &th = *t[0].first;
      if(fabs(h) < 1e-8 * sh.e->length) {
        ComplexVector3 K21[3];
        for(unsigned int j=0; j<s.size(); j++) {
          K21[j] = IL13 + ISn11 * (v - *s[j].first->p0);
        }
        int ij = 0;
        for(unsigned int i=0; i<t.size(); i++) {
          Vector3 phi = 0.5 * (v - *t[i].first->p0);
          phi = th.n * phi;
          for(unsigned int j=0; j<s.size(); j++) {
            L[ij] += w * (g->ksquared * dot(phi,inner_smooth.Gphi[j] + K21[j]));
            if(t[i].first->t == s[j].first->t) {
              //K[ij] += w * Ai2 * (g->isExternal()?TWOPI:-TWOPI) * dot(phi,n * (v - *s[j].first->p0));
            }
            ij++;
          }
        }
      } else {
        ComplexVector3 K21[3];
        ComplexVector3 K41[3];
        Vector3 rho = v - h * n;
        ComplexVector3 K31 = Ai2 * (-ISn3h * n + k2 * (IL1 - h * ISn1 * n));
        for(unsigned int j=0; j<s.size(); j++) {
          K21[j] = IL13 + ISn11 * (rho - *s[j].first->p0);
          K41[j] = K31 * (v - *s[j].first->p0);
        }
        int ij = 0;
        for(unsigned int i=0; i<t.size(); i++) {
          Vector3 phi = 0.5 * (v - *(t[i].first)->p0);
          phi = th.n * phi;
          for(unsigned int j=0; j<s.size(); j++) {
            L[ij] += w * (g->ksquared * dot(phi,inner_smooth.Gphi[j] + K21[j]));
            K[ij] -= w * dot(phi,inner_smooth.DGxphi[j] + K41[j]);
            ij++;
          }
        }
      }
    }

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
      ComplexVector3 K21[3];
      ComplexVector3 K41[3];
      Vector3 rho = v - h * n;
      // XXX Isn33h term sign
      ComplexVector3 K31 = Ai2 * (-ISn3h * n + k2 * (IL1 - h * ISn1 * n));
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
          K[ij] -= w * dot(phi,inner_smooth.DGxphi[j] + K41[j]);
          ij++;
        }
      }
    }
  }
};

class LKEdgeIntegrand : public Integrand {
public:
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  Vector3 *K;

  LKEdgeIntegrand(vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, int basisType, Vector3 *K) : t(t), s(s), K(K) {}

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

class nCrossLKEdgeIntegrand : public Integrand {
public:
  GreensFunction *g;
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  ComplexVector3 *L;
  Vector3 *K;
  Complex *Kn;
  int basis;
  LKSmoothInnerIntegrand inner_smooth;
  TriangleQuadratureGauss91 innerquad;

  nCrossLKEdgeIntegrand(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, ComplexVector3 *L, Vector3 *K, Complex *Kn) : t(t), s(s), L(L), K(K), Kn(Kn), inner_smooth(g,s,nCrossRWGBasis)  {}

  void weighPoint(Real w, Vector3 &v) {
    inner_smooth.setTestPoint(v);
    HalfRWG &sh = *s[0].first;
    innerquad.integrate(inner_smooth,*sh.p0,sh.v01,sh.v02);

    Real ISn1;
    Real IS1;
    Vector3 IL1;
    HalfRWG &th = *t[0].first;
    LKEdgeIntegral(th,v,&ISn1,&IL1,&IS1);

    Complex &k2 = g->k2;
    Real Ai = sh.Ai;
    Real Ai2 = 0.5 * Ai;
    Complex ISn11 = Ai2 * (ISn1 + k2 * IS1);
    Complex L2 = 2.0 * ISn11 + inner_smooth.G;      

    Vector3 vpn = v - *(sh.p0);

    Vector3 &n = th.n;
    Real h = dot(v - *th.p0,n);
    Vector3 rho = v - h * n;
    for(unsigned int i=0; i<t.size(); i++) {
      Vector3 rho_p0 = (rho - *(t[i].first)->p0);
      Vector3 vpm = v - *(t[i].first)->p0;
      L[i] += w * L2 * vpn;
      K[i] += w * (IL1 + ISn1 * rho_p0);
      // XXX Isn1 normal vector
      Kn[i] += w * (IS1 - norm(vpm) * ISn1 + 2. * dot(vpm, IL1 + rho_p0 + ISn1 * n));
    }
  }
};


void LKTriangle(GreensFunction *g, Triangle *trit, Triangle *tris, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, int basis, ComplexDenseMatrix &Z, Options &options)
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
    LKCloseOuterIntegrand outerintegrand(g,t,s,basis,L,K);
    outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
  } else {
    int shared = sharedpoints(trit,tris);
    if(shared > 0) {
#ifdef DEBUG
      printf("Neighbors\n");
#endif
      TriangleQuadratureGauss105 outerquad;
      LKCloseOuterIntegrand outerintegrand(g,t,s,basis,L,K);
      outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
 
      QuadratureGauss15 quad;
      ComplexVector3 L01[3];
      ComplexVector3 L12[3];
      ComplexVector3 L20[3];
      memset(L01,0,9*sizeof(Complex));
      memset(L12,0,9*sizeof(Complex));
      memset(L20,0,9*sizeof(Complex));
      Vector3 K01[3];
      Vector3 K12[3];
      Vector3 K20[3];
      memset(K01,0,9*sizeof(Real));
      memset(K12,0,9*sizeof(Real));
      memset(K20,0,9*sizeof(Real));
      Complex Kn01[3];
      Complex Kn12[3];
      Complex Kn20[3];
      memset(Kn01,0,3*sizeof(Complex));
      memset(Kn12,0,3*sizeof(Complex));
      memset(Kn20,0,3*sizeof(Complex));

      HalfRWG &sh = *(s[0].first);

      if(basis == RWGBasis) {
        LKEdgeIntegrand integrand01(t,s,basis,K01);
        quad.integrate(integrand01,*sh.p0,sh.v01);
        LKEdgeIntegrand integrand12(t,s,basis,K12);
        quad.integrate(integrand12,*sh.p1,sh.v12);
        LKEdgeIntegrand integrand20(t,s,basis,K20);
        quad.integrate(integrand20,*sh.p2,sh.v20);
      } else if(basis == nCrossRWGBasis){ 
        nCrossLKEdgeIntegrand integrand01(g,t,s,L01,K01,Kn01);
        quad.integrate(integrand01,*sh.p0,sh.v01);
        nCrossLKEdgeIntegrand integrand12(g,t,s,L12,K12,Kn12);
        quad.integrate(integrand12,*sh.p1,sh.v12);
        nCrossLKEdgeIntegrand integrand20(g,
t,s,L20,K20,Kn20);
        quad.integrate(integrand20,*sh.p2,sh.v20);
      }

      int ij = 0;
      Real c = 0.25 * t[0].first->Ai * s[0].first->Ai;
      for(unsigned int i=0; i<t.size(); i++) {
        for(unsigned int j=0; j<s.size(); j++) {
          Vector3 v = *t[i].first->p0 - *s[j].first->p0;
          if(basis == RWGBasis) {
            K[ij] += c * (dot(v * sh.m01,K01[i]) / sh.normv01i + 
                          dot(v * sh.m12,K12[i]) / sh.normv12i + 
                          dot(v * sh.m20,K20[i]) / sh.normv20i);
          } else {
            L[ij] += c * (dot(sh.v01, L01[i]) + dot(sh.v12, L12[i]) + dot(sh.v20, L20[i]));
            K[ij] += c * ((dot(v * sh.m01,th.n * K01[i]) - dot(sh.m01,th.n) * Kn01[i]) / sh.normv01i +
                          (dot(v * sh.m12,th.n * K12[i]) - dot(sh.m12,th.n) * Kn12[i]) / sh.normv12i + 
                          (dot(v * sh.m20,th.n * K20[i]) - dot(sh.m20,th.n) * Kn20[i]) / sh.normv20i);
          }
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
        LKFarOuterIntegrand outerintegrand(g,t,s,basis,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh10 > emax2) {
        TriangleQuadratureGauss10 innerquad;
        TriangleQuadratureGauss10 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,basis,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh15 > emax2) {
        TriangleQuadratureGauss15 innerquad;
        TriangleQuadratureGauss15 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,basis,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh21 > emax2) {
        TriangleQuadratureGauss21 innerquad;
        TriangleQuadratureGauss21 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,basis,L,K,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else {
        TriangleQuadratureGauss36 innerquad;
        TriangleQuadratureGauss36 outerquad;
        LKFarOuterIntegrand outerintegrand(g,t,s,basis,L,K,innerquad);
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
      if(isnan(Lij) || isnan(Kij)) {
        cout << *(th.t) << " " << *(sh.t);
        cout << Lij << " " << Kij << "\n";
        abort();
      }
      
      Z(ii,jj) += g->Zeps * Lij;
      Z(ii,jj+1) += Kij;
      Z(ii+1,jj) -= Kij;
      Z(ii+1,jj+1) += g->Zmu * Lij;
      
      ij++;
    
    }
  }
}

void Zdense(RWGTree *t, 
            RWGTree *s,
            ComplexDenseMatrix *dense,
            GreensFunction *g,
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
      LKTriangle(g,i->first,j->first,i->second,j->second,RWGBasis,*dense,options);
      if(bConnected) {
        //LKTriangle(t->g,i->first,j->first,i->second,j->second,*dense,options);
        if(t->g->gLong) {
          LKTriangle(t->g->gLong,i->first,j->first,i->second,j->second,nCrossRWGBasis,*dense,options);
          //LKTriangleLongitudinal(t->g->gLong,i->first,j->first,i->second,j->second,nCrossRWGBasis,*dense,options);
        }
      }
    }
  }
}

// get an array of 3d points distributed over the surface described by the surface parameterization
void getSurfacePoints(int n, SurfaceParameterization &u, Point **point, int *points)
{
  *points = SurfaceParameterization :: hexagonalPoints(n);
  *point = new Point[*points];

  int k = 0;
  for(int p = -n+1; p < n; p++) {
    int q0;
    int q1;
    if(p < 0) q0 = -n - p;
    else q0 = -n;
    if(p > 0) q1 = n - p;
    else q1 = n;
    for(int q = q0+1; q < q1; q++) {
      Point2 p2 = SurfaceParameterization :: hexagonalPoint(n,p,q);
      (*point)[k++] = u.inverse(p2);
    }
  }
}

// get all the edge midpoints
void getSurfacePoints(RWGTree *t, Point **point, int  *points)
{
  *points = (int) t->rwgs;
  *point = new Point[*points];
  
  for(int i=0; i<t->rwgs; i++) {
    RWG *rwg = t->rwg[i];
    Point p = *(rwg->e->p1);
    p += *(rwg->e->p2);
    p *= 0.5;
    (*point)[i] = p;
  }
}

void Zrk(RWGTree *t,
         RWGTree *s,
         ComplexRk *rk,
         GreensFunction *gExt,
         SurfaceParameterizationCache *cache,
         //         omp_lock_t *cacheLock,
         Options &options,
         Real gap,
         int basis)
{
  TriangleRWGMap tt;
  TriangleRWGMap ss;

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
  GreensFunction *gInt = t->g == s->g?t->g:NULL;
    
  // choose values for nt, ns, which determine the number of points in the interpolation
  int nt;
  int ns;

  Real k2 = Square(options.gapNormalization / gap);
  Real dt = options.interpolationDensity * pow(t->totalArea * k2, options.interpolationPowerLaw);
  Real ds = options.interpolationDensity * pow(s->totalArea * k2, options.interpolationPowerLaw);
  
  nt = 3 + (int)dt;
  ns = 3 + (int)ds;

  Point *pt;
  int npt;
  Point *ps;
  int nps;

  if(t->rwgs < options.minParameterizableSegmentSize ||
     (Real)SurfaceParameterization::hexagonalPoints(nt) > options.parameterizationThresh * (Real)t->rwgs) {
    getSurfacePoints(t,&pt,&npt);
  } else {
    SurfaceParameterization *ut;
    if(!t->bParameterizable) {
      cout << "RWGTree not parameterizable\n";
      abort();
    }
    if(cache) {
      //omp_set_lock(cacheLock);
      SurfaceParameterizationCache::iterator uti = cache->find(t);
      if(uti != cache->end()) {
        ut = uti->second;
      } else {
        ut = t->parameterize();
        (*cache)[t] = ut;
      }
      //omp_unset_lock(cacheLock);
    } else {
      ut = t->parameterize();
    }
    getSurfacePoints(nt,*ut,&pt,&npt);
    if(!cache) delete ut;    
  }

  if(s->rwgs < options.minParameterizableSegmentSize ||
     (Real)SurfaceParameterization::hexagonalPoints(ns) > options.parameterizationThresh * (Real)s->rwgs) {
    getSurfacePoints(s,&ps,&nps);
  } else {
    SurfaceParameterization *us;
    if(!s->bParameterizable) {
      cout << "RWGTree not parameterizable\n";
      abort();
    }
    if(cache) {
      //omp_set_lock(cacheLock);
      SurfaceParameterizationCache::iterator usi = cache->find(s);
      if(usi != cache->end()) {
        us = usi->second;
      } else {
        us = s->parameterize();
        (*cache)[s] = us;
      }
      //omp_unset_lock(cacheLock);
    } else {
      us = s->parameterize();
    }
    getSurfacePoints(ns,*us,&ps,&nps);
    if(!cache) delete us;
  }

  ComplexDenseMatrix G(2*npt,2*nps);
  for(int i=0; i<npt; i++) {
    int ii = i<<1;
    for(int j=0; j<nps; j++) {
      int jj = j<<1;
      G(ii,jj) = gExt->scalar(pt+i,ps+j);
      G(ii,jj+1) = Complex(0.,0.);
      G(ii+1,jj) = Complex(0.,0.);
      G(ii+1,jj+1) = gInt?gInt->scalar(pt+i,ps+j):0.0;
    }
  }

  vector<int> r;
  vector<int> c;
  aca(G,r,c,options.acaTolerance);

  int m = t->rwgs;
  int n = s->rwgs;
  int rows = (int) r.size();
  int cols = (int) c.size();

  ComplexVector3Matrix Av00(2*m,cols,0);
  ComplexVector3Matrix Av01(2*m,cols,0);
  ComplexVector3Matrix Av10(2*m,cols,0);
  ComplexVector3Matrix Av11(2*m,cols,0);

  ComplexVector3Matrix Bv00(rows,2*n,0);
  ComplexVector3Matrix Bv01(rows,2*n,0);
  ComplexVector3Matrix Bv10(rows,2*n,0);
  ComplexVector3Matrix Bv11(rows,2*n,0);
  for(TriangleRWGMap::iterator i = tt.begin(); i != tt.end(); ++i) {
    vector<RWGIndexPair> &t = i->second;
    for(int j=0; j<cols; j++) {
      int jj = c[j];
      int sindex = jj>>1;
      GreensFunction *g = jj%2==0?gExt:gInt;
      LK1Integrand integrand(ps[sindex], g, t, basis);
      TriangleQuadratureGauss10 quad;
      HalfRWG &th = *t[0].first;
      quad.integrate(integrand,*th.p0,th.v01,th.v02);
      for(unsigned int k=0; k<t.size(); k++) {
        int ii = t[k].second<<1;
        Real c = t[k].first->e->length;
        if(!t[k].first->bPositive) {
          c = -c;
        }
        ComplexVector3 Gphi = c * integrand.Gphi[k];
        Av00(ii,j) += g->Zeps * Gphi;
        Av01(ii,j) += Gphi;
        Av10(ii+1,j) -= Gphi;
        Av11(ii+1,j) += g->Zmu * Gphi;
      }
    }
  }

  for(TriangleRWGMap::iterator i = ss.begin(); i != ss.end(); ++i) {
    vector<RWGIndexPair> &s = i->second;
    for(int j=0; j<rows; j++) {
      int jj = r[j];
      int tindex = jj>>1;
      TriangleQuadratureGauss10 quad;
      LK2Integrand integrand(pt[tindex], jj%2==0?gExt:gInt, s);
      HalfRWG &sh = *s[0].first;
      quad.integrate(integrand,*sh.p0,sh.v01,sh.v02);
      for(unsigned int k = 0; k<s.size(); k++) {
        int ii = s[k].second<<1;
        Real c = ONEOVER16PI * s[k].first->e->length;
        if(!s[k].first->bPositive) {
          c = -c;
        }
        ComplexVector3 L2 = c * integrand.L[k];
        ComplexVector3 K2 = c * integrand.K[k];
        Bv00(j,ii) += L2;
        Bv01(j,ii+1) += K2;
        Bv10(j,ii) += K2;
        Bv11(j,ii+1) += L2;
      }
    }
  }

  delete [] pt;
  delete [] ps;

  ComplexDenseMatrix S = G(r,c);

  ComplexDenseMatrix A00(2*m,3*cols);
  ComplexDenseMatrix A01(2*m,3*cols);
  ComplexDenseMatrix A10(2*m,3*cols);
  ComplexDenseMatrix A11(2*m,3*cols);

  for(int i=0; i<m; i++) {
    int ii = 2 * i;
    int ii1 = ii + 1;
    for(int j=0; j<cols; j++) {
      int jj = 3 * j;
      ComplexVector3 v;

      v = Av00(ii,j);
      A00(ii,jj) = v.x;
      A00(ii,jj+1) = v.y;
      A00(ii,jj+2) = v.z;
      v = Av00(ii1,j);
      A00(ii1,jj) = v.x;
      A00(ii1,jj+1) = v.y;
      A00(ii1,jj+2) = v.z;

      v = Av01(ii,j);
      A01(ii,jj) = v.x;
      A01(ii,jj+1) = v.y;
      A01(ii,jj+2) = v.z;
      v = Av01(ii1,j);
      A01(ii1,jj) = v.x;
      A01(ii1,jj+1) = v.y;
      A01(ii1,jj+2) = v.z;

      v = Av10(ii,j);
      A10(ii,jj) = v.x;
      A10(ii,jj+1) = v.y;
      A10(ii,jj+2) = v.z;
      v = Av10(ii1,j);
      A10(ii1,jj) = v.x;
      A10(ii1,jj+1) = v.y;
      A10(ii1,jj+2) = v.z;

      v = Av11(ii,j);
      A11(ii,jj) = v.x;
      A11(ii,jj+1) = v.y;
      A11(ii,jj+2) = v.z;
      v = Av11(ii1,j);
      A11(ii1,jj) = v.x;
      A11(ii1,jj+1) = v.y;
      A11(ii1,jj+2) = v.z;
    }
  }

  S.LU();
  solveLU(S,Bv00);
  solveLU(S,Bv01);
  solveLU(S,Bv10);
  solveLU(S,Bv11);

  ComplexDenseMatrix B00(3*rows,2*n);
  ComplexDenseMatrix B01(3*rows,2*n);
  ComplexDenseMatrix B10(3*rows,2*n);
  ComplexDenseMatrix B11(3*rows,2*n);
  for(int i=0; i<rows; i++) {
    int ii = 3 * i;
    for(int j=0; j<n; j++) {
      int jj = 2 * j;
      int jj1 = jj + 1;
      ComplexVector3 v;

      v = Bv00(i,jj);
      B00(ii,jj) = v.x;
      B00(ii+1,jj) = v.y;
      B00(ii+2,jj) = v.z;
      v = Bv00(i,jj1);
      B00(ii,jj1) = v.x;
      B00(ii+1,jj1) = v.y;
      B00(ii+2,jj1) = v.z;

      v = Bv01(i,jj);
      B01(ii,jj) = v.x;
      B01(ii+1,jj) = v.y;
      B01(ii+2,jj) = v.z;
      v = Bv01(i,jj1);
      B01(ii,jj1) = v.x;
      B01(ii+1,jj1) = v.y;
      B01(ii+2,jj1) = v.z;

      v = Bv10(i,jj);
      B10(ii,jj) = v.x;
      B10(ii+1,jj) = v.y;
      B10(ii+2,jj) = v.z;
      v = Bv10(i,jj1);
      B10(ii,jj1) = v.x;
      B10(ii+1,jj1) = v.y;
      B10(ii+2,jj1) = v.z;

      v = Bv11(i,jj);
      B11(ii,jj) = v.x;
      B11(ii+1,jj) = v.y;
      B11(ii+2,jj) = v.z;
      v = Bv11(i,jj1);
      B11(ii,jj1) = v.x;
      B11(ii+1,jj1) = v.y;
      B11(ii+2,jj1) = v.z;
    }
  }         

  ComplexRk rk00(A00,B00,options.rkEps);
  ComplexRk rk01(A01,B01,options.rkEps);
  ComplexRk rk10(A10,B10,options.rkEps);
  ComplexRk rk11(A11,B11,options.rkEps);

  ComplexRk rk0 = rk00 + rk01;
  ComplexRk rk1 = rk10 + rk11;

  *rk = rk0 + rk1;
  cout << rk->m << " " << rk->n << " " << rk->k << " " << rk00.k << " " << rk01.k << " " << rk10.k << " " << rk11.k << "\n";
}

Real projectedDistance(RWGTree *t, RWGTree *s, Vector3 &vst)
{
  Real maxd = 0.;
  for(int i=0;i<s->rwgs;i++) {
    RWG *rwg = s->rwg[i];

    Triangle *t;
    Vector3 v;
    Real d;

    t = rwg->pos.t;

    v = *(t->p1) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p2) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p3) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    t = rwg->neg.t;

    v = *(t->p1) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p2) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p3) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;
  }
  return maxd;
}

Real admissable(RWGTree *t, RWGTree *s, Options &options)
{
  return -1.0;
  if(t==s) return -1.0;
  if(t->rwgs < options.minAdmissableSegmentSize || s->rwgs < options.minAdmissableSegmentSize) return false;
  if(!t->bParameterizable || !s->bParameterizable) return false;
  Vector3 v = t->centroid - s->centroid;
  Real d = sqrt(norm(v));
  v *= 1. / d;
  Real d1 = projectedDistance(t,s,v);
  v *= -1.;
  Real d2 = projectedDistance(s,t,v);
  Real gap = d - (d1 + d2);
  Real r = max(t->radius,s->radius);
  if(gap > r * options.admissableThresh)
    return gap;
  else 
    return -1.;
}


ComplexDenseMatrix *Zts(RWGTree *t, 
                        RWGTree *s,
                        GreensFunction *g,
                        Options &options,
                        int basis)
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
      LKTriangle(g,i->first,j->first,i->second,j->second,basis,*dense,options);
      if(g->gLong) {
        LKTriangle(t->g->gLong,i->first,j->first,i->second,j->second,basis,*dense,options);
        //LKTriangleLongitudinal(t->g->gLong,i->first,j->first,i->second,j->second,*dense,options);
      }
    }
  }
  return dense;
}


ComplexDenseMatrix *ZtsNcross(RWGTree *t, 
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
      LKTriangle(g,i->first,j->first,i->second,j->second,nCrossRWGBasis,*dense,options);
    }
  }
  return dense;
}

ComplexDenseMatrix *ZtsDot(RWGTree *t, 
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
    for(TiangleRWGMap::iterator j = ss.begin(); j != ss.end(); ++j) {
      // XXX imp?
      //LKTriangleDot(g,i->first,j->first,i->second,j->second,*dense,options);
    }
  }
  return dense;
}


ComplexMatrix * Zh(RWGTree *t,
                   RWGTree *s,
                   GreensFunction *gExt,
                   ImpedanceJobVector &jobs,
                   Options &options)
{
  Real gap = admissable(t,s,options);
  if(gap > 0.) {
    ComplexRk *rk = new ComplexRk;
    ImpedanceJob job(t,s,rk,gap);
    jobs.push_back(job);
    return rk;
  }
  if(!t->sub && !s->sub) {
    ComplexDenseMatrix *dense = new ComplexDenseMatrix;
    ImpedanceJob job(t,s,dense,gap);
    jobs.push_back(job);
    return dense;
  }
  if(!t->sub) {
    ComplexHMatrix *Hp = new ComplexHMatrix(1,s->subs);
    ComplexHMatrix &H = *Hp;
    H.r[0] = 2*t->rwgs;
    for(int j=0; j<s->subs; j++)
      H.c[j] = 2*s->sub[j]->rwgs;
    for(int j=0; j<s->subs; j++) {
      ComplexMatrix *M = Zh(t,s->sub[j],gExt,jobs,options);
      H.set(0,j,M);
    }
    return Hp;
  }
  if(!s->sub) {
    ComplexHMatrix *Hp = new ComplexHMatrix(t->subs,1);
    ComplexHMatrix &H = *Hp;
    for(int i=0; i<t->subs; i++)
      H.r[i] = 2*t->sub[i]->rwgs;
    H.c[0] = 2*s->rwgs;
    for(int i=0; i<t->subs; i++) {
      ComplexMatrix *M = Zh(t->sub[i],s,gExt,jobs,options);
      H.set(i,0,M);
    }
    return Hp;
  }

  ComplexHMatrix *Hp = new ComplexHMatrix(t->subs,s->subs);
  ComplexHMatrix &H = *Hp;
  for(int i=0; i<t->subs; i++)
    H.r[i] = 2*t->sub[i]->rwgs;
  for(int j=0; j<s->subs; j++)
    H.c[j] = 2*s->sub[j]->rwgs;
  for(int i=0; i<t->subs; i++) {
    for(int j=0; j<s->subs; j++) {
      ComplexMatrix *M = Zh(t->sub[i],s->sub[j],gExt,jobs,options);
      H.set(i,j,M);
    }
  }
  return Hp;
}

ComplexMatrix * Z(RWGTree *t,
                  RWGTree *s,
                  GreensFunction *gExt,
                  Progress &progress,
                  Options &options)
{
  ImpedanceJobVector jobs;
  ComplexMatrix *Zhp = Zh(t,s,gExt,jobs,options);

  SurfaceParameterizationCache cache;
  //omp_lock_t cacheLock;
  //omp_lock_t progressLock;
  //omp_init_lock(&cacheLock);
  //omp_init_lock(&progressLock);

  //int chunk = 1;
  int n = jobs.size();
  int i;
  #pragma omp parallel for default(shared) private(i) schedule(static,chunk) 
  for(i=0; i<n; i++) {
    ImpedanceJob &job = jobs[i];
    //job.perform(gExt,&cache,&cacheLock,progress,&progressLock,options);
    job.perform(gExt,&cache,progress,options);
  }

  return Zhp;
}

ImpedanceJob :: ImpedanceJob(RWGTree *t, RWGTree *s, ComplexMatrix *M, Real gap)
{
  this->t = t;
  this->s = s;
  this->M = M;
  this->gap = gap;
}

void ImpedanceJob :: perform(GreensFunction *g,
                             SurfaceParameterizationCache *cache,
                             //omp_lock_t *cacheLock,
                             Progress &progress,
                             //omp_lock_t *progressLock,
                             Options &options)
{
  if(M->type == ComplexMatrix::dense) {
    ComplexDenseMatrix *dense = (ComplexDenseMatrix*)M;
    //Zdense(t,s,dense,gExt,options);
    Zdense(t,s,dense,g,options);
    //omp_set_lock(progressLock);
    progress.finishedMatrix(*dense);
    //omp_unset_lock(progressLock);
  } else {
    abort();
    ComplexRk *rk = (ComplexRk*)M;
    //Zrk(t,s,rk,g,cache,cacheLock,options,gap,basis);
    Zrk(t,s,rk,g,cache,options,gap,RWGBasis);
    // XXX add longitudinal?
    if(t->g->gLong) {
      Zrk(t,s,rk,g,cache,options,gap,nCrossRWGBasis);
    }

    //omp_set_lock(progressLock);
    progress.finishedMatrix(*rk);
    //omp_unset_lock(progressLock);
  }
}

  
