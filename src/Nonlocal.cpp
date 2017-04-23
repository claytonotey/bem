#include "Nonlocal.h"
#include "Source.h"

class LKInnerIntegrandLongitudinal : public TriangleIntegrand {
public:
  Complex G;
  GreensFunction *g;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  LKInnerIntegrandLongitudinal(GreensFunction *g, vector<RWGIndexPair> &s) : g(g), s(s) {}
  
  void setTestPoint(Vector3 &v) {
    rt = v;
    G = 0.0;
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs = g->scalar(&rt,&v);
    gs *= w;
    G += gs;
  }
};

class LKFarOuterIntegrandLongitudinal : public TriangleIntegrand {
public:
  GreensFunction *g;
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  Complex *L;
  LKInnerIntegrandLongitudinal inner;
  TriangleQuadrature &innerquad;

  LKFarOuterIntegrandLongitudinal(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, Complex *L, TriangleQuadrature &innerquad) 
    : g(g),t(t),s(s),L(L),inner(g,s),innerquad(innerquad) {}

  void weighPoint(Real w, Vector3 &v) {
    inner.setTestPoint(v);
    HalfRWG &sh = *s[0].first;
    innerquad.integrate(inner,*sh.p0,sh.v01,sh.v02);
    int ij = 0;

    for(unsigned int i=0; i<t.size(); i++) {
      for(unsigned int j=0; j<s.size(); j++) {
        L[ij] += w * -inner.G;
        ij++;
      }
    }
  }  
};

class LKSmoothInnerIntegrandLongitudinal : public TriangleIntegrand {
public:
  ComplexVector3 Gphi[3];
  GreensFunction *g;
  Complex G;
  vector<RWGIndexPair> &s;
  Vector3 rt;

  LKSmoothInnerIntegrandLongitudinal(GreensFunction *g, vector<RWGIndexPair> &s) : g(g), s(s) {}
  
  void setTestPoint(Vector3 &v) {
    rt = v;
    G = 0.0;
  }

  void weighPoint(Real w, Vector3 &v) { 
    Complex gs = g->scalar_smooth(&rt,&v);
    gs *= w;
    G += gs;
  }
};

class LKCloseOuterIntegrandLongitudinal : public TriangleIntegrand {
public:
  GreensFunction *g;
  vector<RWGIndexPair> &t;
  vector<RWGIndexPair> &s;
  Complex *L;
  LKSmoothInnerIntegrandLongitudinal inner_smooth;
  TriangleQuadratureGauss91 innerquad;
  Vector3 vts[9];
  bool bInPlane;

  LKCloseOuterIntegrandLongitudinal(GreensFunction *g, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, Complex *L) 
    : g(g),t(t),s(s),L(L),inner_smooth(g,s) {
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
    Complex ISn11 = Ai2 * (ISn1 + k2 * IS1);
    Complex L2 = 2.0 * ISn11 + inner_smooth.G;
    int ij = 0;
    for(unsigned int i=0; i<t.size(); i++) {
      for(unsigned int j=0; j<s.size(); j++) {
        L[ij] += w * -L2;
        ij++;
      }
    }
  }
};

/* Only the grad g(k,r) term survives */
void LKTriangleLongitudinal(GreensFunction *g, Triangle *trit, Triangle *tris, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, ComplexDenseMatrix &Z, Options &options)
{
  Complex L[9];
  memset(L,0,9*sizeof(Complex));

  HalfRWG &th = *(t[0].first);

  if(trit == tris) {
#ifdef DEBUG
    printf("Same\n");
#endif
    TriangleQuadratureGauss105 outerquad;
    LKCloseOuterIntegrandLongitudinal outerintegrand(g,t,s,L);
    outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
  } else {
    int shared = sharedpoints(trit,tris);
    if(shared > 0) {
#ifdef DEBUG
      printf("Neighbors\n");
#endif
      TriangleQuadratureGauss105 outerquad;
      LKCloseOuterIntegrandLongitudinal outerintegrand(g,t,s,L);
      outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
    } else {
#ifdef DEBUG
      printf("Far\n");
#endif
      Point pt = center(trit);
      Point ps = center(tris);
      Vector3 v = pt - ps;
      Real d = norm2(v);
      Real temax = max(max(trit->e12->length,trit->e23->length),trit->e31->length);
      Real semax = max(max(tris->e12->length,tris->e23->length),tris->e31->length);
      Real emax2 = Square(max(temax,semax));
      if(d * options.farThresh6 > emax2) {
        TriangleQuadratureGauss6 innerquad;
        TriangleQuadratureGauss6 outerquad;
        LKFarOuterIntegrandLongitudinal outerintegrand(g,t,s,L,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh10 > emax2) {
        TriangleQuadratureGauss10 innerquad;
        TriangleQuadratureGauss10 outerquad;
        LKFarOuterIntegrandLongitudinal outerintegrand(g,t,s,L,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh15 > emax2) {
        TriangleQuadratureGauss15 innerquad;
        TriangleQuadratureGauss15 outerquad;
        LKFarOuterIntegrandLongitudinal outerintegrand(g,t,s,L,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else if(d * options.farThresh21 > emax2) {
        TriangleQuadratureGauss21 innerquad;
        TriangleQuadratureGauss21 outerquad;
        LKFarOuterIntegrandLongitudinal outerintegrand(g,t,s,L,innerquad);
        outerquad.integrate(outerintegrand,*th.p0,th.v01,th.v02);
      } else {
        TriangleQuadratureGauss36 innerquad;
        TriangleQuadratureGauss36 outerquad;
        LKFarOuterIntegrandLongitudinal outerintegrand(g,t,s,L,innerquad);
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
      if(isnan(Lij)) {
        cout << *(th.t) << " " << *(sh.t);
        cout << Lij << "\n";
        abort();
      }
      //cout << Z(ii,jj) << " " << g->Zeps * Lij << "\n";
      //cout << sqrt(norm2(Z(ii,jj) / (g->Zeps * Lij))) << "\n";
      //Z(ii,jj) -= g->Zeps * Lij;
      //Z(ii+1,jj+1) += g->eps0 / g->mu * g->Zeps * Lij;
      ij++;
    
    }
  }
}

