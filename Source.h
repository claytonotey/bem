#ifndef SOURCE_H
#define SOURCE_H

#include "Progress.h"
#include "defs.h"
#include "Vector3.h"
#include "GreensFunction.h"
#include "RWG.h"
#include "ComplexPair.h"
#include "Options.h"
#include "FunctionArgs.h"
#include "TriangleQuadrature.h"

#include <vector>
#include <utility>
using namespace std;

class Source {
 public:
  Source() {};
  virtual ~Source() {};
  virtual void EHTriangle(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options)=0;
  virtual ComplexVector3Pair getField(const Point &p)=0;
  virtual bool isInRegion(GreensFunction *g)=0;
};

class DipolePointSource : public Source {
 public:
  GreensFunction *g;
  Point p;
  Vector3 j;

  DipolePointSource(GreensFunction *g, const Point &q, Vector3 &j);
  void EHTriangle(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options);
  void EHTriangleLongitudinal(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options);
  bool isInRegion(GreensFunction *g);
  ComplexVector3Pair getField(const Point &p);
};

class PlaneWaveSource : public Source {
 public:
  GreensFunction *g;
  Complex k;
  Vector3 n;
  ComplexVector3 E0;
  ComplexVector3 H0;

  PlaneWaveSource(GreensFunction *g, const Vector3 &n, const ComplexVector3 &E0);
  void EHTriangle(GreensFunction *g, Triangle *trit, vector<RWGIndexPair> &t, ComplexDenseMatrix &B, Options &options);
  bool isInRegion(GreensFunction *g);
  ComplexVector3Pair getField(const Point &p);
};

pair<ComplexVector3,ComplexVector3> EH(GreensFunction *green, Vector3 &j, Vector3 *rt, Vector3 *rs);

typedef UnorderedSet<Source*> SourceSet;

void EHNonSmoothIntegral(HalfRWG &sh, const Vector3 &rt, Real *ISn3h, Real *ISn1, Real *IS1, Vector3 *ILn1, Vector3 *IL1, Vector3 *IL3);
ComplexDenseMatrix *Bdense(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options);
ComplexMatrix * BH(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options);

Real extinctionCrossSection(RWGTree *t, const ComplexDenseMatrix &JM, GreensFunction *g, PlaneWaveSource *plane);

ComplexVector3Pair getField(const SourceSet &sources, RWGTree *t, const ComplexDenseMatrix &JM, GreensFunction *g, const Point &p, bool bLong, Options &options);
#endif
