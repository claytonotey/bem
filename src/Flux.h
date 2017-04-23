#ifndef FLUX_H
#define FLUX_H

#include "Util.h"
#include "VolumeMesh.h"
#include "ShapeMesh.h"
#include "Matrix.h"
#include "RWGTree.h"
#include "Progress.h"
#include "Tetrahedron.h"
#include "FunctionArgs.h"

class TetQuadratureFunctionFactory {
 public:
  virtual ~TetQuadratureFunctionFactory() {}
  virtual TetQuadratureFunction<Real,RealValuedFunction>::Type getQuadFunc(Tetrahedron *tet)=0;
};

void enumFluxes(RWGTree *top, RWGTree *t);
Real flux(RWGTree *top, RWGTree *t, const ComplexDenseMatrix &X);
Real totalFlux(VolumeMesh &vMesh, ComplexHMatrix *ZLU, GreensFunction *greenExt, RWGTree *rwgTree, RWGTree *s, RWGTree *t, TetQuadratureFunctionFactory &tetquadfuncFactory, Progress &progress);

#endif
