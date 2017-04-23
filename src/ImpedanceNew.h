#ifndef IMPEDANCENEW_H
#define IMPEDANCENEW_H


enum RWGBasisType {
  RWGBasis,
  nCrossRWGBasis
};


#include "GreensFunction.h"
#include "Vector3.h"
#include "ComplexVector3.h"
#include "FunctionArgs.h"
#include "RWG.h"
#include "ComplexPair.h"
#include "RWGTree.h"
#include "SurfaceParameterization.h"
#include "Options.h"
#include "DenseMatrix.h"
#include "Rk.h"
#include "HMatrix.h"
#include "Progress.h"
#include "Interaction.h"


#include <vector>
using namespace std;


typedef UnorderedMap< RWGTree*, SurfaceParameterization* > SurfaceParameterizationCache;
typedef DenseMatrix<Real,ComplexVector3> ComplexVector3Matrix;

class ImpedanceJob {
public:
  ImpedanceJob(RWGTree *t, RWGTree *s, ComplexMatrix *M, Real gap = 0.);
  void perform(GreensFunction *gExt,
               SurfaceParameterizationCache *cache,
               //omp_lock_t *cacheLock,
               Progress &progress,
               //omp_lock_t *progressLock,
               Options &options);

  RWGTree *t;
  RWGTree *s;
  ComplexMatrix *M;
  Real gap;
};

typedef vector<ImpedanceJob> ImpedanceJobVector;

void Zdense(RWGTree *t, RWGTree *s, ComplexDenseMatrix *dense, GreensFunction *gExt, Options &options);
ComplexMatrix * Z(RWGTree *t, RWGTree *s, GreensFunction *gExt, Progress &progress, Options &options);

#endif
