#ifndef IMPEDANCE_H
#define IMPEDANCE_H

#include "GreensFunction.h"
#include "Vector3.h"
#include "ComplexVector3.h"
#include "FunctionArgs.h"
#include "RWG.h"
#include "ComplexPair.h"
#include "RWGTree.h"
#include "Options.h"
#include "DenseMatrix.h"
#include "Progress.h"

#include <vector>
using namespace std;

typedef DenseMatrix<Real,ComplexVector3> ComplexVector3Matrix;

void Zdense(RWGTree *t, RWGTree *s, ComplexDenseMatrix *dense, GreensFunction *gExt, Options &options);
ComplexDenseMatrix *Zts(RWGTree *t, RWGTree *s, GreensFunction *gExt, Options &options);

#endif
