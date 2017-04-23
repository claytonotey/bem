#ifndef NONLOCAL_H
#define NONLOCAL_H

#include "GreensFunction.h"
#include "Vector3.h"
#include "ComplexVector3.h"
#include "FunctionArgs.h"
#include "RWG.h"
#include "Options.h"
#include "DenseMatrix.h"
#include "TriangleQuadrature.h"

void LKTriangleLongitudinal(GreensFunction *g, Triangle *trit, Triangle *tris, vector<RWGIndexPair> &t, vector<RWGIndexPair> &s, ComplexDenseMatrix &Z, Options &options);

#endif
