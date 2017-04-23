#ifndef FIELD_H
#define FIELD_H

#include "RWGTree.h"
#include "GreensFunction.h"
#include "Source.h"
#include "Options.h"
#include "Matrix.h"
#include "Progress.h"
#include "DenseMatrix.h"

ComplexDenseMatrix *Bdense(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options);
ComplexMatrix * BH(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options);

#endif
