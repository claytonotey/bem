#ifndef HEATTRANSFER_H
#define HEATTRANSFER_H

#include "ComplexVector3.h"
#include "Vector3.h"
#include "MathUtils.h"
#include "Impedance.h"
#include "GreensFunction.h"
#include "RWG.h"
#include "ShapeParser.h"
#include "ShapeMesh.h"
#include "Options.h"
#include "MatrixOperations.h"
#include "SI.h"
#include "FileUtils.h"
#include <pthread.h>

class Point2 {
public:
  double x,y;
};

inline istream& operator >> (istream& is, Point2 &ps)
{
  return is >> ps.x >> ps.y;
}

inline bool refineEdgeFuncBody(Edge *e, void *args)
{
  return false;
}

Real CalcHeatTransfer(Real a, Real eV, Real T1, Real T2, RWGTree *top, RWGTree *s, RWGTree *t, Options &options);

void parseShapes(RWGTree *top, RWGTree *s, RWGTree *t, char *body1x3d, char *body2x3d, const Vector3 &gap, Options &options);

#endif
