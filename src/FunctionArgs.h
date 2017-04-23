#ifndef FUNCTIONARGS_H
#define FUNCTIONARGS_H

#include "GreensFunction.h"
#include "Vector3.h"
#include "ComplexPair.h"
#include "ComplexVector3Pair.h"
#include "RWG.h"
#include "Interaction.h"
#include "RWGTree.h"
#include "Triangle.h"
#include "Tetrahedron.h"

class FunctionArgs;
typedef void (*RealValuedFunction)(Vector3 *r, Real *x, FunctionArgs *args);
typedef void (*VectorValuedFunction)(Vector3 *r, Vector3 *x, FunctionArgs *args);
typedef void (*ComplexVectorValuedFunction)(Vector3 *r, ComplexVector3 *x, FunctionArgs *args);
typedef void (*ComplexVectorPairValuedFunction)(Vector3 *r, ComplexVector3Pair *x, FunctionArgs *args);
typedef void (*ComplexValuedFunction)(Vector3 *r, Complex *x, FunctionArgs *args);
typedef void (*ComplexPairValuedFunction)(Vector3 *r, ComplexPair *x, FunctionArgs *args);

template<class T, class F>
struct QuadratureFunction {
  typedef T (*Type)(F f,
                    FunctionArgs *args,
                    Vector3 &v1, 
                    Vector3 &v2, 
                    Vector3 *v0in,
                    void *quadargs,
                    Real reltol,
                    Real abstol);

};

template<class T, class F>
struct TetQuadratureFunction {
  typedef T (*Type)(F f,
                    FunctionArgs *args,
                    Vector3 *v0in,
                    Vector3 &v1, 
                    Vector3 &v2, 
                    Vector3 &v3);

}; 
                    
class FunctionArgs {
public:
  RWGTree *rwgTree;
  RWGTree *t;
  RWGTree *s;
  Tetrahedron *tet;
  Vector3 *jSource;
  ComplexHMatrix *ZLU;
  GreensFunction *green;
  GreensFunction *greenExt;
};


#endif
