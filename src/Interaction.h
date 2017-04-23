#ifndef INTERACTION_H
#define INTERACTION_H

#include "Vector3.h"
#include "ComplexVector3.h"
#include "Array.h"
#include "DenseMatrix.h"
#include "MatrixOperations.h"

#include <vector>
#include <utility>
using namespace std;

Key getKey(Key t, Key s);

class TrianglePointInteraction {
 public:
  TrianglePointInteraction(Vector3& p0_, Vector3& v01_, Vector3& v02, int n, int useTotal);
  ~TrianglePointInteraction();
  void getScalar(Complex &z);
  void getGradient(ComplexVector3 &v);
  void setScalar(Complex &z);
  void setGradient(ComplexVector3 &v);
  void start();
  void end();
  void next();
  
  Vector3 p0,v01,v02;
  int n;
  bool bSet;
  int i;
  int useCount;
  int useTotal;
  Array<Complex> scalar;
  Array<ComplexVector3> gradient;
};

class TriangleTriangleInteraction {
 public:
  TriangleTriangleInteraction(Vector3& p0t_, Vector3& v01t_, Vector3& v02t_, Vector3& p0s_, Vector3& v01s_, Vector3& v02s_, int nt, int ns, int useTotal);
  ~TriangleTriangleInteraction();
  TriangleTriangleInteraction();
  void getScalar(Complex &z);
  void getGradient(ComplexVector3 &v);
  void setScalar(Complex &z);
  void setGradient(ComplexVector3 &v);

  void start();
  void end();
  void next();
  
  Vector3 p0t, v01t, v02t;
  Vector3 p0s, v01s, v02s;
  int nt,ns;
  bool bSet;
  int i;
  int useCount;
  int useTotal;
  Array<Complex> scalar;
  Array<ComplexVector3> gradient;
};


typedef UnorderedMap< Key, TriangleTriangleInteraction* > TriangleTriangleInteractionCache;
typedef UnorderedMap< Key, char > TriangleTriangleInteractionCount;
typedef UnorderedMap< Key, TrianglePointInteraction* > TrianglePointInteractionCache;
typedef UnorderedMap< Key, char > TrianglePointInteractionCount;

#endif
