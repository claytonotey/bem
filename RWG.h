#ifndef RWG_H
#define RWG_H

#include "Vector3.h"
#include "ZRotation.h"
#include "Edge.h"
#include "Triangle.h"

class HalfRWG {
 public:
  HalfRWG(Edge *e, Triangle *t, bool bPositive);
  Real Ai;
  Real normv01i, normv20i, normv12i;
  Vector3 *p0;
  Vector3 *p1;
  Vector3 *p2;
  Vector3 v01, v12, v20, v02;
  Vector3 n;
  Vector3 m01, m12, m20;
  Edge *e;
  Triangle *t;
  bool bPositive;
};

class RWG {
 public:
  RWG(Edge *e);
  Edge *e;
  HalfRWG pos;
  HalfRWG neg;
};

typedef pair<HalfRWG*,int> RWGIndexPair;
typedef UnorderedMap< Triangle*,vector<RWGIndexPair> > TriangleRWGMap;

#endif
