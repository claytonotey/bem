#ifndef RWGTREE_H
#define RWGTREE_H

#include "ShapeMesh.h"
#include "RWG.h"
#include "Dielectric.h"
#include "GreensFunction.h"

#include <vector>
#include <utility>
#include <list>
using namespace std;

class RWGTree;

struct ltrwg {
  bool operator()(RWG *r1, RWG *r2) const {
    return r1->e->index < r2->e->index;
  }
};

typedef set<RWG*, ltrwg> RWGSet;
typedef UnorderedMap<Segment*, Real> SegmentRealMap;
typedef pair<ShapeMesh*, Dielectric*> ShapeMeshDielectricPair;
typedef list<ShapeMeshDielectricPair> ShapeMeshDielectricPairList;

class RWGTree {
 public:
  RWGTree(ShapeMeshDielectricPairList &meshDiel,
          PointEdgeVectorMap &edgesAtPoint,
          int minSize,
          Real minRadius);
  void output(ostream& os, int *k);
  friend ostream& operator<<(ostream& os, RWGTree *t);
  RWGTree(RWGSet &rwgSet,
          Segment *s,
          PointEdgeVectorMap &edgesAtPoint,
          Dielectric *diel,
          GreensFunction *g,
          int minSize,
          Real minRadius,
          bool bTerminal);
 protected:
  RWGTree(int rwgs,
          PointEdgeVectorMap &edgesAtPoint,
          Dielectric *diel,
          GreensFunction *g);

  void partition(int minSize,
                 Real minRadius);
  void partition(Real *angles,
                 int minSize = 0,
                 Real minRadius = 0.);
  void unshuffle();
 public:
  bool bParameterizable;
  RWG **rwg;
  int rwgs;
  vector<int> index;
  RWGTree **sub;
  int subs;
  EdgeSet boundary;
  TriangleSet triangle;
  PointEdgeVectorMap &edgesAtPoint;
  Dielectric *diel;
  GreensFunction *g;
  Point centroid;
  Vector3 n;
  Real radius;
  Real totalArea;
};

#endif
