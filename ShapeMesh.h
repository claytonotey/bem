#ifndef SHAPEMESH_H
#define SHAPEMESH_H

#include "Shape.h"
#include "defs.h"
#include "Vector3.h"
#include "Edge.h"
#include "Triangle.h"
#include "Array.h"
#include "VolumeMesh.h"

#include <vector>
#include <set>
using namespace std;

struct ltedge {
  bool operator()(Edge *e1, Edge *e2) const {
    return e1->index < e2->index;
  }
};

struct lttriangle {
  bool operator()(Triangle *t1, Triangle *t2) const {
    return t1->index < t2->index;
  }
};

struct ltpoint {
  bool operator()(Point *p1, Point *p2) const {
    if(p1->x == p2->x) {
      if(p1->y == p2->y) {
        return p1->z < p2->z;
      }
      return p1->y < p2->y;
    }
    return p1->x < p2->x;
  }
};

typedef set<Point*,ltpoint> PointSet;
typedef set<Edge*,ltedge> EdgeSet;
typedef set<Triangle*,lttriangle> TriangleSet;
typedef UnorderedMap<Point*, vector<Edge*> > PointEdgeVectorMap;
typedef UnorderedMap<Point*, Real > PointRealMap;
typedef UnorderedMap<Edge*, Edge* > EdgeEdgeMap;

struct Cut {
  Point *p1;
  Point *p2;
  vector<Edge*> e1;
  vector<Edge*> e2;
  vector<Point*> point;
  vector<Edge*> edge;
};

struct ltcut {
  bool operator()(Cut *c1, Cut *c2) const {
    return (c1->edge.size()>0 && c2->edge.size()>0) &&
      ((*(c1->edge.begin()))->index < (*(c2->edge.begin()))->index);
  }
};

typedef set<Cut*,ltcut> CutSet;

class Segment {
 public:
  EdgeSet edge;
  EdgeSet boundary;
  TriangleSet triangle;
};

struct ltsegment {
  bool operator()(Segment *s1, Segment *s2) const {
    return (s1->edge.size()>0 && s2->edge.size()>0) &&
      ((*(s1->edge.begin()))->index < (*(s2->edge.begin()))->index);
  }
};

typedef set<Segment*,ltsegment> SegmentSet;

class ShapeMesh {
 public:
  ShapeMesh(Shape &s, Real sharpThresh, Real concaveThresh, Real relativeCost, bool (*refineEdgeFunc)(Edge *, void *), void *args);
  ~ShapeMesh();

  void generateVolumeMesh(VolumeMesh &vMesh);

  vector< Array<Point>* > refinedPoint;
  vector< Array<Edge>* > refinedEdge;
  vector< Array<Triangle>* > refinedTriangle;
  
  unsigned long getTriangleCount();
  unsigned long getPointCount();

  SegmentSet segment;
};

void getEdgesAtPoints(ShapeMesh &m, PointEdgeVectorMap& edgesAtPoint);

ostream& operator<<(ostream& os, const ShapeMesh &s);

Real angle(Edge *e1, Edge *e2);
void traverse(EdgeSet &edge,
              EdgeSet &boundary,
              TriangleSet &triangle,
              Triangle *t,
              TriangleSet &triangleTraversed,
              EdgeSet &edgeTraversed,
              EdgeSet &boundaryEdges,
              TriangleSet *triangleAllowed = NULL);

#endif
