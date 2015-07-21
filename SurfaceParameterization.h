#ifndef SURFACEPARAMETERIZATION_H
#define SURFACEPARAMETERIZATION_H

#include "defs.h"
#include "MathUtils.h"
#include "Vector3.h"
#include "ShapeMesh.h"
#include "SparseMatrix.h"

class Vector2 {
 public:
  Real x;
  Real y;
  Vector2() {}
  Vector2(Real x, Real y) { this->x = x; this->y = y; }
};

typedef Vector2 Point2;

#include <stack>
#include <vector>
#include <Map>
#include <utility>
using namespace std;

typedef pair<int,int> IntPair;
typedef UnorderedMap< Triangle*, Real> TriangleRealMap;
typedef UnorderedMap<IntPair, Real> IntPairRealMap;
typedef vector<Edge*> EdgeVector;
typedef UnorderedSet<PointSet*> PointSetSet;
typedef UnorderedMap<Point*, int> PointIntMap;
typedef UnorderedMap<Point*, Point2> PointPoint2Map;

inline Real norm(Point2 &p1)
{
  return Square(p1.x) + Square(p1.y);
}

inline Real dot(Point2 &p1, Point2 &p2)
{
  return p1.x * p2.x + p1.y * p2.y;
}

inline Point2 operator+(Point2 &p1, Point2 &p2)
{
  Point2 r;
  r.x = p1.x + p2.x;
  r.y = p1.y + p2.y;
  return r;
}

inline Point2 operator-(Point2 &p1, Point2 &p2)
{
  Point2 r;
  r.x = p1.x - p2.x;
  r.y = p1.y - p2.y;
  return r;
}

inline void operator*=(Point2 &p, Real c)
{
  p.x *= c;
  p.y *= c;
}

inline Point2 operator*(Real c, Point2 &p)
{
  Point2 r;
  r.x = c * p.x;
  r.y = c * p.y;
  return r;
}

inline void operator+=(Point2 &p1, Point2 &p2)
{
  p1.x += p2.x;
  p1.y += p2.y;
}

inline void operator-=(Point2 &p1, Point2 &p2)
{
  p1.x -= p2.x;
  p1.y -= p2.y;
}

class TriangleQuadTree {
 public:
  TriangleQuadTree(TriangleSet &triangle, PointPoint2Map &u, Real x0, Real y0, Real w, Real h);
  TriangleQuadTree(TriangleSet &triangle, PointPoint2Map &u);
  void init(TriangleSet &triangle, PointPoint2Map &u, Real x0, Real y0, Real w, Real h);
  ~TriangleQuadTree();
  Triangle* find(Real x, Real y);

 protected:
  TriangleQuadTree *sub00;
  TriangleQuadTree *sub01;
  TriangleQuadTree *sub10;
  TriangleQuadTree *sub11;
  Triangle *t;
  Real x0, y0;
};

class SurfaceParameterization {
 public:
  SurfaceParameterization(PointSet &point,
                          PointSet &boundary,
                          PointEdgeVectorMap &edgesAtPoint,
                          PointEdgeVectorMap &boundaryEdgesAtPoint,
                          TriangleSet &triangle);
  ~SurfaceParameterization();
  Point inverse(Point2 &p2);
  friend ostream& operator<<(ostream& os, SurfaceParameterization &u);
  Point2& operator[](Point *p);
  static int hexagonalPoints(int n);
  static Point2 hexagonalPoint(int n, int p, int q);

 protected:
  TriangleQuadTree *quad;
  PointPoint2Map u;
  TriangleSet &triangle;
};


ostream& operator<<(ostream& os, PointPoint2Map &u);
ostream& operator<<(ostream& os, Point2 &p2);

#endif
