#ifndef EDGE_H
#define EDGE_H

#include "defs.h"
#include "Vector3.h"

#include <vector>
using namespace std;

struct Triangle;

struct Edge {
  Point *p1;
  Point *p2;
  Real length;
  Triangle *t1;
  Triangle *t2;
  Key index;
};

typedef UnorderedMap<Edge*, int> EdgeIntMap;

#endif
