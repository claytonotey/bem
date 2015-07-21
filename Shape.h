#ifndef SHAPE_H
#define SHAPE_H

#include <vector>
using namespace std;

#include "Vector3.h"
#include "FaceIndex.h"

class Shape {
 public:
  Shape();

  vector<FaceIndex> faceIndex;
  vector<Vector3> coord;
};

ostream& operator<<(ostream& os, const Shape &s);

#endif
