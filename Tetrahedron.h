#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "Vector3.h"

class Tetrahedron {
 public:
  Point *p1;
  Point *p2;
  Point *p3;
  Point *p4;
};

Point center(Tetrahedron *t);
Real volume(Tetrahedron *t);
ostream& operator<<(ostream& os, const Tetrahedron &t);

#endif
