#ifndef VOLUMEMESH_H
#define VOLUMEMESH_H

#include "Tetrahedron.h"
#include "Vector3.h"
#include "Array.h"

class VolumeMesh {
 public:
  Array<Tetrahedron> tetrahedra;
  Array<Point> point;
};

ostream& operator<<(ostream& os, const VolumeMesh &m);

#endif
