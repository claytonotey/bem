#ifndef ZROTATION_H
#define ZROTATION_H

#include "Vector3.h"

class ZRotation {
public:
  ZRotation() {}
  ZRotation(Vector3 &z);
  void rotate(Vector3 &v, Real *x0, Real *y0, Real *z0_ = NULL);
  Vector3 rotate(Vector3 &v);

  Real lxz, lyz, zz;
  Real cosphiy, sinphiy;
  Real cosphix, sinphix;
};

#endif
