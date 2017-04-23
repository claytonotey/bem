#include "Util.h"
#include "ZRotation.h"

ZRotation :: ZRotation(Vector3 &z)
{
  lxz = sqrt(z.x*z.x+z.z*z.z);
  cosphiy = z.z/lxz;
  sinphiy = z.x/lxz;
  if(lxz == 0.0)
    zz = z.z;
  else
    zz = sinphiy * z.x + cosphiy * z.z;
  lyz = sqrt(z.y*z.y+zz*zz);
  cosphix = zz/lyz;
  sinphix = z.y/lyz;
}

// rotate v so z = (0,0,1).  The result is in x0, y0, z0.
void ZRotation :: rotate(Vector3 &v, Real *x0, Real *y0, Real *z0_)
{
  //rotate about y
  Real z0;
  if(lxz == 0.) {
    *x0 = v.x;
    z0 = v.z;
  } else {
    *x0 = cosphiy * v.x - sinphiy * v.z;
    z0 = sinphiy * v.x + cosphiy * v.z;
  }

  //rotate about x
  if(lyz == 0.) {
    *y0 = v.y;
  } else {
    *y0 = cosphix * v.y - sinphix * z0;
    z0 = sinphix * v.y + cosphix * z0;
  }
  if(z0_) *z0_ = z0;
}

Vector3 ZRotation :: rotate(Vector3 &v)
{
  Real x,y,z;
  rotate(v,&x,&y,&z);
  Vector3 r(x,y,z);
  return r;
}
