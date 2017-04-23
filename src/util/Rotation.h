#ifndef ROTATION_H
#define ROTATION_H

#include "Vector3.h"

class Rotation {
 public:
  Rotation(char *str);
  Rotation(Vector3 &v, Real s);
  Rotation(Real s);
  Rotation(Real x, Real y, Real z, Real s);
  Rotation operator-();
  Vector3 rotate(Vector3 &v);
  friend ostream& operator<<(ostream& os, const Rotation &r);

  Real phi;
  Real c,s,t;
  Real x,y,z;

 protected:
  void init();
};


class RotationException: public exception
{
 public:
  RotationException() {};
  virtual const char* what() const throw()
  {
    return "Error parsing Rotation";
  }
};

#endif
