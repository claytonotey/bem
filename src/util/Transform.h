#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Util.h"
#include "Vector3.h"
#include "Rotation.h"

using namespace std;
#include <vector>

class Transform {
 public:
  Transform();

  Vector3 transform(Vector3 &v);
  friend Vector3 transform(vector<Transform> t, Vector3 &v);

  void setC(const Vector3 &C);
  void setR(const Rotation &R);
  void setS(const Vector3 &S);
  void setSR(const Rotation &SR);
  void setT(const Vector3 &T);

 protected:
  bool bNull;
  Vector3 C; //center
  Rotation R; //rotation
  Vector3 S; //scale
  Rotation SR; //scaleOrientation;
  Vector3 T; //translation;
};

#endif
