#include "Transform.h"

Transform :: Transform() : C(0.,0.,0.), R(0.), S(1.,1.,1.), SR(0.), T(0.,0.,0.)
{
  bNull = true;
}

void Transform :: setC(Vector3 &C)
{
  this->C = C;
  if(C.x != 0. || C.y != 0. || C.z != 0.)
    bNull = false;
}

void Transform :: setR(Rotation &R)
{
  this->R = R;
  if(R.phi != 0.)
    bNull = false;
}

void Transform :: setS(Vector3 &S)
{
  this->S = S;
  if(S.x != 1. || S.y != 1. || S.z != 1.)
    bNull = false;
}

void Transform :: setSR(Rotation &SR)
{
  this->SR = SR;
  if(SR.phi != 0.)
    bNull = false;
}

void Transform :: setT(Vector3 &T)
{
  this->T = T;
  if(T.x != 0. || T.y != 0. || T.z != 0.)
    bNull = false;
}

Vector3 transform(vector<Transform> t, Vector3 &v)
{
  Vector3 r = v;
  int size = t.size();
  for(int i=0;i<size;i++) {
    r = t[i].transform(r);
  }
  return r;
}

// From x3d standard 
// P' = T * C * R * SR * S * -SR * -C * P
Vector3 Transform :: transform(Vector3 &p)
{
  if(bNull)
    return p;
  else {
    Vector3 r;
    r = p - C;
    r = (-SR).rotate(r);
    r = S.scale(r);
    r = SR.rotate(r);
    r = R.rotate(r);
    r = r + C;
    r = r + T;
    return r;
  }
}
