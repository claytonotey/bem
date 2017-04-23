#include "Util.h"
#include "Rotation.h"

void Rotation :: init()
{
  s = sin(phi);
  c = cos(phi);
  t = 1. - c;
}

Rotation :: Rotation(char *str) 
{
  char *s = str;

  while(*s == ' ') s++; 
  if(*s == 0) throw RotationException();
  x = atof(s);
  while(*s != ' ') s++; 

  while(*s == ' ') s++; 
  if(*s == 0) throw RotationException();
  y = atof(s);
  while(*s != ' ') s++; 

  while(*s == ' ') s++; 
  if(*s == 0) throw RotationException();
  z = atof(s);
  while(*s != ' ') s++; 

  while(*s == ' ') s++; 
  if(*s == 0) throw RotationException();
  phi = atof(s);

  init();
}

Rotation :: Rotation(Real x, Real y, Real z, Real s)
{
  this->x = x;
  this->y = y;
  this->z = z;
  phi = s;
  init();
}

Rotation :: Rotation(Vector3 &v, Real s)
{
  x = v.x;
  y = v.y;
  z = v.z;
  phi = s;
  init();
}

Rotation :: Rotation(Real s)
{
  x = 0.;
  y = 0.;
  z = 1.;
  phi = s;
  init();
}

Rotation Rotation :: operator-()
{
  Rotation r(x,y,z,-phi);
  return r;
}

Vector3 Rotation :: rotate(Vector3 &v)
{
  Vector3 r;
  r.x = (t * x * x + c) * v.x + (t * x * y + s * z) * v.y + (t * x * z - s * y) * v.z;
  r.y = (t * x * y - s * z) * v.x + (t * y * y + c) * v.y + (t * y * z + s * x) * v.z;
  r.z = (t * x * z + s * y) * v.x + (t * y * z - s * x) * v.y + (t * z * z + c) * v.z;
  return r;
}

ostream& operator<<(ostream& os, const Rotation &v) {
  os << "[" << v.x << " " << v.y << " " << v.z << " " << v.phi << "]" ;
  return os;
}
