#include "Triangle.h"
#include "Edge.h"
#include "MathUtils.h"

Vector3 normalUnnormalized(Triangle *t)
{
  Vector3 *v0 = t->p1;
  Vector3 v1 = (t->e12->p1 == v0) ? *(t->e12->p2) - *(t->e12->p1) : *(t->e12->p1) - *(t->e12->p2);
  Vector3 v2 = (t->e31->p1 == v0) ? *(t->e31->p2) - *(t->e31->p1) : *(t->e31->p1) - *(t->e31->p2);
  Vector3 n = v1 * v2;
  return n;
}

Vector3 normal(Triangle *t)
{
  Vector3 n = normalUnnormalized(t);
  n.normalize();
  return n;
}

Real area(Triangle *t)
{
  Vector3 *v0 = t->p1;
  Vector3 v1 = (t->e12->p1 == v0) ? *(t->e12->p2) - *(t->e12->p1) : *(t->e12->p1) - *(t->e12->p2);
  Vector3 v2 = (t->e31->p1 == v0) ? *(t->e31->p2) - *(t->e31->p1) : *(t->e31->p1) - *(t->e31->p2);
  Vector3 n = v1 * v2;
  return 0.5 * sqrt(norm(n));
}

Point center(Triangle *t)
{
  Point p = *(t->p1) + *(t->p2) + *(t->p3); p *= ONETHIRD;
  return p;
}

ostream& operator<<(ostream& os, const Triangle &t)
{
  os << "[" << *(t.p1) << " " << *(t.p2) << " " << *(t.p3) << "]" ;
  return os;
}
