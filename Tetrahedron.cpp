#include "Tetrahedron.h"
#include "MathUtils.h"

Real volume(Tetrahedron *t)
{
  return ONESIXTH * sqrt(norm( dot( *(t->p2) - *(t->p1), (*(t->p3) - *(t->p1)) * (*(t->p4) - *(t->p1)))));
}

Point center(Tetrahedron *t)
{
  Point p = 0.25 * (*(t->p1) + *(t->p2) + *(t->p3) + *(t->p4));
  return p;
}

ostream& operator<<(ostream& os, const Tetrahedron &t)
{
  os << "[ " << *(t.p1) << " " << *(t.p2) << " " << *(t.p3) << " " << *(t.p4) << "]";
  return os;
}
