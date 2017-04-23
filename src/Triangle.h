#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Vector3.h"

struct Edge;

struct Triangle {
  Point *p1;
  Point *p2;
  Point *p3;
  Edge *e12;
  Edge *e23;
  Edge *e31;
  Key index;
};

Vector3 normalUnnormalized(Triangle *t);
Vector3 normal(Triangle *t);
Real area(Triangle *t);
Point center(Triangle *t);

inline int sharedpoints(Triangle *t1, Triangle *t2)
{
  int c = 0;
  if((t1->p1 == t2->p1) || (t1->p1 == t2->p2) || (t1->p1 == t2->p3)) c++;
  if((t1->p2 == t2->p1) || (t1->p2 == t2->p2) || (t1->p2 == t2->p3)) c++;
  if((t1->p3 == t2->p1) || (t1->p3 == t2->p2) || (t1->p3 == t2->p3)) c++;
  return c;
}

ostream& operator<<(ostream& os, const Triangle &t);

#endif
