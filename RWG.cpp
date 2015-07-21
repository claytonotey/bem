#include "RWG.h"
#include "Triangle.h"
#include "MathUtils.h"

RWG :: RWG(Edge *e) : e(e), pos(e,e->t1,true), neg(e,e->t2,false) {}

HalfRWG :: HalfRWG(Edge *e, Triangle *t, bool bPositive) : e(e), t(t), bPositive(bPositive)
{
  if((t->p1 != e->p1) && (t->p1 != e->p2)) {
    p0 = t->p1;
    p1 = (t->e12->p1 == p0) ? (t->e12->p2) : (t->e12->p1);
    p2 = (t->e31->p1 == p0) ? (t->e31->p2) : (t->e31->p1);
  } else if((t->p2 != e->p1) && (t->p2 != e->p2)) {
    p0 = t->p2;
    p1 = (t->e23->p1 == p0) ? (t->e23->p2) : (t->e23->p1);
    p2 = (t->e12->p1 == p0) ? (t->e12->p2) : (t->e12->p1);
  } else {
    p0 = t->p3;
    p1 = (t->e31->p1 == p0) ? (t->e31->p2) : (t->e31->p1);
    p2 = (t->e23->p1 == p0) ? (t->e23->p2) : (t->e23->p1);
  }
  v01 = *p1 - *p0;
  v12 = *p2 - *p1;  
  v20 = *p0 - *p2;
  v02 = *p2 - *p0;
  normv01i = 1. / sqrt(norm(v01));
  normv12i = 1. / sqrt(norm(v12));
  normv20i = 1. / sqrt(norm(v20));
  n = v01 * v02; 
  Real A2 = sqrt(norm(n));
  n *= 1. / A2;
  Ai = 2.0 / A2;
  m01 = v01 * n;
  m01.normalize();
  m12 = v12 * n;
  m12.normalize();
  m20 = v20 * n;
  m20.normalize();

}
