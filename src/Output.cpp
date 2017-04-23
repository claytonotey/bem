#include "Output.h"

void outputJMMagnitude(RWGTree *top, RWGTree *t, ComplexDenseMatrix &X)
{
  int m = t->rwgs;
  EdgeIntMap edge;

  for(int i=0; i<m; i++) {
    RWG *rwg = t->rwg[i];
    if(t->index.size())
      edge[rwg->e] = t->index[i];
    else
      edge[rwg->e] = i;
  }

  cout << "[\n";
  for(TriangleSet::iterator j = t->triangle.begin(); j != t->triangle.end(); ++j) {
    Triangle *t1 = *j;
    Real A2i = 1. / (2. * area(t1));

    ComplexVector3 j1(0.,0.,0.);
    ComplexVector3 j2(0.,0.,0.);
    ComplexVector3 j3(0.,0.,0.);
    ComplexVector3 m1(0.,0.,0.);
    ComplexVector3 m2(0.,0.,0.);
    ComplexVector3 m3(0.,0.,0.);

    int i12 = edge[t1->e12];
    RWG *rwg12 = top->rwg[i12];
    int ii12 = 2*i12;
    Real c12 = (t1==rwg12->pos.t?rwg12->e->length:-rwg12->e->length) * A2i;
    Complex j12 = c12 * X(ii12,0);
    Complex m12 = c12 * X(ii12+1,0);
    j1 += j12 * (*(t1->p1) - *(t1->p3));
    j2 += j12 * (*(t1->p2) - *(t1->p3));
    m1 += m12 * (*(t1->p1) - *(t1->p3));
    m2 += m12 * (*(t1->p2) - *(t1->p3));

    int i23 = edge[t1->e23];
    RWG *rwg23 = top->rwg[i23];
    int ii23 = 2*i23;
    Real c23 = (t1==rwg23->pos.t?rwg23->e->length:-rwg23->e->length) * A2i;
    Complex j23 = c23 * X(ii23,0);
    Complex m23 = c23 * X(ii23+1,0);
    j2 += j23 * (*(t1->p2) - *(t1->p1));
    j3 += j23 * (*(t1->p3) - *(t1->p1));
    m2 += m23 * (*(t1->p2) - *(t1->p1));
    m3 += m23 * (*(t1->p3) - *(t1->p1));

    int i31 = edge[t1->e31];
    RWG *rwg31 = top->rwg[i31];
    int ii31 = 2*i31;
    Real c31 = (t1==rwg31->pos.t?rwg31->e->length:-rwg31->e->length) * A2i;
    Complex j31 = c31 * X(ii31,0);
    Complex m31 = c31 * X(ii31+1,0);
    j3 += j31 * (*(t1->p3) - *(t1->p2));
    j1 += j31 * (*(t1->p1) - *(t1->p2));
    m3 += m31 * (*(t1->p3) - *(t1->p2));
    m1 += m31 * (*(t1->p1) - *(t1->p2));

    Real j1m = sqrt(norm(j1));
    Real j2m = sqrt(norm(j2));
    Real j3m = sqrt(norm(j3));
    Real m1m = sqrt(norm(m1));
    Real m2m = sqrt(norm(m2));
    Real m3m = sqrt(norm(m3));

    cout << "[" << *(t1->p1) << " " << *(t1->p2) << " " << *(t1->p3) << " [" << j1m << " " << j2m << " " << j3m << "] [" << m1m << " " << m2m << " " << m3m << "]];\n";
  }
  cout << "]\n";
}
