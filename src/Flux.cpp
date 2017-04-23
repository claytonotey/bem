#include "Flux.h"
#include "DenseMatrix.h"
#include "MathUtils.h"
#include "Vector3.h"
#include "ComplexVector3.h"
#include "FunctionArgs.h"
#include "TriangleQuadrature.h"
#include "TetrahedronQuadrature.h"
#include "Array.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Options.h"
#include "Source.h"
#include "Progress.h"
#include "Field.h"
#include "Edge.h"


#if HAVE_OMP
#include <omp.h>
#endif

void enumFluxes(RWGTree *top, RWGTree *t)
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

  for(TriangleSet::iterator j = t->triangle.begin(); j != t->triangle.end(); ++j) {
    Triangle *t1 = *j;

    Vector3 n = normalUnnormalized(t1);
    Real A2 = sqrt(norm(n));
    n *= 1. / A2;

    Vector3 v12 = (*(t1->p2) - *(t1->p1)); 
    Vector3 v23 = (*(t1->p3) - *(t1->p2)); 
    Vector3 v31 = (*(t1->p1) - *(t1->p3)); 

    int i12 = edge[t1->e12];
    RWG *rwg12 = top->rwg[i12];
    int ij12 = 2*i12;
    int im12 = ij12 + 1;
    Real c12 = (t1==rwg12->pos.t?rwg12->e->length:-rwg12->e->length);

    int i23 = edge[t1->e23];
    RWG *rwg23 = top->rwg[i23];
    int ij23 = 2*i23;
    int im23 = ij23 + 1;
    Real c23 = (t1==rwg23->pos.t?rwg23->e->length:-rwg23->e->length);

    int i31 = edge[t1->e31];
    RWG *rwg31 = top->rwg[i31];
    int ij31 = 2*i31;
    int im31 = ij31 + 1;
    Real c31 = (t1==rwg31->pos.t?rwg31->e->length:-rwg31->e->length);
    
    Real c = 1.0 / (48. * A2);

    Real s1223 = c * (c12 * c23) * dot(n,v31 * (v12 - v23));
    cout << ij12 << " " << im23 << " " << s1223  << "\n";
    Real s3112 = c * (c31 * c12) * dot(n,(v12 - v31) * v23);
    cout << ij31 << " " << im12 << " " << s3112  << "\n";

    Real s1231 = c * (c12 * c31) * dot(n,v23 * (v12 - v31));
    cout << ij12 << " " << im31 << " " <<  s1231  << "\n";
    Real s2312 = c * (c23 * c12) * dot(n,(v12 - v23) * v31);
    cout << ij23 << " " << im12 << " " <<  s2312  << "\n";

    Real s2331 = c * (c23 * c31) * dot(n,v12 * (v23 - v31));
    cout << ij23 << " " << im31 << " " << s2331  << "\n";
    Real s3123 = c * (c31 * c23) * dot(n,(v23 - v31) * v12);
    cout << ij31 << " " << im23 << " " << s3123  << "\n";


    //j1 m1
    s1231 = -2. * c * (c12 * c31) * dot(n,v31 * v12);
    cout << ij12 << " " << im31 << " " << s1231  << "\n";
    cout << ij31 << " " << im12 << " " <<  -s1231  << "\n";
    
    //j2 m2 
    s2312 = -2. * c * (c23 * c12) * dot(n,v12 * v23);
    cout << ij23 << " " << im12 << " " << s2312  << "\n";
    cout << ij12 << " " << im23 << " " <<  -s2312  << "\n";

    //j3 m3
    s3123 = -2. * c * (c31 * c23) * dot(n,v23 * v31);
    cout << ij31 << " " << im23 << " " << s3123  << "\n";
    cout << ij23 << " " << im31 << " " <<  -s3123  << "\n";
  }
}

Real flux(RWGTree *top, RWGTree *t, const ComplexDenseMatrix &X)
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

  Real S = 0.;
  for(TriangleSet::iterator j = t->triangle.begin(); j != t->triangle.end(); ++j) {
    Triangle *t1 = *j;

    Vector3 n = normalUnnormalized(t1);
    Real A2 = sqrt(norm(n));
    n *= 1. / A2;

    Vector3 v12 = (*(t1->p2) - *(t1->p1)); 
    Vector3 v23 = (*(t1->p3) - *(t1->p2)); 
    Vector3 v31 = (*(t1->p1) - *(t1->p3)); 

    ComplexVector3 j1(0.,0.,0.);
    ComplexVector3 j2(0.,0.,0.);
    ComplexVector3 j3(0.,0.,0.);
    ComplexVector3 m1(0.,0.,0.);
    ComplexVector3 m2(0.,0.,0.);
    ComplexVector3 m3(0.,0.,0.);

    int i12 = edge[t1->e12];
    RWG *rwg12 = top->rwg[i12];
    int ij12 = 2*i12;
    int im12 = ij12 + 1;
    Real c12 = (t1==rwg12->pos.t?rwg12->e->length:-rwg12->e->length);
    Complex j12 = c12 * X(ij12,0);
    Complex m12 = c12 * X(im12,0);

    j1 += j12 * v31;
    j2 -= j12 * v23;
    m1 += m12 * v31;
    m2 -= m12 * v23;

    int i23 = edge[t1->e23];
    RWG *rwg23 = top->rwg[i23];
    int ij23 = 2*i23;
    int im23 = ij23 + 1;
    Real c23 = (t1==rwg23->pos.t?rwg23->e->length:-rwg23->e->length);
    Complex j23 = c23 * X(ij23,0);
    Complex m23 = c23 * X(im23,0);
    j2 += j23 * v12;
    j3 -= j23 * v31;
    m2 += m23 * v12;
    m3 -= m23 * v31;

    int i31 = edge[t1->e31];
    RWG *rwg31 = top->rwg[i31];
    int ij31 = 2*i31;
    int im31 = ij31 + 1;
    Real c31 = (t1==rwg31->pos.t?rwg31->e->length:-rwg31->e->length);
    Complex j31 = c31 * X(ij31,0);
    Complex m31 = c31 * X(im31,0);
    j3 += j31 * v23;
    j1 -= j31 * v12;
    m3 += m31 * v23;
    m1 -= m31 * v12;
    
    Complex omegai = Complex(1.,0.) / t->g->omega;
    m1 *= omegai;
    m2 *= omegai;
    m3 *= omegai;
    j1 = conj(j1);
    j2 = conj(j2);
    j3 = conj(j3);
    
    Real dS = dot(n,real(j1 * m1) + real(j2 * m2) + real(j3 * m3) + real((j1+j2+j3) * (m1+m2+m3))) / (48. * A2);

    Real dS2 = 0.0;
    Real c = 1.0 / (48. * A2 * real(t->g->omega));

    Real s1223 = c * (c12 * c23) * dot(n,v31 * (v12 - v23));
    dS2 += s1223 * real(conj(X(ij12,0)) * X(im23,0));
    Real s3112 = c * (c31 * c12) * dot(n,(v12 - v31) * v23);
    dS2 += s3112 * real(conj(X(ij31,0)) * X(im12,0));

    Real s1231 = c * (c12 * c31) * dot(n,v23 * (v12 - v31));
    dS2 += s1231 * real(conj(X(ij12,0)) * X(im31,0));
    Real s2312 = c * (c23 * c12) * dot(n,(v12 - v23) * v31);
    dS2 += s2312 * real(conj(X(ij23,0)) * X(im12,0));

    Real s2331 = c * (c23 * c31) * dot(n,v12 * (v23 - v31));
    dS2 += s2331 * real(conj(X(ij23,0)) * X(im31,0));
    Real s3123 = c * (c31 * c23) * dot(n,(v23 - v31) * v12);
    dS2 += s3123 * real(conj(X(ij31,0)) * X(im23,0));


    //j1 m1
    s1231 = -2. * c * (c12 * c31) * dot(n,v31 * v12);
    dS2 += s1231 * real(conj(X(ij12,0)) * X(im31,0));
    dS2 += -s1231 * real(conj(X(ij31,0)) * X(im12,0));
    
    //j2 m2 
    s2312 = -2. * c * (c23 * c12) * dot(n,v12 * v23);
    dS2 += s2312 * real(conj(X(ij23,0)) * X(im12,0));
    dS2 += -s2312 * real(conj(X(ij12,0)) * X(im23,0));

    //j3 m3
    s3123 = -2. * c * (c31 * c23) * dot(n,v23 * v31);
    dS2 += s3123 * real(conj(X(ij31,0)) * X(im23,0));
    dS2 += -s3123 * real(conj(X(ij23,0)) * X(im31,0));

    cout << dS << " " << dS2 << "\n";

    S += dS;
  }
  
  return S;
}

Real flux(ComplexHMatrix &ZLU, Vector3 &j, GreensFunction *green, Point &p, RWGTree *rwgTree, RWGTree *t, GreensFunction *greenExt, Progress &progress, Options &options)
{
  SourceSet sources;
  DipolePointSource source(green,p,j);
  sources.insert(&source);
  ComplexHMatrix *B;
  B = (ComplexHMatrix*) BH(rwgTree,greenExt,sources,progress,options);
  ComplexDenseMatrix Bx = *B;
  ComplexHMatrix Bhx(Bx,NULL,NULL);
  //solveLU(ZLU,*B);
  solveLU(ZLU,Bhx);
  //ComplexDenseMatrix Bxx(*B);
  ComplexDenseMatrix Bxx(Bhx);
  Real S = flux(rwgTree,t,Bxx);
  delete B;
  return S;
}

void flux_integrand_tet(Vector3 *v, Real *x, FunctionArgs *args)
{
  RWGTree *rwgTree = args->rwgTree;
  RWGTree *t = args->t;
  Tetrahedron &tet = *(args->tet);
  ComplexHMatrix &ZLU = *(args->ZLU);
  GreensFunction *green = args->green;
  GreensFunction *greenExt = args->greenExt;

  Vector3 p = *(tet.p1) + (*v);

  Options options;
  Progress progress;
  Real S = 0;

  Vector3 jx(1.,0.,0.);
  Vector3 jy(0.,1.,0.);
  Vector3 jz(0.,0.,1.);
  S += flux(ZLU,jx,green,p,rwgTree,t,greenExt,progress,options);
  S += flux(ZLU,jy,green,p,rwgTree,t,greenExt,progress,options);
  S += flux(ZLU,jz,green,p,rwgTree,t,greenExt,progress,options);
  
  cerr << "S = " << S << "\n";
  *x = S;
}

Real flux(Tetrahedron *tet, ComplexHMatrix *ZLU, GreensFunction *greenExt, RWGTree *rwgTree, RWGTree *s, RWGTree *t, TetQuadratureFunctionFactory &tetquadfuncFactory)
{
  FunctionArgs args;
  args.tet = tet;
  args.ZLU = ZLU;
  args.rwgTree = rwgTree;
  args.green = s->g;
  args.greenExt = greenExt;
  args.t = t;

  Vector3 v2 = *(tet->p2)-*(tet->p1);
  Vector3 v3 = *(tet->p3)-*(tet->p1);
  Vector3 v4 = *(tet->p4)-*(tet->p1);

  TetQuadratureFunction<Real,RealValuedFunction>::Type tetquadfunc;
  tetquadfunc = tetquadfuncFactory.getQuadFunc(tet);
  Real S = tetquadfunc(&flux_integrand_tet, &args, NULL, v2, v3, v4);

  return volume(tet) * S;
}

Real totalFlux(VolumeMesh &vMesh, ComplexHMatrix *ZLU, GreensFunction *greenExt, RWGTree *rwgTree, RWGTree *s, RWGTree *t, TetQuadratureFunctionFactory &tetquadfuncFactory, Progress &progress)
{
  Real S = 0.;

  progress.startFluxCalculation(vMesh.tetrahedra.n);

  int i;
  int chunk = 128;

#pragma omp parallel for default(shared) private(i) schedule(static,chunk) reduction(+:S)
  for(i=0; i<vMesh.tetrahedra.n; i++) {
    Tetrahedron *tet = vMesh.tetrahedra + i;
    S = S + flux(tet,ZLU,greenExt,rwgTree,s,t,tetquadfuncFactory);
    progress.finishedFluxVolumeElement();
  }
  return S;
}
