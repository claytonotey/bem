#include "ComplexVector3.h"
#include "Vector3.h"
#include "defs.h"
#include "MathUtils.h"
#include "Impedance.h"
#include "GreensFunction.h"
#include "RWG.h"
#include "ShapeParser.h"
#include "ShapeMesh.h"
#include "VolumeMesh.h"
#include <time.h>
#include "TriangleQuadrature.h"
#include "TetrahedronQuadrature.h"
#include "AnalyticIntegrals.h"
#include "Options.h"
#include "MatrixOperations.h"
#include "Source.h"
#include "Output.h"
#include "Flux.h"
#include "FrequencySweep.h"
#include "SI.h"

bool refineEdgeFuncBody(Edge *e, void *args)
{
  return false;
}


int main(int argc, char **argv)
{
  initReferenceCountingPointerLock();
  ShapeParser parser;
  char *body1x3d = argv[1];
  char *body2x3d = argv[2];
  Real a = atof(argv[3]);
  char *sweepFile = argv[4];

  Shape body1Shape;
  Dielectric *diel1;
  parser.parse(body1x3d,body1Shape,&diel1);
  Shape body2Shape;
  Dielectric *diel2;
  parser.parse(body2x3d,body2Shape,&diel2);

  fprintf(stderr,"Refining surface mesh\n");

  Options options;

  ShapeMesh body1Mesh(body1Shape,
                      options.sharpThresh,
                      options.concaveThresh,
                      options.relativeCostAngleConcave,
                      NULL,
                      NULL);

  Real maxEdgeLength;
  ShapeMesh refinedBody1Mesh(body1Shape,
                             options.sharpThresh,
                             options.concaveThresh,
                             options.relativeCostAngleConcave,
                             &refineEdgeFuncBody,
                             &maxEdgeLength);
  
  fprintf(stderr,"Generating volume mesh\n");
  VolumeMesh vMesh;
  body1Mesh.generateVolumeMesh(vMesh);
  
  ShapeMesh refinedBody2Mesh(body2Shape,
                             options.sharpThresh,
                             options.concaveThresh,
                             options.relativeCostAngleConcave,
                             &refineEdgeFuncBody,
                             &maxEdgeLength);


  fprintf(stderr,"body 1: %ld original -> %ld refined surface points\n", body1Mesh.getPointCount(),refinedBody1Mesh.getPointCount());
  fprintf(stderr,"body 2: %ld refined surface points\n", refinedBody2Mesh.getPointCount());
  fprintf(stderr,"%d volume points\n", vMesh.point.n);
  fprintf(stderr,"%d volume elements\n", vMesh.tetrahedra.n);


  ShapeMeshDielectricPairList meshDiel;
  PointEdgeVectorMap eap;
  getEdgesAtPoints(refinedBody1Mesh,eap);
  getEdgesAtPoints(refinedBody2Mesh,eap);
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody1Mesh,diel1));
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody2Mesh,diel2));
  Real minSegmentRadius = 0.0;
  RWGTree rwgTree(meshDiel,eap,Options::minSegmentSize,minSegmentRadius);
  RWGTree *s = rwgTree.sub[0];
  RWGTree *t = rwgTree.sub[1];


  fprintf(stderr,"%d unknowns\n", 2*rwgTree.rwgs);

  FrequencySweep sweep(sweepFile);

  GreensFunction :: setLengthScale(a);

  Vacuum vac;
  GreensFunction greenExt(&vac,true);

  for(vector<Real>::iterator i=sweep.f.begin(); i != sweep.f.end(); ++i) {
    Real eV = *i;


    greenExt.init(eV);
    s->g->init(eV);
    t->g->init(eV);

    cerr << "eV = " << eV << ", omega = " << t->g->omega << "\n";
    cerr << "eps1 = " << t->g->eps << ", eps2 = " << s->g->eps << "\n";
    cerr << "k_T = " << s->g->k << "\n";
    cerr << "k_L = " << s->g->gLong->k << "\n";

    Progress progress;
    progress.startMatrixConstruction(2*rwgTree.rwgs,2*rwgTree.rwgs);



    ComplexMatrix *WH = Z(&rwgTree,&rwgTree,&greenExt,progress,options);
    ComplexMatrix *G1H = Z(s,s,s->g,progress,options);
    ComplexMatrix *G2H = Z(t,t,t->g,progress,options);

    ComplexDenseMatrix W(*WH);
    ComplexDenseMatrix G1(*G1H);
    ComplexDenseMatrix G2(*G2H);
 
    for(int i=0; i<2*s->rwgs;i++)
      for(int j=0; j<2*s->rwgs;j++)
        (W)(i,j) += (G1)(i,j);

    int i0 = 2*s->rwgs;
    for(int i=0; i<2*t->rwgs;i++)
      for(int j=0; j<2*t->rwgs;j++)
        (W)(i0+i,i0+j) += (G2)(i,j);



    inverse(W);
    ComplexDenseMatrix W1(W,0,0,2*s->rwgs,2*(s->rwgs+t->rwgs));

    symm(G1);
    ComplexDenseMatrix B1 = G1 * W1;
    B1.resize(0,2*s->rwgs,2*s->rwgs,2*t->rwgs);
    //delete G1;

    ComplexDenseMatrix W2(W,0,2*s->rwgs,2*(s->rwgs+t->rwgs),2*t->rwgs);
    symm(G2);
    ComplexDenseMatrix B2 = W2 * G2;
    B2.resize(0,0,2*s->rwgs,2*t->rwgs);
    //delete G2;

    multiply_Adagger_B(B2,B1,W);

    Complex S = 0.25*trace(W);


    delete G1H;
    delete G2H;
    delete WH;

    cout << eV << " " << S << "\n";
  }
  destroyReferenceCountingPointerLock();
}

