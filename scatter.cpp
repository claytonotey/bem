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
#include <xmmintrin.h>

bool refineEdgeFuncBody(Edge *e, void *args)
{
  return false;
}


int main(int argc, char **argv)
{
  _mm_setcsr( _mm_getcsr() | 0x8040 );
  initReferenceCountingPointerLock();
  ShapeParser parser;
  char *body1x3d = argv[1];
  Real a = atof(argv[2]);
  char *sweepFile = argv[3];

  Shape body1Shape;
  Dielectric *diel1;
  parser.parse(body1x3d,body1Shape,&diel1);

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
  
  ShapeMeshDielectricPairList meshDiel;
  PointEdgeVectorMap eap;
  getEdgesAtPoints(refinedBody1Mesh,eap);
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody1Mesh,diel1));
  Real minSegmentRadius = 0.0;
  RWGTree rwgTree(meshDiel,eap,Options::minSegmentSize,minSegmentRadius);
  RWGTree *s = &rwgTree;

  fprintf(stderr,"%d unknowns\n", 2*rwgTree.rwgs);

  FrequencySweep sweep(sweepFile);

  GreensFunction :: setLengthScale(a);

  ConstantDielectric ext(2.3,1.0);
  //ConstantDielectric ext(1.0,1.0);
  //Vacuum vac;
  GreensFunction greenExt(&ext,true);

  for(vector<Real>::iterator i=sweep.f.begin(); i != sweep.f.end(); ++i) {
    Real eV = *i;


    greenExt.init(eV);
    s->g->init(eV);

    cerr << "eV = " << eV << ", omega = " << s->g->omega << "\n";
    cerr << "eps1 = " << s->g->eps << "\n";
    cerr << "k_T = " << s->g->k << "\n";
    if(s->g->gLong) cerr << "k_L = " << s->g->gLong->k << "\n";
    cerr << "Zeps(T) = " << s->g->Zeps << "\n";
    if(s->g->gLong) cerr << "Zeps(L) = " << s->g->gLong->Zeps << "\n";

    Progress progress;
    progress.startMatrixConstruction(2*rwgTree.rwgs,2*rwgTree.rwgs);
     ComplexMatrix *G1H = Z(s,s,s->g,progress,options);
    //ComplexMatrix *G1H = Zts(s,s,s->g,options);

    progress.startMatrixConstruction(2*rwgTree.rwgs,2*rwgTree.rwgs);
    ComplexMatrix *G0H = Z(s,s,&greenExt,progress,options);
    //ComplexMatrix *G0H = Zts(s,s,&greenExt,options);



    ComplexDenseMatrix G0(*G0H);
    ComplexDenseMatrix G1(*G1H); 
    ComplexDenseMatrix W(2*s->rwgs,2*s->rwgs);

    for(int i=0; i<2*s->rwgs;i++)
      for(int j=0; j<2*s->rwgs;j++)
        (W)(i,j) = (G0)(i,j) + (G1)(i,j);

    SourceSet sources;
    Vector3 k(0,0,1);
    ComplexVector3 E0(1,0,0);
    PlaneWaveSource plane(&greenExt,k,E0);
    sources.insert(&plane);
    ComplexDenseMatrix *B = Bdense(s, &greenExt, sources, progress, options);
    W.LU();
    solveLU(W,*B);
    //ComplexDenseMatrix JM = W * (*B);

    Real sigma = extinctionCrossSection(s, *B, &greenExt, &plane) * Square(a);

    //Vector3 r(.3,.5,0.812403);
    //r.normalize();
    Vector3 r = center(s->rwg[0]->pos.t);
    Vector3 n = r;
    Real rr = n.normalize();

    cerr << r << " " << s->rwg[0]->pos.n << "\n";

    //Vector3 r(0,1,0);
    
    for(int i=0; i<48; i++) {
      Real d = 0.9 + 0.2 * i / 47.0;
      Vector3 p = d * r;
      if(sqrt(norm(p)) < rr) {
        ComplexVector3Pair EHt = getField(sources, s, *B, s->g, p, false, options);
        ComplexVector3Pair EHl = getField(sources, s, *B, s->g, p, true, options);
        cout << eV << " " << d << " " << EHt.first << " " << EHl.first << s->g->eps << " " << n << "" << s->g->eps0 << "\n";
        cout << eV << " " << d << " " << EHt.second << " " << EHl.second << s->g->eps << " " << n << " " << s->g->eps0 << "\n";
      } else {
        ComplexVector3Pair EHt = getField(sources, s, *B, &greenExt, p, false, options);
        ComplexVector3Pair EHl = getField(sources, s, *B, &greenExt, p, true, options);
        cout << eV << " " << d << " " << EHt.first << " " << EHl.first << s->g->eps << " " << n << " " << s->g->eps0 << "\n";
        cout << eV << " " << d << " " << EHt.second << " " << EHl.second << s->g->eps << " " << n << " " << s->g->eps0 << "\n";
      }
    }
    
    /*
    Vector3 pn = 0.995 * r;
    Vector3 pp = 1.005 * r;


    Real dr = .001;
    Real dri = 1.0 / dr;
    Vector3 pn2 = pn-dr*r;
    Vector3 pp2 = pp+dr*r;

    ComplexVector3Pair EHnt = getField(sources, s, *B, s->g, pn, false, options);
    ComplexVector3Pair EHpt = getField(sources, s, *B, &greenExt, pp, false, options);

    ComplexVector3Pair EHnl = getField(sources, s, *B, s->g, pn, true, options);
    ComplexVector3Pair EHpl = getField(sources, s, *B, &greenExt, pp, true, options);

    ComplexVector3Pair EHnt2 = getField(sources, s, *B, s->g, pn2, false, options);
    ComplexVector3Pair EHpt2 = getField(sources, s, *B, &greenExt, pp2, false, options);

    ComplexVector3Pair EHnl2 = getField(sources, s, *B, s->g, pn2, true, options);
    ComplexVector3Pair EHpl2 = getField(sources, s, *B, &greenExt, pp2, true, options);

    ComplexVector3Pair EHn = EHnt;
    EHn += EHnl;

        cout << "E0 = " << EHpt.first - dot(n,EHpt.first) * n << "\n";
        cout << "E1t = " << EHnt.first - dot(n,EHnt.first) * n << "\n";
        cout << "E1l = " << EHnl.first - dot(n,EHnl.first) * n << "\n";
        cout << "E1 = " << EHn.first - dot(n,EHn.first) * n << "\n";

        //cout << "eps E1n = " << s->g->eps0 * dot(r,EHn.first) << "\n";
        //cout << "eps0, eps =  " << s->g->eps0 << ", " << s->g->eps << "\n";

    //H1 = " << EHn.second << "\n";
    //cout << "    E0 = " << (r * EHp.first) << "\n";
        //cout << "    E0n / E1n= " << dot(r,EHp.first) / dot(r,EHn.first) << "\n";
    //cout << "    H0 = " << EHp.second << "\n";
    printf("\n");
    
    */
    //cout << eV << " " << sigma << "\n";
    delete G0H;
    delete G1H;
    delete B;
  }
  destroyReferenceCountingPointerLock();
}
