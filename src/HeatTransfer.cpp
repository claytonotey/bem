#include "HeatTransfer.h"

Real Theta(Real omega, Real T)
{
  return SI::hbar * omega / (exp(SI::hbar * omega / (SI::kb * T)) - 1);
}

Real CalcHeatTransfer(Real a, Real eV, Real T1, Real T2, RWGTree *top, RWGTree *s, RWGTree *t, Options &options)
{
  Real omegaSI = eV * SI::eV / SI::hbar;
  Real omega = omegaSI * a / SI::c;
  
  GreensFunction greenExt(1.,1.,omega,true);
  
  s->diel->setTemperature(T1);
  t->diel->setTemperature(T2);
  
  Complex eps1 = s->diel->eps(eV);
  Complex mu1 = s->diel->mu(eV);
  Complex eps2 = t->diel->eps(eV);
  Complex mu2 = t->diel->mu(eV);
  
  s->g->init(eps1,mu1,omega);
  t->g->init(eps2,mu2,omega);
  
  ComplexDenseMatrix *W = Zts(top,top,&greenExt,options);
  
  ComplexDenseMatrix *G1 = Zts(s,s,s->g,options);
  ComplexDenseMatrix *G2 = Zts(t,t,t->g,options);
  for(int i=0; i<2*s->rwgs;i++)
    for(int j=0; j<2*s->rwgs;j++)
      (*W)(i,j) += (*G1)(i,j);
  
  int i0 = 2*s->rwgs;
  for(int i=0; i<2*t->rwgs;i++)
    for(int j=0; j<2*t->rwgs;j++)
      (*W)(i0+i,i0+j) += (*G2)(i,j);
  
  inverse(*W);
  ComplexDenseMatrix W1(*W,0,0,2*s->rwgs,2*(s->rwgs+t->rwgs));
  
  symm(*G1);
  ComplexDenseMatrix B1 = *G1 * W1;
  B1.resize(0,2*s->rwgs,2*s->rwgs,2*t->rwgs);
  delete G1;
  
  ComplexDenseMatrix W2(*W,0,2*s->rwgs,2*(s->rwgs+t->rwgs),2*t->rwgs);
  symm(*G2);
  ComplexDenseMatrix B2 = W2 * *G2;
  B2.resize(0,0,2*s->rwgs,2*t->rwgs);
  delete G2;
  
  multiply_Adagger_B(B2,B1,*W);
  
  Complex S = 0.25*trace(*W);
  
  delete W;
  return real(S) * (Theta(omega, T1) - Theta(omega, T2));
}

void parseShapes(RWGTree *top, RWGTree *s, RWGTree *t, char *body1x3d, char *body2x3d, const Vector3 &gap, Options &options) {
  
  ShapeParser parser;    
  Shape body1Shape;
  Dielectric *diel1;
  parser.parse(body1x3d,body1Shape,&diel1);
  Transform tform;
  tform.setT(gap);
  body1Shape.transform(tform);
  
  Shape body2Shape;
  Dielectric *diel2;
  parser.parse(body2x3d,body2Shape,&diel2);
  
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
  
  ShapeMesh refinedBody2Mesh(body2Shape,
                             options.sharpThresh,
                             options.concaveThresh,
                             options.relativeCostAngleConcave,
                             &refineEdgeFuncBody,
                             &maxEdgeLength);
  

  ShapeMeshDielectricPairList meshDiel;
  PointEdgeVectorMap eap;
  getEdgesAtPoints(refinedBody1Mesh,eap);
  getEdgesAtPoints(refinedBody2Mesh,eap);
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody1Mesh,diel1));
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody2Mesh,diel2));
  Real minSegmentRadius = 0.0;
  RWGTree rwgTree(meshDiel,eap,options.minSegmentSize,minSegmentRadius);
  s = rwgTree.sub[0];
  t = rwgTree.sub[1];    
  top = &rwgTree;
}
