#include "SurfaceParameterization.h"
#include "ConjugateGradient.h"
#include "ZRotation.h"

#include <queue>
using namespace std;

TriangleQuadTree :: TriangleQuadTree(TriangleSet &triangle, PointPoint2Map &u) 
{
  init(triangle, u, 0., 0., 2., 2.);
}

TriangleQuadTree :: TriangleQuadTree(TriangleSet &triangle, PointPoint2Map &u, Real x0, Real y0, Real w, Real h)
{
  init(triangle,u,x0,y0,w,h);
}

inline bool inside(Point2 &p, Real x0, Real y0, Real x1, Real y1)
{
  return (p.x >= x0 && p.x < x1) && (p.y >= y0 && p.y < y1);
}

void TriangleQuadTree :: init(TriangleSet &triangle, PointPoint2Map &u, Real x0, Real y0, Real w, Real h)
{
  TriangleSet triangle00;
  TriangleSet triangle01;
  TriangleSet triangle10;
  TriangleSet triangle11;

  this->x0 = x0;
  this->y0 = y0;

  if(triangle.size() == 1)
    t = *(triangle.begin());
  else {
    t = NULL;
    
    // split into subsets
    for(TriangleSet::const_iterator i = triangle.begin(); i != triangle.end(); ++i) {
      Triangle *t1 = *i;

      Point2 p1 = u[t1->p1];
      Point2 p2 = u[t1->p2];
      Point2 p3 = u[t1->p3];
      Point2 pc = p1; pc += p2; pc += p3;
      pc *= ONETHIRD;

      if(inside(pc,x0-w,y0-h,x0,y0))
        triangle00.insert(t1);
      else if(inside(pc,x0-w,y0,x0,y0+h))
        triangle01.insert(t1);
      else if(inside(pc,x0,y0-h,x0+w,y0))
        triangle10.insert(t1);
      else if(inside(pc,x0,y0,x0+w,y0+h))
        triangle11.insert(t1);
    }
  }

  w *= 0.5;
  h *= 0.5;

  // create subtrees
  if(triangle00.empty())
    sub00 = NULL;
  else
    sub00 = new TriangleQuadTree(triangle00, u, x0 - w, y0 - w, w, h);

  if(triangle01.empty())
    sub01 = NULL;
  else
    sub01 = new TriangleQuadTree(triangle01, u, x0 - w, y0 + w, w, h);

  if(triangle10.empty())
    sub10 = NULL;
  else
    sub10 = new TriangleQuadTree(triangle10, u, x0 + w, y0 - w, w, h);

  if(triangle11.empty())
    sub11 = NULL;
  else
    sub11 = new TriangleQuadTree(triangle11, u, x0 + w, y0 + w, w, h);

}

TriangleQuadTree :: ~TriangleQuadTree()
{
  if(sub00) delete sub00;
  if(sub01) delete sub01;
  if(sub10) delete sub10;
  if(sub11) delete sub11;
}

Triangle* TriangleQuadTree :: find(Real x, Real y)
{
  //printf("finding %g %g %g %g\n",x,y,x0,y0);
  if(t) {
    //printf("found\n");
    return t;
  }
  if(x < x0) {
    if(y < y0) {
      if(sub00) return sub00->find(x,y);
    } else {
      if(sub01) return sub01->find(x,y);
    }
  } else {
    if(y < y0) {
      if(sub10) return sub10->find(x,y);
    } else {
      if(sub11) return sub11->find(x,y);
    }
  }

  // not found yet, find first one in available subtrees
  if(sub00) return sub00->find(x,y);
  if(sub01) return sub01->find(x,y);
  if(sub10) return sub10->find(x,y);
  if(sub11) return sub11->find(x,y);
  
  return NULL;
}

Real cot(Vector3 &v1, Vector3 &v2)
{
  Vector3 z = v1 * v2;
  z.normalize();
  ZRotation R(z);
  Real x1, y1;
  R.rotate(v1,&x1,&y1);
  Real x2, y2;
  R.rotate(v2,&x2,&y2);
  Real normv1i = 1. / sqrt(Square(x1) + Square(y1));
  Real cosphi = normv1i * x1;
  Real sinphi = normv1i * y1;
  Real x0 = cosphi * x2 + sinphi * y2;
  Real y0 = -sinphi * x2 + cosphi * y2;
  return x0 / fabs(y0);
}

Real diam2(PointSet &point)
{
  Real x0 = 0.;
  Real x1 = 0.;
  Real y0 = 0.;
  Real y1 = 0.;
  Real z0 = 0.;
  Real z1 = 0.;
  bool bInit = false;
  for(PointSet::const_iterator i = point.begin(); i != point.end(); ++i) {
    Point *p = *i;
    if(!bInit || p->x < x0) x0 = p->x;
    if(!bInit || p->x > x1) x1 = p->x;
    if(!bInit || p->y < y0) y0 = p->y;
    if(!bInit || p->y > y1) y1 = p->y;
    if(!bInit || p->z < y0) z0 = p->z;
    if(!bInit || p->z > y1) z1 = p->z;
    bInit = true;
  }
  if(!bInit) return 0.;
  return Square(x1-x0) + Square(y1-y0) + Square(z1-z0);
}

class SparseMatrixVectorPoint2Multiply : public MatrixVectorMultiply<Point2> {
public:
  SparseMatrix<Real> &A;
  SparseMatrixVectorPoint2Multiply(SparseMatrix<Real> &A_) : A(A_) {}
  void multiply(Array<Point2> &x, Array<Point2> &y) {
    y.zero();
    A.MatrixVectorMultiply((Real*)x.p,2,(Real*)y.p,2);
    A.MatrixVectorMultiply(((Real*)x.p)+1,2,((Real*)y.p)+1,2);
  }
};

class DistortionArgs {
public:
  DistortionArgs(TriangleRealMap &triangleArea_,
                 PointPoint2Map &uA_,
                 PointPoint2Map &uchi_)
    : triangleArea(triangleArea_),
      uA(uA_),
      uchi(uchi_) {}

  TriangleRealMap &triangleArea;
  PointPoint2Map &uA;
  PointPoint2Map &uchi;
};

Real distortion(Real lambda,
                void *args_)
{
  DistortionArgs *args = (DistortionArgs*)args_;
  TriangleRealMap &triangleArea = args->triangleArea;
  PointPoint2Map &uA = args->uA;
  PointPoint2Map &uchi = args->uchi;

  Real distortion = 0.;
  Real mu = 1. - lambda;
  Real area3sum = 0.;
  for(TriangleRealMap::const_iterator i = triangleArea.begin(); i != triangleArea.end(); ++i) {
    Real area3 = i->second;
    area3sum += area3;
  }

  // area of unit hexagon
  Real area2sum = THREERHALVESSQRT3;

  for(TriangleRealMap::const_iterator i = triangleArea.begin(); i != triangleArea.end(); ++i) {
    Triangle *t = i->first;    
    Real area3 = i->second;

    Point2 p1A = lambda * uA[t->p1];
    Point2 p1chi = mu * uchi[t->p1];
    Point2 p1 = p1A + p1chi; 

    Point2 p2A = lambda * uA[t->p2];
    Point2 p2chi = mu * uchi[t->p2];
    Point2 p2 = p2A + p2chi; 

    Point2 p3A = lambda * uA[t->p3];
    Point2 p3chi = mu * uchi[t->p3];
    Point2 p3 = p3A + p3chi; 

    Point2 v2 = p2 - p1;
    Point2 v3 = p3 - p1;
    Real area2 = 0.5 * fabs(v2.x * v3.y - v3.x * v2.y);
    
    Real dist = Square((area2 * area3sum) / (area3 * area2sum) - 1.);
    distortion += dist;
  }
  return distortion;
}

SurfaceParameterization :: SurfaceParameterization(PointSet &point, 
                                                   PointSet &boundary,
                                                   PointEdgeVectorMap &edgesAtPoint,
                                                   PointEdgeVectorMap &boundaryEdgesAtPoint,
                                                   TriangleSet &triangle_)
  : triangle(triangle_)
{
  // enumerate all closed boundary curves
  PointSetSet boundariesPoints;

  for(PointSet::const_iterator i = boundary.begin(); i != boundary.end(); ++i) {
    Point *p1 = *i;
    bool bAssigned = false;
    for(PointSetSet::const_iterator j = boundariesPoints.begin(); j != boundariesPoints.end(); ++j) {
      PointSet *boundaryPoints = *j;
      if(boundaryPoints->find(p1) != boundaryPoints->end()) {
        bAssigned = true;
        break;
      }
    }
    if(bAssigned) continue;

    Point *p = p1;
    Point *pprev = NULL;
    PointSet *boundaryPoints;
    stack<PointSet*> boundaryPointsStack;
    stack<Point*> pStack;
    stack<Point*> pprevStack;
    bool bMore;

    boundaryPoints = new PointSet;

    do {
      bool bClosed = false;
      while(p) {
        boundaryPoints->insert(p);
        Point *pnext = NULL;
        vector<Edge*> &eap = boundaryEdgesAtPoint[p];
        int edges = 0;
        bClosed = false;
        for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
          Edge *e = *j;
          Point *p2 = (e->p1 == p) ? e->p2 : e->p1;

          // only follow the boundary
          if(boundary.find(p2) == boundary.end()) continue;

          // the boundary is closed if there is a boundary edge connecting to the starting point
          if(p2 == p1 && pprev != p1) bClosed = true;          

          // don't retraverse points
          if(boundaryPoints->find(p2) != boundaryPoints->end()) continue;

          // if there is more than one path, put everything in stacks
          if(edges) {
            PointSet *pushBoundaryPoints = new PointSet;
            pushBoundaryPoints->insert(boundaryPoints->begin(),
                                       boundaryPoints->end());
            boundaryPointsStack.push(pushBoundaryPoints);
            pprevStack.push(p);
            pStack.push(p2);
          } else {
            pnext = p2;
          }
          edges++;
        }
        pprev = p;
        p = pnext;
      }

      // add this boundary if it is closed
      if(bClosed) {
        boundariesPoints.insert(boundaryPoints);
      } else {
        delete boundaryPoints;
      }

      // follow the next path in the stack
      if(!pStack.empty()) {
        p = pStack.top(); pStack.pop();
        pprev = pprevStack.top(); pprevStack.pop();
        boundaryPoints = boundaryPointsStack.top(); boundaryPointsStack.pop();
        bMore = true;
      } else {
        bMore = false;
      }
    } while(bMore);
  }

  // now find the boundary with the largest diameter
  PointSet *outerBoundaryPoints = NULL;
  if(boundariesPoints.size() == 1)
    outerBoundaryPoints = *(boundariesPoints.begin());
  else {
    Real dmax = -1.;
    for(PointSetSet::const_iterator j = boundariesPoints.begin(); j != boundariesPoints.end(); ++j) {
      PointSet *boundaryPoints = *j;
      Real d = diam2(*boundaryPoints);
      if(dmax == -1. || d > dmax) {
        outerBoundaryPoints = boundaryPoints;
        dmax = d;
      }
    }
  }
  
  Point *pc0, *pc1, *pc2, *pc3, *pc4, *pc5;
  pc0 = NULL;
  Real dmax = -1.;
  for(PointSet::const_iterator i = outerBoundaryPoints->begin(); i != outerBoundaryPoints->end(); ++i) {
    Point *pi = *i;
    for(PointSet::const_iterator j = outerBoundaryPoints->begin(); j != outerBoundaryPoints->end(); ++j) {
      Point *pj = *j;
      Vector3 v = *pi - *pj;
      Real d = norm(v);
      if(d > dmax) {
        dmax = d;
        pc0 = pi;
      }
    }
  }

  // convert the outer boundary point set to a vector of points
  EdgeVector outerBoundaryEdges;
  Point *p = pc0;
  Point *pprev = NULL;
  Real l012345 = 0.;
  while(p) {
    vector<Edge*> &eap = boundaryEdgesAtPoint[p];
    Point *pnext = NULL;
    for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
      Edge *e = *j;
      Point *p2 = (e->p1 == p) ? e->p2 : e->p1;
      if(p2 == pprev) continue;
      if(outerBoundaryPoints->find(p2) == outerBoundaryPoints->end()) continue;
      pnext = p2;
      outerBoundaryEdges.push_back(e);
      l012345 += e->length;
      break;
    }
    pprev = p;
    p = pnext;
    if(p == pc0)
      break;
  }

  // the point that most closely bisects the length of the boundary will be corner 3
  // find the length of the sides 012 and 345 of the parameterization square
  pc3 = NULL;
  bool bSide345 = false;
  Real l012 = 0.;
  Real l345 = 0.;
  Real l = 0.;
  p = pc0;
  pprev = NULL;
  Real dprev = -1.;
  for(EdgeVector::const_iterator j = outerBoundaryEdges.begin(); j != outerBoundaryEdges.end(); ++j) {
    Edge *e = *j;
    Point *pnext = (e->p1 == p) ? e->p2 : e->p1;
    l = e->length;
    pprev = p;
    p = pnext;
    Real d = l012 + l - 0.5*l012345;
    if(!bSide345 && dprev <= 0. && d > 0.) {
      if(fabs(dprev) < fabs(d)) {
        pc3 = pprev;
        l345 += l;
      } else {
        pc3 = p;
        l012 += l;
      }
      bSide345 = true;
    } else {
      if(bSide345) 
        l345 += l;
      else
        l012 += l;
    }
    dprev = d;
  }

  // choose the other corners
  // and find the length of each side
  p = pc0;
  pprev = NULL;
  pc1 = NULL;
  pc2 = NULL;
  pc4 = NULL;
  pc5 = NULL;
  int side = 0;
  Real l0 = 0.;
  Real l1 = 0.;
  Real l2 = 0.;
  Real l3 = 0.;
  Real l4 = 0.;
  Real l5 = 0.;
  Real dmin = -1.;
  for(EdgeVector::const_iterator j = outerBoundaryEdges.begin(); j != outerBoundaryEdges.end(); ++j) {
    Edge *e = *j;
    Point *pnext = (e->p1 == p) ? e->p2 : e->p1;
    l = e->length;
    pprev = p;
    p = pnext;

    Real d;
    if(side == 0) {
      d = fabs(l0 + l - ONETHIRD * l012);
      if(dmin == -1. || d < dmin) {
        l0 += l;
        dmin = d;
      } else {
        pc1 = pprev;
        l1 = l;
        side = 1;
        dmin = fabs(l0 + l1 - TWOTHIRDS * l012);
      }
    } else if(side == 1) {
      d = fabs(l0 + l1 + l - TWOTHIRDS * l012);
      if(d < dmin) {
        l1 += l;
        dmin = d;
      } else {
        pc2 = pprev;
        l2 = l;
        side = 2;
      }
    } else if(side == 2) {
      if(pprev == pc3) {
        side = 3;
        l3 = l;
        dmin = fabs(l3 - ONETHIRD * l345);
      } else {
        l2 += l;
      }
    } else if(side == 3) {
      d = fabs(l3 + l - ONETHIRD * l345);
      if(d < dmin) {
        l3 += l;
        dmin = d;
      } else {
        pc4 = pprev;
        l4 = l;
        side = 4;
        dmin = fabs(l3 + l4 - TWOTHIRDS * l345);
      }
    } else if(side == 4) {
      d = fabs(l3 + l4 + l - TWOTHIRDS * l345);
      if(d < dmin) {
        l4 += l;
        dmin = d;
      } else {
        pc5 = pprev;
        l5 = l;
        side = 5;
        dmin = -1.;
      }
    } else if(side == 5) {
      l5 += l;
    }
  }

  if(!l0 || !l1 || !l2 || !l3 || !l4 || !l5)
    abort();

  //cout << pc0 << " " << pc1 << " " << pc2 << " " << pc3 << " " << pc4 << " " << pc5 << "\n";
  //printf("%g %g / %g %g %g %g %g %g\n",l012,l345,l0,l1,l2,l3,l4,l5);

  // assign positions to outer boundary
  p = pc0;
  pprev = NULL;
  u[p] = Point2(-1.,0.);
  Real ls = 0.;
  Real lt = l0;
  side = 0;
  for(EdgeVector::const_iterator j = outerBoundaryEdges.begin(); j != outerBoundaryEdges.end(); ++j) {
    Edge *e = *j;
    Point *pnext = (e->p1 == p) ? e->p2 : e->p1;
    l = e->length;
    pprev = p;
    p = pnext;
    
    if(p == pc0) {
      break;
    } else if(p == pc1) {
      u[p] = Point2(-0.5,-SQRT3OVER2);
      ls = 0.;
      lt = l1;
      side = 1;
    } else if(p == pc2) {
      u[p] = Point2(0.5,-SQRT3OVER2);
      ls = 0.;
      lt = l2;
      side = 2;
    } else if(p == pc3) {
      u[p] = Point2(1.,0);
      ls = 0.;
      lt = l3;
      side = 3;
    } else if(p == pc4) {
      u[p] = Point2(0.5,SQRT3OVER2);
      ls = 0.;
      lt = l4;
      side = 4;
    } else if(p == pc5) {
      u[p] = Point2(-0.5,SQRT3OVER2);
      ls = 0.;
      lt = l5;
      side = 5;
    } else {    
      ls += l;
      Real a = ls/lt;
      if(side == 0) {
        u[p] = Point2(-1. + 0.5 * a, -SQRT3OVER2 * a);
      } else if(side == 1) {
        u[p] = Point2(-0.5 + a, -SQRT3OVER2);
      } else if(side == 2) {
        u[p] = Point2(0.5 + 0.5 * a, -SQRT3OVER2 + SQRT3OVER2 * a);
      } else if(side == 3) {
        u[p] = Point2(1. - 0.5 * a, SQRT3OVER2 * a);
      } else if(side == 4) {
        u[p] = Point2(0.5 - a, SQRT3OVER2);
      } else if(side == 5) {
        u[p] = Point2(-0.5 - 0.5 * a, SQRT3OVER2 - SQRT3OVER2 * a);
      } 
    }
  }

  // index the points
  PointIntMap index;
  int n = 0;
  for(PointSet::const_iterator i = point.begin(); i != point.end(); ++i) {
    Point *p = *i;
    index[p] = n++;
  }

  // M*U=C  
  Array<Point2> UA(n);
  Array<Point2> Uchi(n);
  Array<Point2> C(n);
  SparseMatrix<Real> MA(n,n);
  SparseMatrix<Real> Mchi(n,n);
  Array<Real> diagA(n);
  Array<Real> diagchi(n);
  PointPoint2Map uA;
  PointPoint2Map uchi;

  for(PointSet::const_iterator i = point.begin(); i != point.end(); ++i) {
    Point *p1 = *i;
    int ii = index[p1];
    if(u.find(p1) != u.end()) {
      MA.insert(ii,ii,1.);
      Mchi.insert(ii,ii,1.);
      diagA[ii] = 1.;
      diagchi[ii] = 1.;
      UA[ii] = u[p1];
      Uchi[ii] = u[p1];
      C[ii] = u[p1];
    } else {
      Real EAsum = 0.;
      Real Echisum = 0.;
      vector<Edge*> &eap = edgesAtPoint[p1];
      for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
        Edge *e = *j;
        Point *p2 = (e->p1 == p1) ? e->p2 : e->p1;
        if(point.find(p2) == point.end()) continue;
        int jj = index[p2];

        Vector3 ve = *p1 - *p2;
        
        Triangle *t1 = e->t1;
        Point *p31;
        if(t1->p1 != p1 && t1->p1 != p2)
          p31 = t1->p1;
        else if(t1->p2 != p1 && t1->p2 != p2)
          p31 = t1->p2;
        else
          p31 = t1->p3;
        Vector3 v31 = *p1 - *p31;
        Vector3 v32 = *p2 - *p31;
        
        Real cotalpha = cot(v31,v32);
        v32 = *p31 - *p2;
        Real cotdelta = cot(v32,ve);
        
        Triangle *t2 = e->t2;
        Point *p32;
        if(t2->p1 != p1 && t2->p1 != p2)
          p32 = t2->p1;
        else if(t2->p2 != p1 && t2->p2 != p2)
          p32 = t2->p2;
        else
          p32 = t2->p3;
        v31 = *p1 - *p32;
        v32 = *p2 - *p32;
        
        Real cotbeta = cot(v32,v31);        
        v32 = *p32 - *p2;
        Real cotgamma = cot(ve,v32);
        
        // energies
        Real EA = cotalpha + cotbeta;
        Real Echi = (cotdelta + cotgamma) / Square(e->length);
        EAsum -= EA;
        Echisum -= Echi;
        //printf("%g %g %g %g\n",cotalpha,cotbeta,cotdelta,cotgamma);
        // matrix elements
        MA.insert(ii,jj,EA);
        Mchi.insert(ii,jj,Echi);
      }
      MA.insert(ii,ii,EAsum);
      Mchi.insert(ii,ii,Echisum);
      diagA[ii] = EAsum;
      diagchi[ii] = Echisum;
      UA[ii] = Point2(0.,0.);
      Uchi[ii] = Point2(0.,0.);
      C[ii] = Point2(0.,0.);
    }
  }
  MA.end();
  Mchi.end();

#define CONJGRADEPS 1e-7

  DiagonalPreconditioner<Real,Point2> PA(diagA);
  SparseMatrixVectorPoint2Multiply MAmvm(MA);
  conjugateGradientSquared<Real,Real,Point2>(MAmvm,PA,UA,C,CONJGRADEPS,1024);

  DiagonalPreconditioner<Real,Point2> Pchi(diagchi);
  SparseMatrixVectorPoint2Multiply Mchimvm(Mchi);
  conjugateGradientSquared<Real,Real,Point2>(Mchimvm,Pchi,Uchi,C,CONJGRADEPS,1024);

#undef CONJGRADEPS

  n = 0;
  for(PointSet::const_iterator i = point.begin(); i != point.end(); ++i) {
    Point *p = *i;
    uA[p] = UA[n];
    uchi[p] = Uchi[n];
    n++;
  }

  // Map triangles to area
  TriangleRealMap triangleArea;
  for(TriangleSet::const_iterator i = triangle.begin(); i != triangle.end(); ++i) {
    Triangle *t = *i;
    triangleArea[t] = area(t);
  }

  // minimize area distortion
  DistortionArgs args(triangleArea,
                      uA,
                      uchi);
  Real lambda;
  Real a = -100.;
  Real b = 0;
  Real c = 100.;
  Real fa,fb,fc;
  mnbrak(&a,&b,&c,&fa,&fb,&fc,&distortion,&args);
  brent(a,b,c,&distortion,&args,1e-2,&lambda);

  // the optimal Map
  for(PointSet::const_iterator i = point.begin(); i != point.end(); ++i) {
    Point *p1 = *i;
    Point2 p1A = lambda * uA[p1];
    Point2 p1chi = (1. - lambda) * uchi[p1];
    u[p1] = p1A + p1chi;
  }

  // delete
  for(PointSetSet::const_iterator j = boundariesPoints.begin(); j != boundariesPoints.end(); ++j) {
    PointSet *boundaryPoints = *j;
    delete boundaryPoints;
  }
  
  // create quad tree
  quad = new TriangleQuadTree(triangle, u);
}

SurfaceParameterization :: ~SurfaceParameterization()
{
  delete quad;
}

Point2& SurfaceParameterization :: operator[](Point *p)
{
  return u[p];
}

Point SurfaceParameterization :: inverse(Point2 &p2)
{
  Triangle *t = quad->find(p2.x,p2.y);
  if(!t) {
    cout << "surface parameterization quadtree doesn't contain point\n";
    abort();
  }

  // the returned triangle may not actually contain the point
  // so check the neighboring triangles
  queue<Triangle*> search;
  TriangleSet searched;
  search.push(t);
  searched.insert(t);
  Real x = -1.;
  Real y = -1.;

  Point *p31;
  Point *p32;
  Point *p33;

  bool bFound = false;
  while(!search.empty()) {

    t = search.front();
    search.pop();

    // 3D points
    p31 = t->p1;
    p32 = t->p2;
    p33 = t->p3;
    
    // 2D points
    Point2 p21 = u[p31];
    Point2 p22 = u[p32];
    Point2 p23 = u[p33];

    // 2D vectors
    Vector2 v210 = p2 - p21; 
    Vector2 v212 = p22 - p21;
    Vector2 v213 = p23 - p21;
    
    // rotate
    Real normv212i = 1. / sqrt(norm(v212));
    Real cosphi = normv212i * v212.x;
    Real sinphi = normv212i * v212.y;
    Real v3x = normv212i * (cosphi * v213.x + sinphi * v213.y);
    Real v3y = normv212i * (-sinphi * v213.x + cosphi * v213.y);
    Real x0 = normv212i * (cosphi * v210.x + sinphi * v210.y);
    Real y0 = normv212i * (-sinphi * v210.x + cosphi * v210.y);

    // triangle coordinates
    y = y0 / v3y;
    x = x0 - v3x * y;
    
    //printf("searching for interior %g %g\n",x,y);

    if(x > -1e-3 && y > -1e-3 && 1. - x - y > -1e-3) {
      bFound = true;
      //printf("found interior %g %g\n",x,y);
      break;
    }

    Triangle *t1;

    t1 = (t->e12->t1 == t) ? t->e12->t2 : t->e12->t1;
    if(triangle.find(t1) != triangle.end() && searched.find(t1) == searched.end()) {
      search.push(t1);
      searched.insert(t1);
    }

    t1 = (t->e23->t1 == t) ? t->e23->t2 : t->e23->t1;
    if(triangle.find(t1) != triangle.end() && searched.find(t1) == searched.end()) {
      search.push(t1);
      searched.insert(t1);
    }

    t1 = (t->e31->t1 == t) ? t->e31->t2 : t->e31->t1;
    if(triangle.find(t1) != triangle.end() && searched.find(t1) == searched.end()) {
      search.push(t1);
      searched.insert(t1);
    }
  }

  if(!bFound) {
    cout << "surface parameterization couldn't find a triangle that contains point\n";
    abort();
  }

  // 3D vectors
  Vector3 v312 = *p32 - *p31;
  Vector3 v313 = *p33 - *p31;
  
  // 3D triangle coordinates
  Vector3 p3x = x * v312;
  Vector3 p3y = y * v313;
  
  // 3D position
  Point p3 = *p31; p3 += p3x; p3 += p3y;
    
  return p3;
}

int SurfaceParameterization :: hexagonalPoints(int n)
{
  return 1 + 3 * n * (n - 1);
}

Point2 SurfaceParameterization :: hexagonalPoint(int n, int p, int q)
{
  if(n == 1) {
    Point2 p2(0.,0.);
    return p2;
  }   
  Real a = 1. / (Real)(n-1);
  Point2 p2(a * (p + 0.5 * q), a * SQRT3OVER2 * q);
  return p2;
}

ostream& operator<<(ostream& os, Point2 &p2) {
  os << "[" << p2.x << " " << p2.y << "]";
  return os;
}

ostream& operator<<(ostream& os, SurfaceParameterization &u)
{
  os << "s = [\n";
  for(TriangleSet::const_iterator i = u.triangle.begin(); i != u.triangle.end(); ++i) {
    Triangle *t = *i;
    os << "[" << u[t->p1] << " " << u[t->p2] << " " << u[t->p3] << "];\n";
  }
  os << "];\n";
  return os;
}
