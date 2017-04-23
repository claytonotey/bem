#include "Util.h"
#include "ShapeMesh.h"
#include "ZRotation.h"
#include "MathUtils.h"

using namespace std;
#include <vector>
#include <stack>

Real sharpness(Edge *e)
{
  Vector3 n1 = normal(e->t1);
  Vector3 n2 = normal(e->t2);
  Vector3 z = n1 * n2;
  Real normz = sqrt(norm(z));
  if(normz / sqrt(norm(n1)) < 1e-6)
    return 0.;
  else {
    z *= 1. / normz;
    Real x = dot(n2,n1);
    Real y = dot(n2,z * n1);
    return atan2(y,x);
  }
}

Real angle(Edge *e1, Edge *e2)
{
  Vector3 v1,v2;
  if(e1->p1 == e2->p1) {
    v1 = *(e1->p2) - *(e1->p1);
    v2 = *(e2->p2) - *(e2->p1);
  } else if(e1->p1 == e2->p2) {
    v1 = *(e1->p2) - *(e1->p1);
    v2 = *(e2->p1) - *(e2->p2);
  } else if(e1->p2 == e2->p1) {
    v1 = *(e1->p1) - *(e1->p2);
    v2 = *(e2->p2) - *(e2->p1);
  } else if(e1->p2 == e2->p2) {
    v1 = *(e1->p1) - *(e1->p2);
    v2 = *(e2->p1) - *(e2->p2);
  } else {
    abort();
  }
  Vector3 z = v1 * v2;
  Real normz = sqrt(norm(z));
  if(normz / sqrt(norm(v1)) < 1e-6)
    return PI;
  z *= 1. / normz;
  Real x = dot(v2,v1);
  Real y = dot(v2,z * v1);
  return atan2(y,x);
}

Point *sharedPoint(Edge *e1, Edge *e2)
{
  if(e1->p1 == e2->p1)
    return e1->p1;
  if(e1->p2 == e2->p1)
    return e1->p2;
  if(e1->p1 == e2->p2)
    return e1->p1;
  if(e1->p2 == e2->p2)
    return e1->p2;
  return NULL;
}

#define DISTANCE_INVALID 1e99

// attempt a cut using the chosen rotation angle theta, and return the distance
Real circumCut(Real theta,
               Point *p1, 
               Point *p2,
               Vector3 &direction,
               ZRotation &ZR,
               Point &O,
               PointEdgeVectorMap &edgesAtPoint,
               PointSet &boundaryPoints,
               Cut *c = NULL,
               bool bCrossingOK = false)
{                        
  Real costheta = cos(theta);
  Real sintheta = sin(theta);

  PointSet pointSet;

  Real distance = 0.;
  Point *p = p1;
  while(p != p2) {
    vector<Edge*> &eap = edgesAtPoint[p];
    Real minDist = -1.;                         
    Edge *bestEdge = NULL;
    Point *bestPoint = NULL;
    for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
      Edge *e = *j;
      Vector3 v;
      Point *popp;
      if(e->p1 == p) {
        v = *(e->p2) - *(e->p1);
        popp = e->p2;
      } else {
        v = *(e->p1) - *(e->p2);
        popp = e->p1;
      }
      // optionally don't cross other boundaries
      if(bCrossingOK || (boundaryPoints.find(popp) == boundaryPoints.end())) {
        // don't go backwards (too much)
        Real forward = dot(direction,v);
        forward = forward>=0. ? 0. : -Square(forward) / norm(v);
        if(forward >= -.01) {
          // choose the point closest to the chosen cross section boundary
          Vector3 Op = *popp - O;
          Vector3 ZRp = ZR.rotate(Op);
          Real y = -sintheta * ZRp.x + costheta * ZRp.y;
          if(minDist == -1. || y < minDist) {
            bestPoint = popp;
            bestEdge = e;
            minDist = y;
          }        
        }
      }
    }
    // never self-cross
    if(bestPoint == p1 || pointSet.find(bestPoint) != pointSet.end())
      return DISTANCE_INVALID;
    if(bestPoint == NULL) {
      return DISTANCE_INVALID;
    } else {
      distance += bestEdge->length;
      pointSet.insert(bestPoint);
      if(c) {
        if(bestPoint == p2) {
          c->edge.push_back(bestEdge);
        } else {
          c->point.push_back(bestPoint);
          c->edge.push_back(bestEdge);
        }
      }
    }
    direction = *bestPoint - *p;
    direction.normalize();
    p = bestPoint;
  }
  return distance;
}

class DistanceArgs {
public:
  DistanceArgs(Point *p1, 
               Point *p2,
               Vector3 &direction_,
               ZRotation &ZR_,
               Point &O_, 
               PointEdgeVectorMap &edgesAtPoint_,
               PointSet &boundaryPoints_) 
    : direction(direction_),
      ZR(ZR_),
      O(O_), 
      edgesAtPoint(edgesAtPoint_), 
      boundaryPoints(boundaryPoints_)
  {
    this->p1 = p1;
    this->p2 = p2;
  }

  Point *p1;
  Point *p2;
  Vector3 &direction;
  ZRotation &ZR;
  Point &O;
  PointEdgeVectorMap &edgesAtPoint;
  PointSet &boundaryPoints;
};

Real cutDistance(Real theta, void *args)
{
  DistanceArgs *dargs = (DistanceArgs*) args;
  return circumCut(theta,
                   dargs->p1, 
                   dargs->p2,
                   dargs->direction,
                   dargs->ZR,
                   dargs->O, 
                   dargs->edgesAtPoint,
                   dargs->boundaryPoints);
}

Real circumCut(Point *p1,
               Point *p2,
               PointEdgeVectorMap &edgesAtPoint,
               PointEdgeVectorMap &boundaryEdgesAtPoint,
               PointSet &boundaryPoints,
               Cut *c = NULL,
               bool bCrossingOK = false)
{
  Vector3 z = *p2 - *p1;

  // find the direction to travel
  Point *p = p1;
  vector<Edge*> &beap = boundaryEdgesAtPoint[p];
  vector<Edge*> &eap = edgesAtPoint[p];
  Vector3 direction;
  Real maxminAngle = 0.;
  for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
    Edge *ej = *j;
    Real minAngle = TWOPI;
    for(vector<Edge*>::const_iterator k = beap.begin(); k != beap.end(); ++k) {
      Edge *ek = *k;
      if(ej == ek) continue;
      Real a = fabs(angle(ej,ek));
      if(a < minAngle) {
        minAngle = a;
      }
    }
    if(minAngle > maxminAngle) {
      maxminAngle = minAngle;
      if(ej->p1 == p)
        direction = *(ej->p2) - *(ej->p1);
      else
        direction = *(ej->p1) - *(ej->p2);
    }
  }
  if(maxminAngle == 0.)
    return DISTANCE_INVALID;

  direction.normalize();
  
  // the origin for the rotation
  Point O = (*p1) + 0.5 * z;
  Vector3 zn = z;
  zn.normalize();
  ZRotation ZR(zn);

  // find the rotation about the z axis that minimizes the length of cut
  Real theta;
  DistanceArgs distanceArgs(p1,p2,direction,ZR,O,edgesAtPoint,boundaryPoints);

  // bracket the minimum
  Real d0 = circumCut(0.,p1,p2,direction,ZR,O,edgesAtPoint,boundaryPoints);
  Real dPIOVER2 = circumCut(PIOVER2,p1,p2,direction,ZR,O,edgesAtPoint,boundaryPoints);
  if(d0 <= dPIOVER2)
    brent(-PIOVER2,0.,PIOVER2,&cutDistance,&distanceArgs,1e-3,&theta);
  else
    brent(0.,PIOVER2,PI,&cutDistance,&distanceArgs,1e-3,&theta);

  // cut using the optimal theta
  Real distance = circumCut(theta,p1,p2,direction,ZR,O,edgesAtPoint,boundaryPoints,c,bCrossingOK);
  if(distance == DISTANCE_INVALID)
    return distance;

  c->p1 = p1;
  c->p2 = p2;
  c->e1.insert(c->e1.begin(),boundaryEdgesAtPoint[p1].begin(),boundaryEdgesAtPoint[p1].end());
  c->e2.insert(c->e2.begin(),boundaryEdgesAtPoint[p2].begin(),boundaryEdgesAtPoint[p2].end());
  return distance;
}

// find the boundary points that are extrema of the distance function
 void findExtremalBoundaryPoints(PointEdgeVectorMap &boundaryEdgesAtPoint,
                                 PointSet &boundaryPoints, 
                                 PointSet &cutPoints, 
                                 Real (*distFunc)(Point*, void *),
                                 void *args,
                                 PointSet &extremalPoints)
{
  // Map boundary points to distances
  PointRealMap dist;
  for(PointSet::const_iterator j = boundaryPoints.begin(); j != boundaryPoints.end(); ++j) {
    Point *p = *j;
    dist[p] = distFunc(p,args);
  }

  // find the extrema in the distance from the point to the boundary points
  for(PointSet::const_iterator j = boundaryPoints.begin(); j != boundaryPoints.end(); ++j) {
    Point *p = *j;

    // don't include cut points
    if(cutPoints.find(p) != cutPoints.end()) continue;

    // don't include the point itself
    Real distj = dist[p];
    if(distj == 0.) continue;

    // don't include endpoints
    vector<Edge*> &beap = boundaryEdgesAtPoint[p];
    int size = beap.size();
    if(size < 2) continue;
    
    bool bMin = true;
    bool bMax = true;
    for(vector<Edge*>::const_iterator k = beap.begin(); k != beap.end(); ++k) {
      Edge *e = *k;
      Real distk;
      if(e->p1 != p)
        distk = dist[e->p1];
      else
        distk = dist[e->p2];
      if(distk <= distj) 
        bMin = false;
      if(distk >= distj)
        bMax = false;
    }
    if(bMin || bMax) {
      extremalPoints.insert(p);
    }
  }
}

Real selectMinima(PointSet &extremalPoints,
                  PointEdgeVectorMap &boundaryEdgesAtPoint,
                  Real (*distFunc)(Point*, void *),
                  void *args,
                  Real d2Thresh,
                  PointSet &minimalPoints,
                  PointRealMap &secondDerivative)
{
  Real minDist = -1.;
  for(PointSet::const_iterator j = extremalPoints.begin(); j != extremalPoints.end(); ++j) {
    Point *p = *j;

    Real distj = distFunc(p,args);

    vector<Edge*> &beap = boundaryEdgesAtPoint[p];
    int size = beap.size();

    if(size < 2) continue;

    Edge *e1;
    Edge *e2;
    Real dist1, dist2;
    
    Edge *e10 = beap[0];
    Edge *e20 = beap[1];
    
    Point *p1;
    if(e10->p1 == p)
      p1 = e10->p2;
    else
      p1 = e10->p1;
    
    Point *p2;
    if(e20->p1 == p)
      p2 = e20->p2;
    else
      p2 = e20->p1;
    
    Real dist10 = distFunc(p1,args);
    Real dist20 = distFunc(p2,args);
    
    if(size == 2) {
      e1 = e10;
      e2 = e20;
      dist1 = dist10;
      dist2 = dist20;
    } else {
      if(dist10 < dist20) {
        dist1 = dist10;
        dist2 = dist20;
        e1 = e10;
        e2 = e20;
      } else {
        dist1 = dist20;
        dist2 = dist10;
        e1 = e20;
        e2 = e10;
      }
      for(int i=2; i<size; i++) {
        Edge *e3 = beap[i];
        Point *p3;
        if(e3->p1 == p)
          p3 = e3->p2;
        else
          p3 = e3->p1;
        Real dist3 = distFunc(p3,args);
        if(dist3 < dist1) {
          e2 = e1;
          e1 = e3;
          dist2 = dist1;
          dist1 = dist3;
        } else if(dist3 < dist2) {
          e2 = e3;
          dist2 = dist3;
        }
      }
    }

    if(dist1 <= distj || dist2 <= distj) continue;

    Real L1 = e1->length;
    Real L2 = e2->length;

    // poor man's 2nd derivative
    // first, normalize
    Real n = 2. / (L1 + L2);
    Real x0 = L1 * n;
    Real x1 = L2 * n;
    Real y0 = dist1 * n;
    Real y1 = distj * n;
    Real y2 = dist2 * n;
    Real d2dd2 = 2. * (x0 * (y2 - y1) + x1 * (y0 - y1)) / (x0 * x1 * (x0 + x1));
  
    // only include this point if the 2nd derivative exceeds a threshold
    if(d2dd2 >= d2Thresh) {
      if(minDist == -1. || distj < minDist)
        minDist = distj;
      minimalPoints.insert(p);
      secondDerivative[p] = d2dd2;
    }
  }
  return minDist;
}

Real distToPoint(Point *p, void *args)
{
  Point *p0 = (Point*)args;
  Vector3 v = *p - *p0;
  return sqrt(norm(v));
}

class MinDistToBoundaryArgs {
public:
  MinDistToBoundaryArgs(PointEdgeVectorMap &edgesAtPoint_,
                        PointEdgeVectorMap &boundaryEdgesAtPoint_,
                        PointSet &boundaryPoints_,
                        PointSet &cutPoints_) :
    edgesAtPoint(edgesAtPoint_),
    boundaryEdgesAtPoint(boundaryEdgesAtPoint_),
    boundaryPoints(boundaryPoints_),
    cutPoints(cutPoints_) {}

  Point *p;
  PointEdgeVectorMap &edgesAtPoint;
  PointEdgeVectorMap &boundaryEdgesAtPoint;
  PointSet &boundaryPoints;
  PointSet &cutPoints;
};

Real surfaceDistToPoint(Point *p, void *args)
{
  MinDistToBoundaryArgs *mdtbargs = (MinDistToBoundaryArgs*) args;
  return circumCut(mdtbargs->p,
                   p,
                   mdtbargs->edgesAtPoint,
                   mdtbargs->boundaryEdgesAtPoint,
                   mdtbargs->boundaryPoints);
}

Real minDistToBoundary(Point *p, void *args)
{
  MinDistToBoundaryArgs *mdtbargs = (MinDistToBoundaryArgs*) args;
  mdtbargs->p = p;
  PointSet extremalPoints;
  findExtremalBoundaryPoints(mdtbargs->boundaryEdgesAtPoint,
                             mdtbargs->boundaryPoints, 
                             mdtbargs->cutPoints, 
                             &distToPoint,
                             p,
                             extremalPoints);

  PointSet minimalPoints;
  PointRealMap secondDerivative;
  return selectMinima(extremalPoints,
                      mdtbargs->boundaryEdgesAtPoint,
                      &surfaceDistToPoint,
                      mdtbargs,
                      0.,
                      minimalPoints,
                      secondDerivative);
}

void findClosestBoundaryPoints(PointEdgeVectorMap &boundaryEdgesAtPoint,
                               PointSet &boundaryPoints,
                               PointSet &cutPoints,
                               Point *p,
                               PointSet &closestPoints)
{
  PointSet extremalPoints;
  findExtremalBoundaryPoints(boundaryEdgesAtPoint,
                             boundaryPoints,
                             cutPoints,
                             &distToPoint,
                             p,
                             extremalPoints);

  PointRealMap secondDerivative;
  selectMinima(extremalPoints,
               boundaryEdgesAtPoint,
               &distToPoint,
               p,
               0.,
               closestPoints,
               secondDerivative);
}

void findConcavePoints(PointEdgeVectorMap &edgesAtPoint,
                       PointEdgeVectorMap &boundaryEdgesAtPoint,
                       PointSet &boundaryPoints,
                       PointSet &cutPoints,
                       Real concaveThresh,
                       PointSet concavePoints,
                       PointRealMap concavity)
{
  PointSet extremalPoints;
  MinDistToBoundaryArgs mdtbargs(edgesAtPoint,boundaryEdgesAtPoint,cutPoints,extremalPoints);

  findExtremalBoundaryPoints(boundaryEdgesAtPoint,
                             boundaryPoints,
                             cutPoints,
                             &minDistToBoundary,
                             &mdtbargs,
                             extremalPoints);

  selectMinima(extremalPoints,
               boundaryEdgesAtPoint,
               &minDistToBoundary,
               &mdtbargs,
               concaveThresh,
               concavePoints,
               concavity);
}

Real curvature(Edge *e1, Edge *e2)
{
  Vector3 v1,v2;
  if(e1->p1 == e2->p1) {
    v1 = *(e1->p2) - *(e1->p1);
    v2 = *(e2->p2) - *(e2->p1);
  } else if(e1->p1 == e2->p2) {
    v1 = *(e1->p2) - *(e1->p1);
    v2 = *(e2->p1) - *(e2->p2);
  } else if(e1->p2 == e2->p1) {
    v1 = *(e1->p1) - *(e1->p2);
    v2 = *(e2->p2) - *(e2->p1);
  } else if(e1->p2 == e2->p2) {
    v1 = *(e1->p1) - *(e1->p2);
    v2 = *(e2->p1) - *(e2->p2);
  } else {
    abort();
  }
  Vector3 z = v1 * v2;
  z.normalize();
  ZRotation R(z);
  Real x10,y10;
  Real x20,y20;
  R.rotate(v1,&x10,&y10);
  R.rotate(v2,&x20,&y20);
  Real xh = 0.5 * (x10 + x20);
  Real yh = 0.5 * (y10 + y20);
  Real normv12 = norm(v1);
  Real normh2 = Square(xh) + Square(yh);
  if(normh2 / normv12 < 1e-6)
    return 0.;
  Real normhi = 1. / sqrt(normh2);
  Real cosphi = normhi * xh;
  Real sinphi = normhi * yh;
  Real x12 = fabs(cosphi * x10 + sinphi * y10);
  Real y12 = fabs(-sinphi * x10 + cosphi * y10);
  Real x22 = fabs(cosphi * x20 + sinphi * y20);
  Real y22 = fabs(-sinphi * x20 + cosphi * y20);
  Real n = 2. / (x12 + x22);
  Real x0 = x12 * n;
  Real x1 = x22 * n;
  Real y0 = y12 * n;
  Real y1 = 0.;
  Real y2 = y22 * n;
  return 2. * (x0 * (y2 - y1) + x1 * (y0 - y1)) / (x0 * x1 * (x0 + x1));
}

// the cost for a set of boundary edges diverging from a point
Real angleCost(vector<Edge*> &edge)
{
  int edges = edge.size();
  if(edges == 1)
    return 0.;
  if(edges == 2)
    return 2. * curvature(edge[0],edge[1]);

  Real cost = 0.;
  EdgeEdgeMap used;
  for(vector<Edge*>::const_iterator i = edge.begin(); i != edge.end(); ++i) {
    Edge *ei = *i;
    Real curvmin = -1.;
    for(vector<Edge*>::const_iterator j = edge.begin(); j != edge.end(); ++j) {
      Edge *ej = *j;

      EdgeEdgeMap::const_iterator pair1 = used.find(ej);
      if(pair1 != used.end() && pair1->second == ei) 
        continue;
      EdgeEdgeMap::const_iterator pair2 = used.find(ei);
      if(pair2 != used.end() && pair2->second == ej) 
        continue;

      Real curv = curvature(ei,ej);
      if(curvmin == -1. || curv < curvmin) {
        curvmin = curv;
      }

      used[ei] = ej;
      used[ej] = ei;
      cost += curvmin;
    }
  }
  return cost;
}

// the change in angle cost due to adding edges from a cut
Real angleCost(Cut &c)
{
  // the cost without the cut
  Real costOld = angleCost(c.e1) + angleCost(c.e2);

  // the cost with the cut
  vector<Edge*> newE1;
  newE1.insert(newE1.end(),c.e1.begin(),c.e1.end());
  newE1.push_back(c.edge[0]);
  Real costNew = angleCost(newE1);

  // 90 degree angles introduced are free
  Real free = 0;
  if(c.e1.size() == 1) 
    free += 2.;
  else 
    free += 4.;
  if(c.e2.size() == 1) 
    free += 2.;
  else 
    free += 4.;

  // the cost of the cut edges
  vector<Edge*> edge = c.edge;
  Edge *e = NULL;
  for(vector<Edge*>::const_iterator i = edge.begin(); i != edge.end(); ++i) {
    Edge *e2 = *i;
    if(e) costNew += 2. * curvature(e,e2);
    e = e2;
  }

  vector<Edge*> newE2;
  newE2.insert(newE2.end(),c.e2.begin(),c.e2.end());
  newE2.push_back(e);
  costNew += angleCost(newE2);

  return costNew - costOld - free;
}

// the sum over each concave point of the 2nd derivatives of the minimal distance to the boundary
Real concavityCost(PointRealMap &concavity)
{
  Real cost = 0.;
  for(PointRealMap::const_iterator i = concavity.begin(); i != concavity.end(); ++i) {
    cost += i->second;
  }
  return cost;
}

// the change in concavity cost due to a cut
Real concavityCost(Cut &c,
                   Real concaveThresh,
                   PointEdgeVectorMap &edgesAtPoint,
                   PointEdgeVectorMap &boundaryEdgesAtPoint,
                   PointSet &boundaryPoints, 
                   PointSet &cutPoints, 
                   PointRealMap &concavity)
{
  Real costOld = concavityCost(concavity);

  // the new concavity is the old concavity ...
  PointRealMap newConcavity;
  newConcavity.insert(concavity.begin(),concavity.end());

  // minus the concavity due to the cut endpoints ...
  newConcavity.erase(c.p1);
  newConcavity.erase(c.p2);
  
  // plus the concavity introduced by new points
  PointSet newConcavePoints;
  PointSet cPoints;
  cPoints.insert(c.point.begin(),c.point.end());

  findConcavePoints(edgesAtPoint,
                    boundaryEdgesAtPoint,
                    cPoints,
                    cutPoints,
                    concaveThresh,
                    newConcavePoints,
                    newConcavity);

  Real costNew = concavityCost(newConcavity);
  
  return costNew - costOld;
}

// the cost for a cut
// the weighted sum of the angle cost and the concavity cost
Real cutCost(Cut &c,
             Real relativeCost,
             Real concaveThresh,
             PointEdgeVectorMap &edgesAtPoint,
             PointEdgeVectorMap &boundaryEdgesAtPoint,
             PointSet &boundaryPoints, 
             PointSet &cutPoints, 
             PointRealMap &concavity)
{
  return relativeCost * angleCost(c) + concavityCost(c,
                                                     concaveThresh,
                                                     edgesAtPoint,
                                                     boundaryEdgesAtPoint,
                                                     boundaryPoints,
                                                     cutPoints,
                                                     concavity);
}

int neighborTriangle(int i, vector<int> &tap1, vector<int> &tap2)
{
  int ntap1 = tap1.size();
  int ntap2 = tap2.size();

  for(int j=0;j<ntap1;j++) {
    int tap1j = tap1[j];
    if(tap1j != i) {
      for(int k=0;k<ntap2;k++) {
        int tap2k = tap2[k];
        if(tap2k != i)
          if(tap1j == tap2k)
            return tap1j;
      }
    }
  }
  return -1;
}

// NULL triangleAllowed means all trianges can be traversed
void traverse(EdgeSet &edge,
              EdgeSet &boundary,
              TriangleSet &triangle,
              Triangle *t,
              TriangleSet &triangleTraversed,
              EdgeSet &edgeTraversed,
              EdgeSet &boundaryEdges,
              TriangleSet *triangleAllowed)
{
  stack<Triangle*> search;
  search.push(t);

  //cout << "t = [\n";
  while(!search.empty()) {
    t = search.top();
    search.pop();   

    //cout << *t << "\n";
    triangle.insert(t);
    triangleTraversed.insert(t);

    Edge *e12 = t->e12;
    bool bTraversed12 = (edgeTraversed.find(e12) != edgeTraversed.end());
    bool bBoundary12 = (boundaryEdges.find(e12) != boundaryEdges.end());
    if(bBoundary12)
      boundary.insert(e12);    
    if(!bTraversed12) {
      edge.insert(e12);
      edgeTraversed.insert(e12);
      if(!bBoundary12) {
        Triangle *t2 = e12->t1==t?e12->t2:e12->t1;
        if(triangleAllowed == NULL || triangleAllowed->find(t2) != triangleAllowed->end())
          if(triangleTraversed.find(t2) == triangleTraversed.end())
            search.push(t2);
      }
    }
    
    Edge *e23 = t->e23;
    bool bTraversed23 = (edgeTraversed.find(e23) != edgeTraversed.end());
    bool bBoundary23 = (boundaryEdges.find(e23) != boundaryEdges.end());    
    if(bBoundary23)
      boundary.insert(e23);      
    if(!bTraversed23) {
      edge.insert(e23);
      edgeTraversed.insert(e23);
      if(!bBoundary23) {
        Triangle *t2 = e23->t1==t?e23->t2:e23->t1;
        if(triangleAllowed == NULL || triangleAllowed->find(t2) != triangleAllowed->end())
          if(triangleTraversed.find(t2) == triangleTraversed.end())
            search.push(t2);
      }
    }

    Edge *e31 = t->e31;
    bool bTraversed31 = (edgeTraversed.find(e31) != edgeTraversed.end());
    bool bBoundary31 = (boundaryEdges.find(e31) != boundaryEdges.end());    
    if(bBoundary31)
      boundary.insert(e31);      
    if(!bTraversed31) {
      edge.insert(e31);
      edgeTraversed.insert(e31);
      if(!bBoundary31) {
        Triangle *t2 = e31->t1==t?e31->t2:e31->t1;
        if(triangleAllowed == NULL || triangleAllowed->find(t2) != triangleAllowed->end())
          if(triangleTraversed.find(t2) == triangleTraversed.end())
            search.push(t2);
      }
    }
  }
  //cout << "];\n";
}


// Step 1: 
// Construct arrays of points, edges, and triangles from vectors of coordinates and face indices.
//
// Step 2: 
// Segment the mesh into a set of bounded regions, which, with possible exceptions,
// 1. have no "sharp" edges in their interior
// 2. are convex
// 3. minimize the "acuteness" at boundary points
//
// Step 3:
// "Refine" the mesh, so that the edge lengths don't exceed a given threshold
//
// sharpThresh is the threshold (in radians) below which an edge is considered "sharp"
// concaveThresh is the threshold (in unitless curvature) above wich a point is considered concave
// relativeCost determines the relative cost of concavity and "acuteness"; higher relativeCost gives a higher cost to "acuteness"
// minEdgeLength is used for refining edges;  it should return the maximum edge length as a function of position
ShapeMesh :: ShapeMesh(Shape &s, Real sharpThresh, Real concaveThresh, Real relativeCost, bool (*refineEdgeFunc)(Edge *, void *), void *refineEdgeFuncArgs)
{

  // Step 1: Construction
  int triangles = s.faceIndex.size();
  int points = s.coord.size();
  // there may be less, because of possible holes in the shape
  int edges = triangles + points - 2;

  printf("%d %d\n",triangles,points);
  Array<Point> *ppoint = new Array<Point>(points); Array<Point> &point = *ppoint;
  Array<Edge> *pedge = new Array<Edge>(edges); Array<Edge> &edge = *pedge;
  Array<Triangle> *ptriangle = new Array<Triangle>(triangles); Array<Triangle> &triangle = *ptriangle;

  // points first
  for(int i=0;i<points;i++) {
    point[i] = s.coord[i];
  }

  // triangles next - just set the triangle points for now
  vector<int> *trianglesAtPoint = new vector<int>[points];

  for(int i=0;i<triangles;i++) {
    FaceIndex face = s.faceIndex[i];
    Triangle *t = triangle + i;
    t->p1 = point + face.i[0];
    t->p2 = point + face.i[1];
    t->p3 = point + face.i[2];
    trianglesAtPoint[face.i[0]].push_back(i);
    trianglesAtPoint[face.i[1]].push_back(i);
    trianglesAtPoint[face.i[2]].push_back(i);
  }

  // finally, do edges and set triangle edges
  bool *triangleTouched = new bool[triangles];
  for(int i=0;i<triangles;i++) {
    triangleTouched[i] = false;
  }

  edges = 0;
  Key tindex = 0;
  Key eindex = 0;

  for(int i=0;i<triangles;i++) {
    FaceIndex face = s.faceIndex[i];
    Triangle *t = triangle + i;
    t->index = tindex++;
    int neighbor12 = neighborTriangle(i,trianglesAtPoint[face.i[0]],trianglesAtPoint[face.i[1]]);
    int neighbor23 = neighborTriangle(i,trianglesAtPoint[face.i[1]],trianglesAtPoint[face.i[2]]);
    int neighbor31 = neighborTriangle(i,trianglesAtPoint[face.i[2]],trianglesAtPoint[face.i[0]]);

    Edge *e = NULL;
    if(neighbor12 != -1 && !triangleTouched[neighbor12]) {
      e = edge + (edges++);
      e->index = eindex++;
      e->t1 = t;
      e->t2 = triangle + neighbor12;
      e->p1 = point + face.i[0];
      e->p2 = point + face.i[1];
      t->e12 = e;
      if((e->t2->p1 != e->p1) && (e->t2->p1 != e->p2)) {
        e->t2->e23 = e;
      } else if((e->t2->p2 != e->p1) && (e->t2->p2 != e->p2)) {
        e->t2->e31 = e;
      } else {
        e->t2->e12 = e;
      }
    }

    if(neighbor23 != -1 && !triangleTouched[neighbor23]) {
      e = edge + (edges++);
      e->index = eindex++;
      e->t1 = t;
      e->t2 = triangle + neighbor23;
      e->p1 = point + face.i[1];
      e->p2 = point + face.i[2];
      t->e23 = e;
      if((e->t2->p1 != e->p1) && (e->t2->p1 != e->p2)) {
        e->t2->e23 = e;
      } else if((e->t2->p2 != e->p1) && (e->t2->p2 != e->p2)) {
        e->t2->e31 = e;
      } else {
        e->t2->e12 = e;
      }
    }

    if(neighbor31 != -1 && !triangleTouched[neighbor31]) {
      e = edge + (edges++);
      e->index = eindex++;
      e->t1 = t;
      e->t2 = triangle + neighbor31;
      e->p1 = point + face.i[2];
      e->p2 = point + face.i[0];
      t->e31 = e;
      if((e->t2->p1 != e->p1) && (e->t2->p1 != e->p2)) {
        e->t2->e23 = e;
      } else if((e->t2->p2 != e->p1) && (e->t2->p2 != e->p2)) {
        e->t2->e31 = e;
      } else {
        e->t2->e12 = e;
      }
    }
    triangleTouched[i] = true;
  }

  edge.n = edges;
  delete [] triangleTouched;
  delete [] trianglesAtPoint;

  // Step 2: Segmentation
  // get a set of sharp edges
  EdgeSet sharpEdges;

  // these sharp edges also initialize the set of boundary edges
  EdgeSet boundaryEdges;

  // Map points to their edges
  PointEdgeVectorMap boundaryEdgesAtPoint;
  PointEdgeVectorMap edgesAtPoint;

  // also initialize some edge fields
  for(int i=0;i<edges;i++) {
    Edge *e = edge + i;
    Point *p1 = e->p1;
    Point *p2 = e->p2;
    Vector3 v = *p2 - *p1;
    e->length = sqrt(norm(v));
    edgesAtPoint[p1].push_back(e);
    edgesAtPoint[p2].push_back(e);
    Real s = sharpness(e);
    if(s > sharpThresh) {
      sharpEdges.insert(e);
      boundaryEdges.insert(e);
      boundaryEdgesAtPoint[e->p1].push_back(e);
      boundaryEdgesAtPoint[e->p2].push_back(e);
    }
  }

  // get the endpoints of the sharpEdges
  PointSet endPoints;

  // initialize the set of boundary points with the sharp edge points
  PointSet boundaryPoints;

  for(EdgeSet::const_iterator i = sharpEdges.begin(); i != sharpEdges.end(); ++i) {
    Edge *e = *i;
    Point *p1 = e->p1;
    boundaryPoints.insert(p1);
    vector<Edge*> &eap1 = edgesAtPoint[p1];
    bool bConnected1 = false;
    for(vector<Edge*>::const_iterator j = eap1.begin(); j != eap1.end(); ++j) {
      Edge *e1 = *j;
      if(e1 != e) {
        Point *shared = sharedPoint(e,e1);
        if(shared) {
          bConnected1 = true;
          break;
        }
      }
    }
    if(!bConnected1) {
      endPoints.insert(p1);
    }

    Point *p2 = e->p2;
    boundaryPoints.insert(p2);
    vector<Edge*> &eap2 = edgesAtPoint[p2];
    bool bConnected2 = false;
    for(vector<Edge*>::const_iterator j = eap2.begin(); j != eap2.end(); ++j) {
      Edge *e2 = *j;
      if(e2 != e) {
        Point *shared = sharedPoint(e,e2);
        if(shared) {
          bConnected2 = true;
          break;
        }
      }
    }
    if(!bConnected2) {
      endPoints.insert(p2);
    }
  }

  // the cut points are the union of endpoints and concave points 
  // the concave points are yet to be determined
  PointSet cutPoints;
  cutPoints.insert(endPoints.begin(),endPoints.end());

  // get the set of "concave" points that are minima in the minimal distance to other boundary points
  PointSet concavePoints;

  // get the "concavity", i.e. the 2nd derivative of the minimal distance to the boundary
  PointRealMap concavity;

  findConcavePoints(edgesAtPoint,
                    boundaryEdgesAtPoint,
                    boundaryPoints,
                    cutPoints,
                    concaveThresh,
                    concavePoints,
                    concavity);
  cutPoints.insert(concavePoints.begin(),concavePoints.end());

  // enumerate candidate cuts between cut points and other cuts points, 
  // or between cut points and closest points on existing boundaries
  CutSet candidateCuts;

  for(PointSet::const_iterator i = cutPoints.begin(); i != cutPoints.end(); ++i) {
    Point *pi = *i;

    // find the closest boundary points
    PointSet closestPoints;
    findClosestBoundaryPoints(boundaryEdgesAtPoint,
                              boundaryPoints,
                              cutPoints,
                              pi,
                              closestPoints);

    // the destination points are the union of the cut points and the closest points
    PointSet destPoints;
    destPoints.insert(cutPoints.begin(),cutPoints.end());
    destPoints.insert(closestPoints.begin(),closestPoints.end());

    // attempt cuts from cut points to destination points
    for(PointSet::const_iterator j = destPoints.begin(); j != destPoints.end(); ++j) {
      Point *pj = *j;
      // don't try to cut to self
      if(pi == pj) continue;
      Cut *c = new Cut;
      Real distance = circumCut(pi,pj,edgesAtPoint,boundaryEdgesAtPoint,boundaryPoints,c);
      if(distance != DISTANCE_INVALID) candidateCuts.insert(c);
      else delete c;
    }
  }

  // make cuts until there are no more candidates or there are no more endpoints to cut and there
  // is no decrease in cost
  Real minCost = 0.;
  while(!candidateCuts.empty() && (!endPoints.empty() || minCost < 0.)) {
    minCost = -1.;
    Cut *bestCut = NULL;

    // find the candidate that minimizes cost
    for(CutSet::const_iterator j = candidateCuts.begin(); j != candidateCuts.end(); ++j) {
      Cut *cj = *j;
      Real cost = cutCost(*cj,
                          relativeCost,
                          concaveThresh,
                          edgesAtPoint,
                          boundaryEdgesAtPoint,
                          boundaryPoints,
                          cutPoints,
                          concavity);

      if(minCost == -1. || cost < minCost) {
        minCost = cost;
        bestCut = cj;
      }

      if(bestCut) {
        
        // add the cut to the boundary
        for(vector<Edge*>::const_iterator k = bestCut->edge.begin(); k != bestCut->edge.end(); ++k) {
          Edge *e = *k;
          boundaryEdges.insert(e);
          boundaryEdgesAtPoint[e->p1].push_back(e);
          boundaryEdgesAtPoint[e->p2].push_back(e);
        }
        
        for(vector<Point*>::const_iterator k = bestCut->point.begin(); k != bestCut->point.end(); ++k) {
          boundaryPoints.insert(*k);
        }

        // cut endpoints are no longer endpoints
        PointSet::const_iterator p1iter = endPoints.find(bestCut->p1);
        if(p1iter != endPoints.end()) endPoints.erase(p1iter);
        PointSet::const_iterator p2iter = endPoints.find(bestCut->p2);
        if(p2iter != endPoints.end()) endPoints.erase(p2iter);

        // cut concave points are no longer concave
        concavity.erase(bestCut->p1);
        concavePoints.erase(bestCut->p1);
        concavity.erase(bestCut->p2);
        concavePoints.erase(bestCut->p2);

        // remove from the set of candidates
        candidateCuts.erase(candidateCuts.find(bestCut));
        
        // find the other candidates that intersect the new boundary
        PointSet newBoundaryPoints;
        newBoundaryPoints.insert(bestCut->point.begin(),bestCut->point.end());
        CutSet intersectingCuts;

        for(CutSet::const_iterator k = candidateCuts.begin(); k != candidateCuts.end(); ++k) {        
          Cut *ck = *k;
          vector<Point*> &cPoints = ck->point;

          for(vector<Point*>::const_iterator m = cPoints.begin(); m != cPoints.end(); ++m) {
            Point *pm = *m;
            if(newBoundaryPoints.find(pm) != newBoundaryPoints.end()) {
              intersectingCuts.insert(ck);
              break;
            }
          }
        }

        // recompute new cuts 
        // between the endpoints of the intersecting cuts and the new boundary
        for(CutSet::const_iterator k = intersectingCuts.begin(); k != intersectingCuts.end(); ++k) {
          Cut *ck = *k;

          Point *p1 = ck->p1;
          PointSet closestPoints1;

          findClosestBoundaryPoints(boundaryEdgesAtPoint,
                                    newBoundaryPoints,
                                    cutPoints,
                                    p1,
                                    closestPoints1);
          
          for(PointSet::const_iterator m = closestPoints1.begin(); m != closestPoints1.end(); ++m) {
            Point *pm = *m;
            // don't try to cut to self
            if(pm == p1) continue;
            Cut *newc = new Cut;
            Real distance = circumCut(p1,pm,edgesAtPoint,boundaryEdgesAtPoint,boundaryPoints,newc);
            if(distance != DISTANCE_INVALID) candidateCuts.insert(newc);
            else delete newc;
          }

          Point *p2 = ck->p2;
          PointSet closestPoints2;
          findClosestBoundaryPoints(boundaryEdgesAtPoint,
                                    newBoundaryPoints,
                                    cutPoints,
                                    p2,
                                    closestPoints2);
          
          for(PointSet::const_iterator m = closestPoints2.begin(); m != closestPoints2.end(); ++m) {
            Point *pm = *m;
            // don't try to cut to self
            if(pm == p2) continue;
            Cut *newc = new Cut;
            Real distance = circumCut(p2,pm,edgesAtPoint,boundaryEdgesAtPoint,boundaryPoints,newc);
            if(distance != DISTANCE_INVALID) candidateCuts.insert(newc);
            else delete newc;
          }
        }

        // find new concave points
        PointSet newConcavePoints;
        findConcavePoints(edgesAtPoint,
                          boundaryEdgesAtPoint,
                          newBoundaryPoints,
                          cutPoints,
                          concaveThresh,
                          newConcavePoints,
                          concavity);
        cutPoints.insert(newConcavePoints.begin(),newConcavePoints.end());
        
        // find new candidates originating from the new concave points
        for(PointSet::const_iterator k = newConcavePoints.begin(); k != newConcavePoints.end(); ++k) {
          Point *pk = *k;

          PointSet closestPointsc;

          findClosestBoundaryPoints(boundaryEdgesAtPoint,
                                    newBoundaryPoints,
                                    cutPoints,
                                    pk,
                                    closestPointsc);

          PointSet destPoints;
          destPoints.insert(cutPoints.begin(),cutPoints.end());
          destPoints.insert(closestPointsc.begin(),closestPointsc.end());

          // attempt cuts from new concave points to destination points
          for(PointSet::const_iterator m = destPoints.begin(); m != destPoints.end(); ++m) {
            Point *pm = *m;
            // don't try to cut to self
            if(pm == pk) continue;
            Cut *newc = new Cut;
            Real distance = circumCut(pk,pm,edgesAtPoint,boundaryEdgesAtPoint,boundaryPoints,newc);
            if(distance != DISTANCE_INVALID) candidateCuts.insert(newc);
            else delete newc;
          }
        }

        // delete the cut
        delete bestCut;
      }
    }
  }

  // delete the cuts
  for(CutSet::const_iterator j = candidateCuts.begin(); j != candidateCuts.end(); ++j) {
    delete *j;
  }

  // Step 3: Refinement
  
  bool bRefine;
  do {
    // add the previous refinement
    refinedPoint.push_back(ppoint);
    refinedEdge.push_back(pedge);
    refinedTriangle.push_back(ptriangle);

    bRefine = false;

    EdgeSet edgesToRefine;
    TriangleSet trianglesToRefine;

    // count new points
    int newPoints = 0;

    for(vector< Array<Edge>* >::const_iterator j = refinedEdge.begin(); j != refinedEdge.end(); ++j) {
      Array<Edge> &a = *(*j);
      int n = a.n;
      for(int i=0;i<n;i++) {
        Edge *e = a + i;
        Point p = 0.5 * (*(e->p1) + *(e->p2));
        if(refineEdgeFunc != NULL && refineEdgeFunc(e,refineEdgeFuncArgs)) {
          edgesToRefine.insert(e);
          trianglesToRefine.insert(e->t1);
          trianglesToRefine.insert(e->t2);
          newPoints++;
          bRefine = true;
        }
      }
    }
            
    if(bRefine) {

      // allocate new arrays
      int newEdges = newPoints;
      ppoint = new Array<Point>(newPoints); Array<Point> &point = *ppoint;
      pedge = new Array<Edge>(newEdges); Array<Edge> &edge = *pedge;
      points = 0;
      edges = 0;
      
      // count new triangles
      int newTriangles = 0;

      // Map edges to new edges
      EdgeEdgeMap edgeRefined;

      // add new points and exterior edges
      for(TriangleSet::const_iterator i = trianglesToRefine.begin(); i != trianglesToRefine.end(); ++i) {
        Triangle *t = *i;

        Point *p1 = t->p1;
        Point *p2 = t->p2;
        Point *p3 = t->p3;

        Edge *e12 = t->e12;

        if(edgesToRefine.find(e12) != edgesToRefine.end()) {

          newTriangles++;

          if(edgeRefined.find(e12) == edgeRefined.end()) {

            Point *p12;
            Edge *e122;

            // new point
            p12 = point + points++;
            *p12 = 0.5 * (*(e12->p1) + *(e12->p2));

            Real length = 0.5 * e12->length;
            
            // new edge
            e122 = edge + edges++;
            e122->index = eindex++;
            e122->p1 = p12;
            e122->p2 = p2;
            e122->length = length;
            e122->t1 = NULL;
            e122->t2 = NULL;
            
            // fix the old edge
            e12->p1 = p1;
            e12->p2 = p12;
            e12->length = length;
            e12->t1 = NULL;
            e12->t2 = NULL;
            
            if(boundaryEdges.find(e12) != boundaryEdges.end()) {
              boundaryEdges.insert(e122);
            }
            
            edgeRefined[e12] = e122;
          }
        }
        
        Edge *e23 = t->e23;
        if(edgesToRefine.find(e23) != edgesToRefine.end()) {

          newTriangles++;
          
          if(edgeRefined.find(e23) == edgeRefined.end()) {

            Point *p23;
            Edge *e232;
            
            // new point
            p23 = point + points++;
            *p23 = 0.5 * (*(e23->p1) + *(e23->p2));

            Real length = 0.5 * e23->length;

            // new edge
            e232 = edge + edges++;
            e232->index = eindex++;
            e232->p1 = p23;
            e232->p2 = p3;
            e232->length = length;
            e232->t1 = NULL;
            e232->t2 = NULL;
            
            // fix the old edge
            e23->p1 = p2;
            e23->p2 = p23;
            e23->length = length;
            e23->t1 = NULL;
            e23->t2 = NULL;
            
            if(boundaryEdges.find(e23) != boundaryEdges.end()) {
              boundaryEdges.insert(e232);
            }
            
            edgeRefined[e23] = e232;
          }
        }

        Edge *e31 = t->e31;

        if(edgesToRefine.find(e31) != edgesToRefine.end()) {

          newTriangles++;

          if(edgeRefined.find(e31) == edgeRefined.end()) {

            Point *p31;
            Edge *e312;

            // new point
            p31 = point + points++;
            *p31 = 0.5 * (*(e31->p1) + *(e31->p2));

            Real length = 0.5 * e31->length;
            
            // new edge
            e312 = edge + edges++;
            e312->index = eindex++;
            e312->p1 = p31;
            e312->p2 = p1;
            e312->length = length;
            e312->t1 = NULL;
            e312->t2 = NULL;
            
            // fix the old edge
            e31->p1 = p3;
            e31->p2 = p31;
            e31->length = length;
            e31->t1 = NULL;
            e31->t2 = NULL;

            if(boundaryEdges.find(e31) != boundaryEdges.end()) {
              boundaryEdges.insert(e312);
            }
            
            edgeRefined[e31] = e312;
          }
        }
      }
   
      // add the new external edges
      refinedEdge.push_back(pedge);

      // more new arrays
      newEdges = newTriangles;
      pedge = new Array<Edge>(newEdges); Array<Edge> &edge2 = *pedge;
      ptriangle = new Array<Triangle>(newTriangles); Array<Triangle> &triangle = *ptriangle;
      edges = 0;
      triangles = 0;

      // index the new triangles
      for(int i=0; i<newTriangles; i++) {
        Triangle *t = triangle + i;
        t->index = tindex++;
      }

      // add new interior edges and sub triangles
      for(TriangleSet::const_iterator i = trianglesToRefine.begin(); i != trianglesToRefine.end(); ++i) {
        Triangle *t = *i;
        
        int newSubTriangles = 0;

        Point *p1 = t->p1;
        Point *p2 = t->p2;
        Point *p3 = t->p3;

        Edge *e12 = t->e12;
        Point *p12 = NULL;
        Edge *e122 = NULL;
        EdgeEdgeMap::const_iterator e12i = edgeRefined.find(e12);
        if(e12i != edgeRefined.end()) {
          e122 = e12i->second;
          p12 = e122->p1;
          if(e12->p1 != p1 && e12->p2 != p1) swap(e12,e122);
          newSubTriangles++;
        }
        
        Edge *e23 = t->e23;
        Point *p23 = NULL;
        Edge *e232 = NULL;
        EdgeEdgeMap::const_iterator e23i = edgeRefined.find(e23);
        if(e23i != edgeRefined.end()) {
          e232 = e23i->second;
          p23 = e232->p1;
          if(e23->p1 != p2 && e23->p2 != p2) swap(e23,e232);
          newSubTriangles++;
        }

        Edge *e31 = t->e31;
        Point *p31 = NULL;
        Edge *e312 = NULL;
        EdgeEdgeMap::const_iterator e31i = edgeRefined.find(e31);
        if(e31i != edgeRefined.end()) {
          e312 = e31i->second;
          p31 = e312->p1;
          if(e31->p1 != p3 && e31->p2 != p3) swap(e31,e312);
          newSubTriangles++;
        }
        
        if(newSubTriangles == 0) {
          abort();
        } else if(newSubTriangles == 1) {
          if(p12) {

            // new edge
            Edge *e1 = edge2 + edges++;
            e1->index = eindex++;
            e1->p1 = p3;
            e1->p2 = p12;
            Vector3 v = *(e1->p1) - *(e1->p2);
            e1->length = sqrt(norm(v));
            
            // new triangle
            Triangle *t1 = triangle + triangles++;
            t1->index = tindex++;
            t1->p1 = p12;
            t1->p2 = p2;
            t1->p3 = p3;
            t1->e12 = e122;
            t1->e23 = e23;
            t1->e31 = e1;

            // fix the edge triangle
            if(e23->t1 == t) e23->t1 = t1; else e23->t2 = t1;

            // fix the old triangle
            t->p2 = p12;
            t->e12 = e12;
            t->e23 = e1;

            // set the edge triangles
            if(e12->t1) e12->t2 = t; else e12->t1 = t;
            if(e122->t1) e122->t2 = t1; else e122->t1 = t1;
            e1->t1 = t;
            e1->t2 = t1;
          } else if(p23) {

            // new edge
            Edge *e1 = edge2 + edges++;
            e1->index = eindex++;
            e1->p1 = p1;
            e1->p2 = p23;
            Vector3 v = *(e1->p1) - *(e1->p2);
            e1->length = sqrt(norm(v));
              
            // new triangle
            Triangle *t1 = triangle + triangles++;
            t1->index = tindex++;
            t1->p1 = p23;
            t1->p2 = p3;
            t1->p3 = p1;
            t1->e12 = e232;
            t1->e23 = e31;
            t1->e31 = e1;

            // fix the edge triangle
            if(e31->t1 == t) e31->t1 = t1; else e31->t2 = t1;
              
            // fix the old triangle
            t->p3 = p23;
            t->e23 = e23;
            t->e31 = e1;

            // set the edge triangles
            if(e23->t1) e23->t2 = t; else e23->t1 = t;
            if(e232->t1) e232->t2 = t1; else e232->t1 = t1;
            e1->t1 = t;
            e1->t2 = t1;
          } else /* if(p31) */ {
            
            // new edge
            Edge *e1 = edge2 + edges++;
            e1->index = eindex++;
            e1->p1 = p2;
            e1->p2 = p31;
            Vector3 v = *(e1->p1) - *(e1->p2);
            e1->length = sqrt(norm(v));
            
            // new triangle
            Triangle *t1 = triangle + triangles++;
            t1->index = tindex++;
            t1->p1 = p31;
            t1->p2 = p1;
            t1->p3 = p2;
            t1->e12 = e312;
            t1->e23 = e12;
            t1->e31 = e1;

            // fix the edge triangle
            if(e12->t1 == t) e12->t1 = t1; else e12->t2 = t1;
              
            // fix the old triangle
            t->p1 = p31;
            t->e31 = e31;
            t->e12 = e1;

            // set the edge triangles
            if(e31->t1) e31->t2 = t; else e31->t1 = t;
            if(e312->t1) e312->t2 = t1; else e312->t1 = t1;
            e1->t1 = t;
            e1->t2 = t1;
          }
        } else if(newSubTriangles == 2) {
          if(!p12) {

            // new edge 1
            Edge *e1 = edge2 + edges++;
            e1->index = eindex++;
            e1->p1 = p23;
            e1->p2 = p31;
            e1->length = 0.5 * e12->length;

            // new triangle 1
            Triangle *t1 = triangle + triangles++;
            t1->index = tindex++;
            t1->p1 = p31;
            t1->p2 = p23;
            t1->p3 = p3;
            t1->e12 = e1;
            t1->e23 = e232;
            t1->e31 = e31;
            
            Real angle1 = angle(e312,e12);
            Real angle2 = angle(e12,e23);
            
            if(angle1 >= angle2) {
 
              // new edge 2
              Edge *e2 = edge2 + edges++;
              e2->index = eindex++;
              e2->p1 = p1;
              e2->p2 = p23;
              Vector3 v = *(e2->p1) - *(e2->p2);
              e2->length = sqrt(norm(v));
              
              // new triangle 2
              Triangle *t2 = triangle + triangles++;
              t2->index = tindex++;
              t2->p1 = p31;
              t2->p2 = p1;
              t2->p3 = p23;
              t2->e12 = e312;
              t2->e23 = e2;
              t2->e31 = e1;
              
              // fix old triangle
              t->p3 = p23;
              t->e23 = e23;
              t->e31 = e2;

              // set the edge triangles
              if(e23->t1) e23->t2 = t; else e23->t1 = t;
              if(e232->t1) e232->t2 = t1; else e232->t1 = t1;
              if(e31->t1) e31->t2 = t1; else e31->t1 = t1;
              if(e312->t1) e312->t2 = t2; else e312->t1 = t2;
              e1->t1 = t1;
              e1->t2 = t2;
              e2->t1 = t;
              e2->t2 = t2;
            } else {
              
              // new edge 2
              Edge *e2 = edge2 + edges++;
              e2->index = eindex++;
              e2->p1 = p2;
              e2->p2 = p31;
              Vector3 v = *(e2->p1) - *(e2->p2);
              e2->length = sqrt(norm(v));
              
              // new triangle 2
              Triangle *t2 = triangle + triangles++;
              t2->index = tindex++;
              t2->p1 = p23;
              t2->p2 = p31;
              t2->p3 = p2;
              t2->e12 = e1;
              t2->e23 = e2;
              t2->e31 = e23;

              // fix old triangle
              t->p3 = p31;
              t->e23 = e2;
              t->e31 = e312;
              
              // set the edge triangles
              if(e23->t1) e23->t2 = t2; else e23->t1 = t2;
              if(e232->t1) e232->t2 = t1; else e232->t1 = t1;
              if(e31->t1) e31->t2 = t1; else e31->t1 = t1;
              if(e312->t1) e312->t2 = t; else e312->t1 = t;
              e1->t1 = t1;
              e1->t2 = t2;
              e2->t1 = t;
              e2->t2 = t2;
            }
          } else if(!p23) {
            
            // new edge 1
            Edge *e1 = edge2 + edges++;
            e1->index = eindex++;
            e1->p1 = p12;
            e1->p2 = p31;
            e1->length = 0.5 * e23->length;
            
            // new triangle 1
            Triangle *t1 = triangle + triangles++;
            t1->index = tindex++;
            t1->p1 = p12;
            t1->p2 = p31;
            t1->p3 = p1;
            t1->e12 = e1;
            t1->e23 = e312;
            t1->e31 = e12;

            Real angle2 = angle(e122,e23);
            Real angle3 = angle(e23,e31);
            
            if(angle2 >= angle3) {
 
              // new edge 2
              Edge *e2 = edge2 + edges++;
              e2->index = eindex++;
              e2->p1 = p2;
              e2->p2 = p31;
              Vector3 v = *(e2->p1) - *(e2->p2);
              e2->length = sqrt(norm(v));
                
              // new triangle 2
              Triangle *t2 = triangle + triangles++;
              t2->index = tindex++;
              t2->p1 = p12;
              t2->p2 = p2;
              t2->p3 = p31;
              t2->e12 = e122;
              t2->e23 = e2;
              t2->e31 = e1;
              
              // fix old triangle
              t->p1 = p31;
              t->e31 = e31;
              t->e12 = e2;

              // set the edge triangles
              if(e31->t1) e31->t2 = t; else e31->t1 = t;
              if(e312->t1) e312->t2 = t1; else e312->t1 = t1;
              if(e12->t1) e12->t2 = t1; else e12->t1 = t1;
              if(e122->t1) e122->t2 = t2; else e122->t1 = t2;
              e1->t1 = t1;
              e1->t2 = t2;
              e2->t1 = t;
              e2->t2 = t2;
            } else {
 
              // new edge 2
              Edge *e2 = edge2 + edges++;
              e2->index = eindex++;
              e2->p1 = p3;
              e2->p2 = p12;
              Vector3 v = *(e2->p1) - *(e2->p2);
              e2->length = sqrt(norm(v));
              
              // new triangle 2
              Triangle *t2 = triangle + triangles++;
              t2->index = tindex++;
              t2->p1 = p31;
              t2->p2 = p12;
              t2->p3 = p3;
              t2->e12 = e1;
              t2->e23 = e2;
              t2->e31 = e31;

              // fix old triangle
              t->p1 = p12;
              t->e12 = e122;
              t->e31 = e2;
              
              // set the edge triangles
              if(e31->t1) e31->t2 = t2; else e31->t1 = t2;
              if(e312->t1) e312->t2 = t1; else e312->t1 = t1;
              if(e12->t1) e12->t2 = t1; else e12->t1 = t1;
              if(e122->t1) e122->t2 = t; else e122->t1 = t;
              e1->t1 = t1;
              e1->t2 = t2;
              e2->t1 = t;
              e2->t2 = t2;
            }
          } else /* if(!p31) */ {
            
            // new edge 1
            Edge *e1 = edge2 + edges++;
            e1->index = eindex++;
            e1->p1 = p23;
            e1->p2 = p12;
            e1->length = 0.5 * e31->length;
            
            // new triangle 1
            Triangle *t1 = triangle + triangles++;
            t1->index = tindex++;
            t1->p1 = p23;
            t1->p2 = p12;
            t1->p3 = p2;
            t1->e12 = e1;
            t1->e23 = e122;
            t1->e31 = e23;
            
            Real angle3 = angle(e232,e31);
            Real angle1 = angle(e31,e12);
            
            if(angle3 >= angle1) {
              
              // new edge 2
              Edge *e2 = edge2 + edges++;
              e2->index = eindex++;
              e2->p1 = p3;
              e2->p2 = p12;
              Vector3 v = *(e2->p1) - *(e2->p2);
              e2->length = sqrt(norm(v));
              
              // new triangle 2
              Triangle *t2 = triangle + triangles++;
              t2->index = tindex++;
              t2->p1 = p23;
              t2->p2 = p3;
              t2->p3 = p12;
              t2->e12 = e232;
              t2->e23 = e2;
              t2->e31 = e1;
              
              // fix old triangle
              t->p2 = p12;
              t->e12 = e12;
              t->e23 = e2;
              
              // set the edge triangles
              if(e12->t1) e12->t2 = t; else e12->t1 = t;
              if(e122->t1) e122->t2 = t1; else e122->t1 = t1;
              if(e23->t1) e23->t2 = t1; else e23->t1 = t1;
              if(e232->t1) e232->t2 = t2; else e232->t1 = t2;
              e1->t1 = t1;
              e1->t2 = t2;
              e2->t1 = t;
              e2->t2 = t2;
            } else {
              
              // new edge 2
              Edge *e2 = edge2 + edges++;
              e2->index = eindex++;
              e2->p1 = p1;
              e2->p2 = p23;
              Vector3 v = *(e2->p1) - *(e2->p2);
              e2->length = sqrt(norm(v));
              
              // new triangle 2
              Triangle *t2 = triangle + triangles++;
              t2->index = tindex++;
              t2->p1 = p12;
              t2->p2 = p23;
              t2->p3 = p1;
              t2->e12 = e1;
              t2->e23 = e2;
              t2->e31 = e12;
              
              // fix old triangle
              t->p2 = p23;
              t->e12 = e2;
              t->e23 = e232;
              
              // set the edge triangles 
              if(e12->t1) e12->t2 = t2; else e12->t1 = t2;
              if(e122->t1) e122->t2 = t1; else e122->t1 = t1;
              if(e23->t1) e23->t2 = t1; else e23->t1 = t1;
              if(e232->t1) e232->t2 = t; else e232->t1 = t;
              e1->t1 = t1;
              e1->t2 = t2;
              e2->t1 = t;
              e2->t2 = t2;
            }
          }
        } else /* if(newSubTriangles == 3) */ {
          
          // new edge 1
          Edge *e1 = edge2 + edges++;
          e1->index = eindex++;
          e1->p1 = p23;
          e1->p2 = p12;
          e1->length = e31->length;
          
          // new edge 2
          Edge *e2 = edge2 + edges++;
          e2->index = eindex++;
          e2->p1 = p31;
          e2->p2 = p23;
          e2->length = e12->length;
          
          // new edge 3
          Edge *e3 = edge2 + edges++;
          e3->index = eindex++;
          e3->p1 = p12;
          e3->p2 = p31;
          e3->length = e23->length;

          // new triangle 1
          Triangle *t1 = triangle + triangles++;
          t1->index = tindex++;
          t1->p1 = p12;
          t1->p2 = p2;
          t1->p3 = p23;
          t1->e12 = e122;
          t1->e23 = e23;
          t1->e31 = e1;
          
          // new triangle 2
          Triangle *t2 = triangle + triangles++;
          t2->index = tindex++;
          t2->p1 = p31;
          t2->p2 = p23;
          t2->p3 = p3;
          t2->e12 = e2;
          t2->e23 = e232;
          t2->e31 = e31;

          // new triangle 3
          Triangle *t3 = triangle + triangles++;
          t3->index = tindex++;
          t3->p1 = p23;
          t3->p2 = p31;
          t3->p3 = p12;
          t3->e12 = e2;
          t3->e23 = e3;
          t3->e31 = e1;

          // fix the old triangle
          t->p2 = p12;
          t->p3 = p31;
          t->e12 = e12;
          t->e23 = e3;
          t->e31 = e312;
          
          // set edge triangles
          if(e12->t1) e12->t2 = t; else e12->t1 = t;
          if(e122->t1) e122->t2 = t1; else e122->t1 = t1;
          if(e23->t1) e23->t2 = t1; else e23->t1 = t1;
          if(e232->t1) e232->t2 = t2; else e232->t1 = t2;
          if(e31->t1) e31->t2 = t2; else e31->t1 = t2;
          if(e312->t1) e312->t2 = t; else e312->t1 = t;
          e1->t1 = t1;
          e1->t2 = t3;
          e2->t1 = t2;
          e2->t2 = t3;
          e3->t1 = t;
          e3->t2 = t3;
        }
      }
    }
  } while(bRefine);

  // determine the segments
  EdgeSet edgeTraversed;
  TriangleSet triangleTraversed;
  for(vector< Array<Triangle>* >::const_iterator i = refinedTriangle.begin(); i != refinedTriangle.end(); ++i) {
    Array<Triangle> &a = *(*i);
    int n = a.n;
    for(int i=0; i<n; i++ ) {
      Triangle *t = a + i;
      if(triangleTraversed.find(t) != triangleTraversed.end()) continue;
      Segment *s = new Segment;
      segment.insert(s);
      traverse(s->edge,s->boundary,s->triangle,
               t,triangleTraversed,edgeTraversed,boundaryEdges);
    }
  }
}

template<typename T>
void deletePointers(vector<T*> &v)
{
  for(typename vector<T*>::const_iterator i = v.begin(); i != v.end(); ++i) {
    delete *i;
  }
}

unsigned long ShapeMesh :: getTriangleCount()
{
  unsigned long t = 0;
  for(vector< Array<Triangle>* >::const_iterator i = refinedTriangle.begin(); i != refinedTriangle.end(); ++i) {
    Array<Triangle> *triangle = *i;
    t += triangle->n;
  }
  return t;
}

unsigned long ShapeMesh :: getPointCount()
{
  unsigned long p = 0;
  for(vector< Array<Point>* >::const_iterator i = refinedPoint.begin(); i != refinedPoint.end(); ++i) {
    Array<Point> *point = *i;
    p += point->n;
  }
  return p;
}

ShapeMesh :: ~ShapeMesh()
{
  deletePointers< Array<Point> >(refinedPoint);
  deletePointers< Array<Edge> >(refinedEdge);
  deletePointers< Array<Triangle> >(refinedTriangle);
}


void getEdgesAtPoints(ShapeMesh &m, PointEdgeVectorMap& edgesAtPoint)
{
  SegmentSet &segment = m.segment;
  for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
    Segment *s = *i;
    for(EdgeSet::const_iterator k = s->edge.begin(); k != s->edge.end(); ++k) {
      Edge *e = *k;
      edgesAtPoint[e->p1].push_back(e);
      edgesAtPoint[e->p2].push_back(e);
    }
  }
}

ostream& operator<<(ostream& os, const ShapeMesh &m)
{
  int k = 1;
  for(SegmentSet::const_iterator i = m.segment.begin(); i != m.segment.end(); ++i) {
    Segment *s = *i;
    os << "s{" << k++ << "} = [\n";
    for(EdgeSet::const_iterator j = s->edge.begin(); j != s->edge.end(); ++j) {
      Edge *e = *j;
      //os << "[" << center(e->t1) << " " << normal(e->t1) << "];\n";
      //os << "[" << center(e->t2) << " " << normal(e->t2) << "];\n";
      os << *(e->t1) << ";\n" << *(e->t2) << ";\n";
    }
    os << "];\n";
  }

  return os;
}

