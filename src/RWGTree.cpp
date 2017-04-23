#include "RWGTree.h"
#include "MathUtils.h"
#include "Rotation.h"

#define DEBUG 0

class CrossSectionBoundaryLengthArgs {
public:
  CrossSectionBoundaryLengthArgs(Point &centroid_, 
                                 EdgeSet &boundary_, 
                                 PointEdgeVectorMap &edgesAtPoint_)
    : centroid(centroid_), boundary(boundary_), edgesAtPoint(edgesAtPoint_) {}
  Point &centroid;
  RWG **rwg;
  int rwgs;
  EdgeSet &boundary;
  PointEdgeVectorMap &edgesAtPoint;
  PointAndVector *pv;
};

void segmentInfo(Segment *s, Real *radius, Real *sarea, Point *centroid, Vector3 *n)
{
  *centroid = Point(0.,0.,0.);
  *n = Vector3(0.,0.,0.);

  Real norma = 0.;
  for(TriangleSet::const_iterator i = s->triangle.begin(); i != s->triangle.end(); ++i) {
    Triangle *t = *i;
    Real a = area(t);
    Point p(0.,0.,0.);
    p += *(t->p1);
    p += *(t->p2);
    p += *(t->p3);
    p *= ONETHIRD * a;
    *centroid += p;
    *n += a * normal(t);
    norma += a;
  }
  
  *centroid *= 1. / norma;
  *n *= 1. / norma;
  if(norm(*n) < 1e-7)
    *n = Vector3(0.,0.,1.);
  else
    n->normalize();

  Real maxR2 = 0.;
  for(EdgeSet::const_iterator i = s->edge.begin(); i != s->edge.end(); ++i) {
    Edge *e = *i;
    Point *p1 = e->p1;
    Point *p2 = e->p2;
    Vector3 v;
    Real r2;
    v = *(p1) - *(centroid);
    r2 = norm(v);
    if(r2 > maxR2)
      maxR2 = r2;
    v = *(p2) - *(centroid);
    r2 = norm(v);
    if(r2 > maxR2)
      maxR2 = r2;
  }
  *radius = sqrt(maxR2);
  *sarea = norma;
}

void RWGTree :: partition(Real *angles,
                          int minSize,
                          Real minRadius)
{
  //printf("Attempting partition with angles [%g %g %g]\n",angles[0],angles[1],angles[2]);
  Real sinphi = sin(angles[0]);
  Real cosphi = cos(angles[0]);
  Real sintheta = sin(angles[1]);
  Real costheta = cos(angles[1]);
  Real x = sinphi * costheta;
  Real y = sinphi * sintheta;
  Real z = cosphi;
  Vector3 axis(x,y,z);
  Rotation R(axis,angles[2]);

  // gather all edges into a set
  EdgeSet edge;
  // find the edges that cross the z=0 plane
  EdgeSet cuts;
  // gather some sets of points
  PointSet nonBoundaryCutPoints;
  PointSet boundaryCutPoints;
  PointSet cutPoints;
  // Map points to their distance to the z=0 plane
  PointRealMap pointDistance;
  
  for(int i=0; i<rwgs; i++) {
    RWG *r = rwg[i];
    Edge *e = r->e;
    edge.insert(e);
    Point p1 = *(e->p1) - centroid;
    Point Rp1 = R.rotate(p1);
    Real z1 = Rp1.z;
    Point p2 = *(e->p2) - centroid;
    Point Rp2 = R.rotate(p2);
    Real z2 = Rp2.z;
    pointDistance[e->p1] = z1;
    pointDistance[e->p2] = z2;

    // if the edge touches the z=0 plane it is a cut
    Real small = 1e-6 * e->length;
        if((z1 <= small && z2 >= -small) || (z1 >= -small && z2 <= small)) {

      cuts.insert(e);

      // check to see if the endpoints are on the boundary      
      vector<Edge*> &eap1 = edgesAtPoint[e->p1];
      bool bBoundary1 = false;
      for(vector<Edge*>::const_iterator j = eap1.begin(); j != eap1.end(); ++j) {
        Edge *e1 = *j;
        if(boundary.find(e1) != boundary.end()) {
          bBoundary1 = true;
          break;
        }
      }
      vector<Edge*> &eap2 = edgesAtPoint[e->p2];
      bool bBoundary2 = false;
      for(vector<Edge*>::const_iterator j = eap2.begin(); j != eap2.end(); ++j) {
        Edge *e2 = *j;
        if(boundary.find(e2) != boundary.end()) {
          bBoundary2 = true;
          break;
        }
      }

      // add to the sets
      if(bBoundary1) boundaryCutPoints.insert(e->p1);
      else nonBoundaryCutPoints.insert(e->p1);
      if(bBoundary2) boundaryCutPoints.insert(e->p2);
      else nonBoundaryCutPoints.insert(e->p2);
      cutPoints.insert(e->p1);
      cutPoints.insert(e->p2);
    }
  }

  // the set of new boundary edges
  EdgeSet newBoundary;

  // the set of new points that intersect the boundary
  PointSet newBoundaryCutPoints;

  // traverse the cuts, starting with the boundary points, and moving on to non boundary points
  vector<Point*> startPoints;
  startPoints.insert(startPoints.end(),boundaryCutPoints.begin(),boundaryCutPoints.end());
  startPoints.insert(startPoints.end(),nonBoundaryCutPoints.begin(),nonBoundaryCutPoints.end());

  if(DEBUG)
    cout << "Starting with " << startPoints.size() << " cut points\n";

  bool bHitBoundary = false;

  for(vector<Point*>::const_iterator i = startPoints.begin(); i != startPoints.end(); ++i) {
    Point *p = *i;

    bool boundaryStart;
    if(boundaryCutPoints.find(p) != boundaryCutPoints.end() &&
       newBoundaryCutPoints.find(p) == newBoundaryCutPoints.end())
      boundaryStart = true;
    else
      boundaryStart = false;

    // first check to see if there is an edge that connects to a new boundary cut point
    // in which case the point is skipped
    bool bSkip = false;
    vector<Edge*> &eap = edgesAtPoint[p];
    for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
      Edge *e = *j;
      Point *p1 = (e->p1 == p) ? e->p2 : e->p1;
      if(newBoundaryCutPoints.find(p1) != newBoundaryCutPoints.end()) {
        bSkip = true;
        if(DEBUG) {
          if(boundaryStart)
            cout << "boundary start: skipping\n";
          else
            cout << "non boundary start: skipping\n";
        }
        break;
      }
    }
    if(bSkip) continue;
   
    // find the best edge to cut through
    Vector3 direction;
    Edge *minePrev = NULL;
    int n = 0;
    //cout << "c = [\n";
    do {
      Real minMerit = Infinity;
      Edge *mine = NULL;
      Point *minp = NULL;

      //cout << *p << "\n";
      vector<Edge*> &eap = edgesAtPoint[p];
      bool bBoundary = (boundaryCutPoints.find(p) != boundaryCutPoints.end());

      // now find the best point
      for(vector<Edge*>::const_iterator j = eap.begin(); j != eap.end(); ++j) {
        Edge *e = *j;

        // dont go backwards
        if(e == minePrev) continue;

        Point *p1 = (e->p1 == p) ? e->p2 : e->p1;

        // don't follow an existing boundary
        bool bBoundary1 = (boundaryCutPoints.find(p1) != boundaryCutPoints.end());
        if(bBoundary && bBoundary1) {
          if(DEBUG) cout << "already a boundary\n";
          continue;
        }
        
        
        // the point must be a cut point
        if(cutPoints.find(p1) == cutPoints.end()) {
          if(DEBUG) cout << "not a cut\n";
          continue;
        }
        

        Real merit;
        Real forward;

        if(minePrev) {
          Vector3 v0 = *(p) + e->length * direction; 
          Vector3 v = v0 - *(p1);
          forward = norm(v);
        } else {
          forward = 0.;
        }

        merit = (forward + Square(pointDistance[p1])) / Square(e->length);
        if(DEBUG) {
          cout << *p1 << " merit (forward) = " << merit << "\n";
        }

        if(merit < minMerit) {
          minMerit = merit;
          mine = e;
          minp = p1;
        }
      }
 
      // traverse
      if(minp) {
        
        newBoundary.insert(mine);
        newBoundaryCutPoints.insert(p);
        boundaryCutPoints.insert(p);
        nonBoundaryCutPoints.erase(p);

        // stop if it hits a boundary 
        if(boundaryCutPoints.find(minp) != boundaryCutPoints.end()) {
          newBoundaryCutPoints.insert(minp);
          if(DEBUG) cout << "hit boundary\n";
          bHitBoundary = true;
          p = NULL;
        } else {
          // move on
          direction = *minp - *p;
          direction *= 1. / mine->length;
          minePrev = mine;
          p = minp;
          n++;
          if(DEBUG) cout << "traversed " << n << "\n";
        }
      } else {
        p = NULL;
      }
    } while(p);
    //cout << "];\n";
    //exit(0);
    if(DEBUG) {
      if(boundaryStart)
        cout << "boundary start: traversed " << n << "\n";
      else
        cout << "non boundary start: traversed " << n << "\n";
    }
  }

  if(bHitBoundary) {
    // prune dangling boundaries
    bool bPruned;
    do {
      bPruned = false;
      EdgeSet pruneMe;
      for(EdgeSet::const_iterator i = newBoundary.begin(); i != newBoundary.end(); ++i) {
        Edge *e = *i;
        Point *p1 = e->p1;
        Point *p2 = e->p2;
        
        vector<Edge*> &eap1 = edgesAtPoint[p1];
        bool bConnected1 = false;
        for(vector<Edge*>::const_iterator j = eap1.begin(); j != eap1.end(); ++j) {
          Edge *e1 = *j;
          if(e == e1) continue;
          if(newBoundary.find(e1) != newBoundary.end() || boundary.find(e1) != boundary.end()) {
            bConnected1 = true;
            break;
          }
        }
        if(!bConnected1) {
          pruneMe.insert(e);
          continue;
        }
        
        vector<Edge*> &eap2 = edgesAtPoint[p2];
        bool bConnected2 = false;
        for(vector<Edge*>::const_iterator j = eap2.begin(); j != eap2.end(); ++j) {
          Edge *e2 = *j;
          if(e == e2) continue;
          if(newBoundary.find(e2) != newBoundary.end() || boundary.find(e2) != boundary.end()) {
            bConnected2 = true;
            break;
          }
        }
        if(!bConnected2) {
          pruneMe.insert(e);
          continue;
        }    
      }
      
      if(DEBUG&&pruneMe.size()) {
        printf("Pruning %d\n",(int)pruneMe.size());
      }
      
      // erase the prunes
      for(EdgeSet::const_iterator i = pruneMe.begin(); i != pruneMe.end(); ++i) {
        bPruned = true;
        newBoundary.erase(*i);
      }
    } while(bPruned);
  }

  // add the old boundary
  newBoundary.insert(boundary.begin(),boundary.end());

  TriangleSet triangleTraversed;
  EdgeSet edgeTraversed;
  SegmentSet segment;

  for(TriangleSet::const_iterator i = triangle.begin(); i != triangle.end(); ++i) {
    Triangle *t = *i;
    if(triangleTraversed.find(t) != triangleTraversed.end()) continue;
    Segment *s = new Segment;
    traverse(s->edge,
             s->boundary,
             s->triangle,
             t,
             triangleTraversed,
             edgeTraversed,
             newBoundary,
             &triangle);
    segment.insert(s);
  }

  bool bSplit = false;
  if(segment.size() > 1) {
    for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
      Segment *s = *i;
      bSplit = true;
      int size = s->edge.size();
      if(size < minSize) {
        bSplit = false;
        break;
      }
      Real r;
      Real a;
      Point sc;
      Vector3 sn;
      segmentInfo(s,&r,&a,&sc,&sn);
      if(r < minRadius) {
        bSplit = false;
        break;
      }
    }
  }

  if(DEBUG) {
    if(bSplit)
      printf("Splitting %d -> [",(int)segment.size());
    else
      printf("Not Splitting %d -> [",(int)segment.size());
    for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
      Segment *s = *i;
      printf("%d ",(int)s->edge.size());
    }
    printf("]\n");
  }
  
  if(bSplit) {
    sub = new RWGTree*[segment.size()];
    subs = 0;
    for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
      Segment *s = *i;
      RWGSet rwgSub;
      for(int j=0; j<rwgs; j++) {
        RWG *r = rwg[j];
        Edge *e = r->e;
        if(s->edge.find(e) != s->edge.end())
          rwgSub.insert(r);
      }
      int size = rwgSub.size();
      if(DEBUG) {
        printf("rwgs: %d->%d\n",rwgs,size);
        printf("boundaries: %d->%d/%d\n",(int)boundary.size(),(int)s->boundary.size(),(int)newBoundary.size());
        printf("triangles: %d->%d\n", (int)triangle.size(), (int)s->triangle.size());
        printf("\n");
      }
      sub[subs++] = new RWGTree(rwgSub,
                                s,
                                edgesAtPoint,
                                diel,
                                g,
                                minSize,
                                minRadius,
                                segment.size()==1);
    }
  } else {
    subs = 0;
    sub = NULL;
  }

  for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
    Segment *s = *i;
    delete s;
  }
  
  if(DEBUG) 
    printf("\n");
}

Real crossSectionBoundaryLength(Real alpha, void *args_)
{
  CrossSectionBoundaryLengthArgs *args = (CrossSectionBoundaryLengthArgs*) args_;

  Point &centroid = args->centroid;
  RWG **rwg = args->rwg;
  int rwgs = args->rwgs;
  PointAndVector *pv = args->pv;
  
  Real phi = pv->p[0] + alpha * pv->xi[0];
  Real theta = pv->p[1] + alpha * pv->xi[1];
  Real psi = pv->p[2] + alpha * pv->xi[2];

  Real sinphi = sin(phi);
  Real cosphi = cos(phi);
  Real sintheta = sin(theta);
  Real costheta = cos(theta);
  Real x = sinphi * costheta;
  Real y = sinphi * sintheta;
  Real z = cosphi;
  Vector3 axis(x,y,z);
  Rotation R(axis,psi);

  Real length = 0.;
  for(int i=0; i<rwgs; i++) {
    RWG *r = rwg[i];
    Edge *e = r->e;
    Point p1 = *(e->p1) - centroid;
    Point Rp1 = R.rotate(p1);
    Real z1 = Rp1.z;
    Point p2 = *(e->p2) - centroid;
    Point Rp2 = R.rotate(p2);
    Real z2 = Rp2.z;

    // if the edge touches the z=0 plane it is a cut
    Real small = 1e-6 * e->length;
    if((z1 <= small && z2 >= -small) || (z1 >= -small && z2 <= small))
      length += sqrt(Square(e->length) - Square(z1-z2));
  }

  return length;
}

void RWGTree :: partition(int minSize, Real minRadius)
{
  PointAndVector pv;
  Real angles[3];
  Real dangles[3];
  pv.p = angles;
  pv.xi = dangles;

  // initialize with the rotation that sends the average normal to the z=0 plane
  Vector3 xy(n.x,n.y,0.);
  Real normxy = sqrt(norm(xy));
  Real psi = (n.z==0. && normxy==0.) ? 0. : atan2(n.z,normxy);
  Vector3 v = n * xy;
  Real normv = sqrt(norm(v));
  Real phi, theta;
  if(normv < 1e-9) {
    phi = PIOVER2;
    theta = 0.;
  } else {
    v *= 1. / normv;
    phi = acos(v.z);
    theta = (v.y == 0. && v.x == 0.) ? 0. : atan2(v.y,v.x);
  }
  Real bestAngles1[3];
  Real bestAngles2[3];
  bestAngles1[0] = phi;
  bestAngles1[1] = theta;
  bestAngles1[2] = psi;

  Real* direction[3];
  direction[0] = new Real[3];
  direction[1] = new Real[3];
  direction[2] = new Real[3];
  direction[0][0] = 1.; direction[0][1] = 0.; direction[0][2] = 0.;
  direction[1][0] = 0.; direction[1][1] = 1.; direction[1][2] = 0.;
  direction[2][0] = 0.; direction[2][1] = 0.; direction[2][2] = 1.;

  CrossSectionBoundaryLengthArgs args(centroid,boundary,edgesAtPoint);
  args.rwg = rwg;
  args.rwgs = rwgs;
  args.pv = &pv;

  int iters;
  Real min1;
  powell(bestAngles1, 
         (Real**)direction,
         3,
         1e-3,
         1e-3,
         &iters,
         &min1,
         &crossSectionBoundaryLength,
         &args,
         &pv);

  // try again with the same axis and a 90 degree rotation
  bestAngles2[0] = bestAngles1[0];
  bestAngles2[1] = bestAngles1[1];
  bestAngles2[2] = bestAngles1[2] + PIOVER2;

  Real min2;
  powell(bestAngles2, 
         (Real**)direction,
         3,
         1e-3,
         1e-3,
         &iters,
         &min2,
         &crossSectionBoundaryLength,
         &args,
         &pv);

  delete [] direction[0];
  delete [] direction[1];
  delete [] direction[2];
  
  Real *bestAngles = (min1 < min2) ? bestAngles1 : bestAngles2;
  partition(bestAngles,minSize,minRadius);
}

RWGTree :: RWGTree(RWGSet &rwgSet,
                   Segment *s,
                   PointEdgeVectorMap &edgesAtPoint_,
                   Dielectric *diel,
                   GreensFunction *g,
                   int minSize,
                   Real minRadius,
                   bool bTerminal) : edgesAtPoint(edgesAtPoint_)
{
  int size = rwgSet.size();
  this->diel = diel;
  this->g = g;
  rwg = new RWG*[size];
  rwgs = 0;
  boundary.insert(s->boundary.begin(),s->boundary.end());
  triangle.insert(s->triangle.begin(),s->triangle.end());
  for(RWGSet::const_iterator i = rwgSet.begin(); i != rwgSet.end(); ++i) {
    RWG *r = *i;
    rwg[rwgs++] = r;
  }
  segmentInfo(s,&radius,&totalArea,&centroid,&n);
  bParameterizable = true;
  if(!bTerminal)
    partition(minSize, minRadius);
  else {
    subs = 0;
    sub = NULL;
  }
}

RWGTree :: RWGTree(int rwgs, PointEdgeVectorMap &edgesAtPoint_, Dielectric *diel, GreensFunction *g) : edgesAtPoint(edgesAtPoint_)
{
  rwg = new RWG*[rwgs];
  this->diel = diel;
  this->g = g;
  this->rwgs = 0;
  this->subs = 0;
  this->sub = NULL;
  this->bParameterizable = false;
}

RWGTree :: RWGTree(ShapeMeshDielectricPairList &meshDiel,
                   PointEdgeVectorMap &edgesAtPoint_,
                   int minSize,
                   Real minRadius)
  : edgesAtPoint(edgesAtPoint_)
{
  bParameterizable = false;

  // count the total number of rwgs, and the number of subtrees
  int tmssize = 0;
  int tsubs = 0;
  int msize = meshDiel.size();
  for(ShapeMeshDielectricPairList::const_iterator j = meshDiel.begin(); j != meshDiel.end(); ++j) {
    ShapeMesh *m = j->first;
    SegmentSet &segment = m->segment;
    if(msize > 1)
      tsubs++;
    for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
      Segment *s = *i;
      tmssize += s->edge.size();
      if(msize == 1)
        tsubs++;
    }
  }

  // allocate the rwg array
  rwg = new RWG*[tmssize];
  rwgs = 0;

  // allocate subtree array if there is more than one subtree
  if(tsubs > 1) {
    sub = new RWGTree*[tsubs];
  }
  subs = 0;
  
  for(ShapeMeshDielectricPairList::const_iterator j = meshDiel.begin(); j != meshDiel.end(); ++j) {

    ShapeMesh *m = j->first;
    Dielectric *diel = j->second;
    GreensFunction *g = new GreensFunction(false);

    // count the rwgs in this mesh
    int mssize = 0;
    // cont the segments in this mesh
    int msubs = 0;
    SegmentSet &segment = m->segment;
    for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
      Segment *s = *i;
      EdgeSet &edge = s->edge;
      mssize += edge.size();
      msubs++;
    }

    // if there is more than one mesh, allocate a new subtree for this mesh
    // and add it to the subtrees
    RWGTree *tm;
    if(msize > 1) {
      tm = new RWGTree(mssize,edgesAtPoint,diel,g);
      sub[subs++] = tm;
      // and if there is more than one segment in this mesh, allocate a subtree array
      if(msubs > 1) {
        tm->sub = new RWGTree*[msubs];
        tm->subs = 0;
      }
    } else {
      this->diel = diel;
      this->g = g;
      tm = this;
    }

    for(SegmentSet::const_iterator i = segment.begin(); i != segment.end(); ++i) {
      Segment *s = *i;

      // count the rwgs in this segment
      EdgeSet &edge = s->edge;
      int ssize = edge.size();

      // if there is more than one segment in this mesh, allocate a new subtree for this segment
      // and add it to the mesh subtrees
      RWGTree *ts;
      if(msubs > 1) {
        ts = new RWGTree(ssize,edgesAtPoint,diel,g);
        tm->sub[tm->subs++] = ts;
      } else {
        ts = tm;
      }

      // construct the rwgs in this segment
      for(EdgeSet::const_iterator k = edge.begin(); k != edge.end(); ++k) {
        Edge *e = *k;
        RWG *r = new RWG(e);

        // if there is more than one segment in this mesh, add the rwg to its subtree
        if(msubs > 1) {
          ts->rwg[ts->rwgs] = r;
          ts->rwgs++;
        }

        // if there is more than one mesh, add the rwg to its subtree
        if(msize > 1) {
          tm->rwg[tm->rwgs] = r;
          tm->rwgs++;
        }
        
        rwg[rwgs] = r;
        rwgs++;
      }
      
      // add segment boundaries/triangles
      ts->boundary.insert(s->boundary.begin(),s->boundary.end());
      ts->triangle.insert(s->triangle.begin(),s->triangle.end());

      // amend segment info
      segmentInfo(s,&ts->radius,&ts->totalArea,&ts->centroid,&ts->n);

      // partition the segment
      ts->partition(minSize, minRadius);
    }
  }
  unshuffle();
}

void RWGTree :: output(ostream& os, int *k)
{
  if(!subs) {
    os << "s{" << (*k) << "} = [\n";
    for(TriangleSet::const_iterator i = triangle.begin(); i != triangle.end(); ++i) {
      Triangle *t = *i;
      os << *t << ";\n";
      //os << "[" << center(t) << " " << normal(t) << "];\n";
    }
    os << "];\n";
    os << "b{" << (*k)++ << "} = [\n";
    for(EdgeSet::const_iterator i = boundary.begin(); i != boundary.end(); ++i) {
      Edge *e = *i;
      os << "[" << *(e->p1) << " " << *(e->p2) << "];\n";
    }
    os << "];\n";
  } else {
    for(int i=0; i<subs; i++) {
      RWGTree *tsub = sub[i];
      tsub->output(os,k);
    }
  }
}

void RWGTree :: unshuffle()
{
  if(!subs) return;
  int k = 0;
  for(int i=0; i<subs; i++) {
    RWGTree *s = sub[i];
    s->unshuffle();
    for(int j=0; j<s->rwgs; j++) {
      rwg[k] = s->rwg[j];
      s->index.push_back(k);
      k++;
    }
  }
}

ostream& operator<<(ostream& os, RWGTree *t)
{
  int k = 1;
  t->output(os, &k);
  return os;
}

#undef DEBUG
