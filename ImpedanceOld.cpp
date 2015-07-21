#include "Impedance.h"
#include "GreensFunction.h"
#include "FunctionArgs.h"
#include "TriangleQuadratureOld.h"
#include "AnalyticIntegrals.h"
#include "Triangle.h"
#include "Options.h"
#include "MatrixOperations.h"
#include "ComplexVector3.h"

//#define DEBUG
#undef DEBUG

#define BORDEREPS 1e-12

void avoidBorders(Real v2x, Real v2y, Real *x0, Real *y0)
{
  Real x,y;
  y = *y0 / v2y;
  x = *x0 - v2x * y;
  while(fabs(x) < BORDEREPS || fabs(y) < BORDEREPS || fabs(1. - x - y) < BORDEREPS) {
    if(fabs(x) < BORDEREPS) {
      *y0 += BORDEREPS * v2x;
      *x0 += BORDEREPS * v2y;
    } else if(fabs(y) < BORDEREPS) {
      *y0 += BORDEREPS;
    } else {
      *x0 += BORDEREPS * v2y;
      *y0 += BORDEREPS * (v2x - 1.);
    }
    y = *y0 / v2y;
    x = *x0 - v2x * y;
  }
  if(fabs(v2x - 1.) < BORDEREPS && fabs(*x0 - 1.) < BORDEREPS) {
    *x0 = (1. + BORDEREPS);
  }
}

void integrand_G3s(Vector3 *f, Vector3 *x, FunctionArgs *args)
{
  Vector3 rs = args->sh.v0 + (*f);
  Vector3 r = (*args->rt - rs);
  Real R = sqrt(norm(r));
  *x = -ONEOVER4PI / (R*R*R) * r;
}

void integrand_G3t(Vector3 *g, Real *x, FunctionArgs *args)
{
  Vector3 rt = args->th.v0 + (*g);
  args->rt = &rt;
  QuadratureFunction<Vector3,VectorValuedFunction>::Type quadfunc = &quad5553<Vector3,VectorValuedFunction>;
  Vector3 Lv = quadfunc(integrand_G3s, args, args->sh.v1, args->sh.v2, NULL, NULL, args->tol, 0.);
  *x = -2. * dot(*g,Lv);
}

void integrand_LK_nonsingular(Vector3 *f, ComplexVector3Pair *x, FunctionArgs *args)
{
  Vector3 rs = args->tti->v0s + (*f);
  Vector3 f2 = rs - args->sh.v0;
  ComplexVector3 gg;
  Complex g;
  if(args->tti->bSet) {
    args->tti->getScalar(g);
    args->tti->getGradient(gg);
    args->tti->next();
  } else {
    args->green->scalar_gradient_nonsingular(args->rt,&rs,&g,&gg);
    args->tti->setScalar(g);
    args->tti->setGradient(gg);
    args->tti->next();
  }
  
  x->first = (args->green->ksquared*g)*f2+2.*gg;
  x->second = gg*f2;
}

void integrand_LK(Vector3 *f, ComplexVector3Pair *x, FunctionArgs *args)
{
  Vector3 rs = args->tti->v0s + (*f);
  Vector3 f2 = rs - args->sh.v0;
  ComplexVector3 gg;
  Complex g;
  if(args->tti->bSet) {
    args->tti->getScalar(g);
    args->tti->getGradient(gg);
    args->tti->next();
  } else {
    args->green->scalar_gradient(args->rt,&rs,&g,&gg);
    args->tti->setScalar(g);
    args->tti->setGradient(gg);
    args->tti->next();
  }
  
  x->first = (args->green->ksquared*g)*f2+2.*gg;
  x->second = gg*f2;
}

void integrand_LK_FarTriangles(Vector3 *g, ComplexPair *x, FunctionArgs *args)
{
  Vector3 rt = args->tti->v0t + (*g);
  args->rt = &rt;
  Vector3 g2 = rt - args->th.v0;
  ComplexVector3Pair lk = args->squadfunc(integrand_LK, args, args->tti->v1s, args->tti->v2s, NULL, NULL, args->tol, 0.);
  *x = ComplexPair(dot(g2,lk.first),dot(g2,lk.second));
}

void integrand_LK1(Vector3 *g, ComplexVector3 *x, FunctionArgs *args)
{
  Vector3 rt = args->tpi->v0 + (*g);
  Vector3 g2 = rt - args->th.v0;
  Complex gs;
  if(args->tpi->bSet) {
    args->tpi->getScalar(gs);
    args->tpi->next();
  } else {
    gs = args->green->scalar(&rt,args->rs);
    args->tpi->setScalar(gs);
    args->tpi->next();
  }
  *x = gs * g2;
}

void integrand_LK2(Vector3 *f, ComplexVector3Pair *x, FunctionArgs *args)
{
  Vector3 rs = args->tpi->v0 + (*f);
  Vector3 f2 = rs - args->sh.v0;
  ComplexVector3 gg;
  Complex gs;
  if(args->tpi->bSet) {
    args->tpi->getScalar(gs);
    args->tpi->getGradient(gg);
    args->tpi->next();
  } else {
    args->green->scalar_gradient(args->rt,&rs,&gs,&gg);
    args->tpi->setScalar(gs);
    args->tpi->setGradient(gg);
    args->tpi->next();
  }
  
  x->first = (args->green->ksquared*gs)*f2+2.*gg;
  x->second = gg*f2;
}

// leave out G3 and C3
void integrand_LK_NeighboringTriangles(Vector3 *g, ComplexPair *x, FunctionArgs *args)
{
  // get source position relative to s origin
  Vector3 gr = args->th.v0 + (*g) - args->sh.v0;

  // collapse to 2D coords
  Vector3 &v1 = args->sh.v1;
  Vector3 &v2 = args->sh.v2;
  Vector3 &n = args->sh.n;
  Real normv1i = args->sh.normv1i;
  Real normv2i = args->sh.normv2i;
  Real normv1i2 = normv1i * normv1i;
  Real normv2i2 = normv2i * normv2i;

  Real z00 = dot(gr,n);
  Real z0 = fabs(z00) * normv1i;
  Real yz0 = fabs(z00) * normv2i;

  // For x
  Real x0 = normv1i2 * dot(gr,v1);
  Real y0 = normv1i2 * dot(gr,n * v1);

  // For y
  Real yx0 = normv2i2 * dot(gr,n * v2);
  Real yy0 = normv2i2 * dot(gr,v2);

  Real v1x = args->sh.v1x;
  Real v1y = args->sh.v1y;
  Real v2x = args->sh.v2x;
  Real v2y = args->sh.v2y;

  avoidBorders(v2x,v2y,&x0,&y0);
  avoidBorders(v1y,v1x,&yy0,&yx0);
 
  // The Grunge
  Real S0, S1, G0x, G0y, G1x, G1y, C0z, C1z;
  if(args->bzSmall) {
    z0_S0_S1_G0_G1(&S0, &S1, &G0x, &G0y, &G1x, &G1y, v2x, v2y, x0, y0, v1x, v1y, yx0, yy0);
    C0z = 0.;
    C1z = 0.;
    z00 = 0.;
  } else {
    S0_S1_G0_G1_C0_C1(&S0, &S1,
                      &G0x, &G0y, &G1x, &G1y,
                      &C0z, &C1z,
                      v2x, v2y, x0, y0, z0,
                      v1x, v1y, yx0, yy0, yz0);
  }
  //printf("%g %g %g %g %g %g %g %g\n",S0,S1,G0x,G0y,G1x,G1y,C0z,C1z);
    
  // get the y values in v1 coords
  Real n21 = normv2i / normv1i;
  G0y = G0y * v2y - v2x*(G0x - G0y*v2x)/v2y;
  G1y *= n21;
  G1y = G1y * v2y - v2x*(G1x - G1y*v2x)/v2y;

  S1 *= normv1i;
  G1x *= normv1i;
  G1y *= normv1i;
  C0z /= -normv1i2;
  C1z /= -normv1i;
  
  Complex& ksquared = args->green->ksquared;
  Complex& ik = args->green->ik;
  Complex k2 = -0.5 * ksquared;

  Complex xLscalar = (x0 * S1 - G1x) + ik * (x0 * S0 - G0x);
  Complex xLgradient = k2 * (ik * G0x + G1x);
  Complex xL = +ksquared * xLscalar + 2.0 * xLgradient;
  Complex yLscalar = (y0 * S1 - G1y) + ik * (y0 * S0 - G0y);
  Complex yLgradient = k2 * (ik * G0y + G1y);
  Complex yL = +ksquared * yLscalar + 2.0 * yLgradient;
  Complex zLscalar = 0.;
  Complex zLgradient = z00 * (k2 * (ik * S0 + S1));
  Complex zL = +ksquared * zLscalar + 2.0 * zLgradient;

  Complex xK = z00 * (k2 * (ik * (G0y - y0 * S0) + (G1y - y0 * S1)));
  Complex yK = -z00 * (k2 * (ik * (G0x - x0 * S0) + (G1x - x0 * S1)));
  Complex zK = k2 * (ik * C0z + C1z);

  ComplexVector3 Lv = expandCoords(args->sh.v1,args->sh.n,xL,yL,zL);
  ComplexVector3 Kv = expandCoords(args->sh.v1,args->sh.n,xK,yK,zK);

  Complex L = dot(*g,Lv);
  Complex K = dot(*g,Kv);

  *x = ONEOVER4PI * ComplexPair(L,K);
}


void integrand_LK_NeighboringTriangles_nonsingular(Vector3 *g, ComplexPair *x, FunctionArgs *args)
{
  Vector3 rt = args->tti->v0t + (*g);
  args->rt = &rt;
  Vector3 g2 = rt - args->th.v0;
  ComplexVector3Pair lk = args->squadfunc(integrand_LK_nonsingular, args, args->tti->v1s, args->tti->v2s, NULL, NULL, args->tol, 0.);
  x->first = dot(g2,lk.first);
  x->second = dot(g2,lk.second);
}

void integrand_LK_NeighboringTriangles_G3(Vector3 *g, Real *x, FunctionArgs *args)
{
  // get source position relative to s origin
  Vector3 gr = args->th.v0 + (*g) - args->sh.v0;

  // collapse to 2D coords
  Vector3 &v1 = args->sh.v1;
  Vector3 &v2 = args->sh.v2;
  Vector3 &n = args->sh.n;
  Real normv1i = args->sh.normv1i;
  Real normv2i = args->sh.normv2i;
  Real normv1i2 = normv1i * normv1i;
  Real normv2i2 = normv2i * normv2i;

  Real z00 = dot(gr,n);
  Real z0 = fabs(z00) * normv1i;
  Real yz0 = fabs(z00) * normv2i;

  // For x
  Real x0 = normv1i2 * dot(gr,v1);
  Real y0 = normv1i2 * dot(gr,n * v1);

  // For y
  Real yx0 = normv2i2 * dot(gr,n * v2);
  Real yy0 = normv2i2 * dot(gr,v2);

  Real v1x = args->sh.v1x;
  Real v1y = args->sh.v1y;
  Real v2x = args->sh.v2x;
  Real v2y = args->sh.v2y;

  avoidBorders(v2x,v2y,&x0,&y0);
  avoidBorders(v1y,v1x,&yy0,&yx0);

  // The Grunge
  Real S3, G3x, G3y;
  if(args->bzSmall) {
    G3(&G3x, &G3y,
       v2x, v2y, x0, y0, z0,
       v1x, v1y, yx0, yy0, yz0);
    z00 = 0.;
    S3 = 0.;
  } else {
    S3_G3(&S3, &G3x, &G3y,
          v2x, v2y, x0, y0, z0,
          v1x, v1y, yx0, yy0, yz0);
  }
    
  // get the y values in v1 coords
  Real n21 = normv2i / normv1i;
  G3y *= n21 * n21 * n21;
  G3y = G3y * v2y - v2x*(G3x - G3y*v2x)/v2y;

  // scale the values
  Real normv1i3 = normv1i2 * normv1i;

  S3 *= normv1i3;
  G3x *= normv1i3;
  G3y *= normv1i3;
  
  Real xL = -2. * G3x;
  Real yL = -2. * G3y;
  Real zL = -2. * z00 * S3;

  Vector3 Lv = expandCoords(args->sh.v1,args->sh.n,xL,yL,zL);

  *x = ONEOVER4PI * dot(*g,Lv);
}

void NeighboringTriangles_G3overlog_edgelimit(Vector3 *g, Real *x, FunctionArgs *args)
{
  // get source position relative to s origin
  Vector3 gr = args->th.v0 + (*g) - args->sh.v0;

  // collapse to 2D coords
  Vector3 &v1 = args->sh.v1;
  Vector3 &n = args->sh.n;
  Real normv1i = args->sh.normv1i;
  Real normv2i = args->sh.normv2i;
  Real normv1i2 = normv1i * normv1i;

  // For x
  Real x0 = normv1i2 * dot(gr,v1);
  Real y0 = normv1i2 * dot(gr,n * v1);

  Real v1x = args->sh.v1x;
  Real v1y = args->sh.v1y;
  Real v2x = args->sh.v2x;
  Real v2y = args->sh.v2y;

  // The Grunge
  Real G3x, G3y;
  G3overlog_edgelimit(&G3x, &G3y,
                      v2x, v2y,
                      x0, y0,
                      v1x, v1y);

  // get the y values in v1 coords
  Real n21 = normv2i / normv1i;
  G3y *= n21 * n21 * n21;
  G3y = G3y * v2y - v2x*(G3x - G3y*v2x)/v2y;

  // scale the values
  Real normv1i3 = normv1i2 * normv1i;

  G3x *= normv1i3;
  G3y *= normv1i3;
  
  Real xL = -2. * G3x;
  Real yL = -2. * G3y;
  Real zL = 0.;

  Vector3 Lv = expandCoords(args->sh.v1,args->sh.n,xL,yL,zL);

  *x = ONEOVER4PI * dot(*g,Lv);
}

void integrand_LK_NeighboringTriangles_C3(Vector3 *g, Real *x, FunctionArgs *args)
{
  // get source position relative to s origin
  Vector3 gr = args->th.v0 + (*g) - args->sh.v0;

  // collapse to 2D coords
  Vector3 &v1 = args->sh.v1;
  Vector3 &v2 = args->sh.v2;
  Vector3 &n = args->sh.n;
  Real normv1i = args->sh.normv1i;
  Real normv2i = args->sh.normv2i;
  Real normv1i2 = normv1i * normv1i;
  Real normv2i2 = normv2i * normv2i;

  Real z00 = dot(gr,n);
  Real z0 = fabs(z00) * normv1i;
  Real yz0 = fabs(z00) * normv2i;

  // For x
  Real x0 = normv1i2 * dot(gr,v1);
  Real y0 = normv1i2 * dot(gr,n * v1);

  // For y
  Real yx0 = normv2i2 * dot(gr,n * v2);
  Real yy0 = normv2i2 * dot(gr,v2);

  Real v1x = args->sh.v1x;
  Real v1y = args->sh.v1y;
  Real v2x = args->sh.v2x;
  Real v2y = args->sh.v2y;

  avoidBorders(v2x,v2y,&x0,&y0);
  avoidBorders(v1y,v1x,&yy0,&yx0);

  // The Grunge
  Real S3, G3x, G3y, C3z;

  if(args->bzSmall) {
    C3(&C3z,v2x, v2y, x0, y0, z0);
    S3 = 0.;
    G3x = 0.;
    G3y = 0.;
    z00 = 0.;
  } else {
    S3_G3_C3(&S3, &G3x, &G3y, &C3z,
             v2x, v2y, x0, y0, z0,
             v1x, v1y, yx0, yy0, yz0);
  }

  // get the y values in v1 coords
  Real n21 = normv2i / normv1i;
  G3y *= n21 * n21 * n21;
  G3y = G3y * v2y - v2x*(G3x - G3y*v2x)/v2y;

  // scale the values
  Real normv1i3 = normv1i2 * normv1i;

  S3 *= normv1i3;
  G3x *= normv1i3;
  G3y *= normv1i3;
  C3z *= -normv1i;

  Real xK = -z00 * (G3y - y0 * S3);
  Real yK = z00 * (G3x - x0 * S3);
  Real zK = -C3z;

  Vector3 Kv = expandCoords(args->sh.v1,args->sh.n,xK,yK,zK);

  *x = ONEOVER4PI * dot(*g,Kv);
}

void NeighboringTriangles_C3overlog_edgelimit(Vector3 *g, Real *x, FunctionArgs *args)
{
  // get source position relative to s origin
  Vector3 gr = args->th.v0 + (*g) - args->sh.v0;

  // collapse to 2D coords
  Vector3 &v1 = args->sh.v1;
  Vector3 &n = args->sh.n;
  Real normv1i = args->sh.normv1i;
  Real normv1i2 = normv1i * normv1i;

  // For x
  Real x0 = normv1i2 * dot(gr,v1);
  Real y0 = normv1i2 * dot(gr,n * v1);

  Real v2x = args->sh.v2x;
  Real v2y = args->sh.v2y;
  
  // The Grunge
  Real C3z;
  C3overlog_edgelimit(&C3z,
                      v2x, v2y, x0, y0);

  C3z *= -normv1i;

  Real xK = 0.;
  Real yK = 0.;
  Real zK = -C3z;

  Vector3 Kv = expandCoords(args->sh.v1,args->sh.n,xK,yK,zK);

  *x = ONEOVER4PI * dot(*g,Kv);
}

// G3 is left out 
void integrand_LK_SameTriangles(Vector3 *g, ComplexPair *x, FunctionArgs *args)
{
  // get source position relative to s origin
  Vector3 gr = args->th.v0 + (*g) - args->sh.v0;
  
  // collapse to 2D coords
  Vector3 &v1 = args->sh.v1;
  Vector3 &v2 = args->sh.v2;
  Vector3 &n = args->sh.n;
  Real normv1i = args->sh.normv1i;
  Real normv2i = args->sh.normv2i;
  Real normv1i2 = normv1i * normv1i;
  Real normv2i2 = normv2i * normv2i;

  // For x
  Real x0 = normv1i2 * dot(gr,v1);
  Real y0 = normv1i2 * dot(gr,n * v1);

  // For y
  Real yx0 = normv2i2 * dot(gr,n * v2);
  Real yy0 = normv2i2 * dot(gr,v2);

  Real v1x = args->sh.v1x;
  Real v1y = args->sh.v1y;
  Real v2x = args->sh.v2x;
  Real v2y = args->sh.v2y;

  // The Grunge
  Real S0, S1, G0x, G0y, G1x, G1y;
  z0_S0_S1_G0_G1(&S0, &S1, &G0x, &G0y, &G1x, &G1y, v2x, v2y, x0, y0, v1x, v1y, yx0, yy0);

  // get the y values in v1 coords
  Real n21 = normv2i / normv1i;
  G0y = G0y * v2y - v2x*(G0x - G0y*v2x)/v2y;
  G1y *= n21;
  G1y = G1y * v2y - v2x*(G1x - G1y*v2x)/v2y;

  // scale the values
  S1 *= normv1i;
  G1x *= normv1i;
  G1y *= normv1i;

  Complex& ksquared = args->green->ksquared;
  Complex& ik = args->green->ik;
  Complex k2 = -0.5 * ksquared;

  Complex xLscalar = (x0 * S1 - G1x) + ik * (x0 * S0 - G0x);
  Complex xLgradient = k2 * (ik * G0x + G1x);
  Complex xL = +ksquared * xLscalar + 2.0 * xLgradient;
  Complex yLscalar = (y0 * S1 - G1y) + ik * (y0 * S0 - G0y);
  Complex yLgradient = k2 * (ik * G0y + G1y);
  Complex yL = +ksquared * yLscalar + 2.0 * yLgradient;
  ComplexVector3 Lv = expandCoords(args->sh.v1,args->sh.n,xL,yL);

  Vector3 Kv = args->sh.n * gr;

  Complex L = ONEOVER4PI * dot(*g,Lv);
  Complex K = (args->green->isExternal()?0.5:-0.5) / args->sh.A2 * dot(*g,Kv);

  x->first = L;
  x->second = K;
}

void integrand_LK_SameTriangles_nonsingular(Vector3 *g, ComplexPair *x, FunctionArgs *args)
{
  Vector3 rt = args->tti->v0t + (*g);
  args->rt = &rt;
  Vector3 g2 = rt - args->th.v0;
  ComplexVector3Pair lk = args->squadfunc(integrand_LK_nonsingular, args, args->tti->v1s, args->tti->v2s, NULL, NULL, args->tol, 0.);
  x->first = dot(g2,lk.first);
  x->second = Complex(0.,0.);
}

Real gdotG3(FunctionArgs *args)
{
  Real normv1i = args->sh.normv1i;
  return ONEOVER4PI * -2.0  * normv1i * z0_gdotG3(args->sh.v2x,args->sh.v2y);
}

ComplexPair LK(RWG *t, RWG *s, FunctionArgs *args, bool btPositive, bool bsPositive, TriangleTriangleInteractionCache &ttiCache, TriangleTriangleInteractionCount &ttiCount, Options &options)
{
  QuadratureFunction<ComplexPair,ComplexPairValuedFunction>::Type tquadfunc;

  Triangle *trit = btPositive?t->t1:t->t2;
  Triangle *tris = bsPositive?s->t1:s->t2;

  Key key = getKey(trit->index,tris->index);
  TriangleTriangleInteractionCache::iterator ttii = ttiCache.find(key);
  if(ttii != ttiCache.end()) {
    args->tti = ttii->second;
  } else {
    args->tti = NULL;
  }  
  int useTotal = ttiCount[key];
  
  ComplexPair lk;

  if(trit == tris) {
#ifdef DEBUG
    printf("Same\n");
#endif
    getTriangle(t,&(args->th),1.,btPositive);
    getTriangle(s,&(args->sh),1.,bsPositive);

    tquadfunc = &quadgauss105<ComplexPair,ComplexPairValuedFunction>;
    lk = tquadfunc(integrand_LK_SameTriangles, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
  
    if(args->tti == NULL) {
      args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,91,105,useTotal);
      ttiCache[key] = args->tti;
    }
    args->tti->start();
    tquadfunc = &quadgauss91<ComplexPair,ComplexPairValuedFunction>;
    args->squadfunc = &quadgauss105<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
    lk += tquadfunc(integrand_LK_SameTriangles_nonsingular, args, args->tti->v1t, args->tti->v2t, NULL, NULL, args->tol, 0.);
    args->tti->end();

    Real LG3 = gdotG3(args);
    
    lk.first += LG3;
#ifdef DEBUG
    cout << LG3 << "\n";
#endif
  } else {
    int shared = sharedpoints(trit,tris);
    if(shared == 2) {
      Point *tv0 = btPositive?t->v01:t->v02;
      if(tv0 != tris->p1 && tv0 != tris->p2 && tv0 != tris->p3) {
        
#ifdef DEBUG 
        printf("Edge Neighbors Singular\n");
#endif
        // the moderately singular contribution from the full triangles
        getTriangle(t,&(args->th),1.,btPositive);
        getTriangle(s,&(args->sh),1.,bsPositive);
        Real A2 = sqrt(norm(args->th.v1 * args->th.v2));
        
        Vector3 dv = center(trit) - args->sh.v0;
        Real z00 = dot(dv,args->sh.n);
        Real z0 = min(fabs(z00*args->sh.normv1i),fabs(z00*args->sh.normv2i));
        args->bzSmall = z0<=options.inPlaneThresh;

        tquadfunc = &quadgauss105<ComplexPair,ComplexPairValuedFunction>;
        lk = tquadfunc(integrand_LK_NeighboringTriangles, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        
        // the nonsingular contribution from the full triangles
        if(args->tti == NULL) {
          args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,36,36,useTotal);
          ttiCache[key] = args->tti;
        }
        args->tti->start();
        tquadfunc = &quadgauss36<ComplexPair,ComplexPairValuedFunction>;
        args->squadfunc = &quadgauss36<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
        lk += tquadfunc(integrand_LK_NeighboringTriangles_nonsingular, args, args->tti->v1t, args->tti->v2t, NULL, NULL, args->tol, 0.);
        args->tti->end();

        // the highly singular G3 contribution        
        Real l,k;
        Real alpha = 0.4;
        
        // leave out the trapezoid at the edge
        getTriangle(t,&(args->th),alpha,btPositive);
        Real A2trunc = sqrt(norm(args->th.v1 * args->th.v2));
        Real c = A2trunc / A2;
        QuadratureFunction<Real,RealValuedFunction>::Type quadfunc = &quadgauss105<Real,RealValuedFunction>;
        l = c * quadfunc(integrand_LK_NeighboringTriangles_G3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        if(args->bzSmall) {
          k = 0.;
        } else {
          k = c * quadfunc(integrand_LK_NeighboringTriangles_C3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        }
        
        
        // now do the trapezoid
        getTriangle(t,&(args->th),1.,btPositive);
        QuadratureFunction<Real,RealValuedFunction>::Type logquadfunc = &quad7774log<Real,RealValuedFunction>;
        AdaptiveQuadLogArgs<Real,RealValuedFunction> quadargs;
        quadargs.limitfunc = &NeighboringTriangles_G3overlog_edgelimit;
        quadargs.alpha = alpha;
        Real minnorm = min(norm(args->green->eps),norm(args->green->mu)) * norm(args->green->omega);
        Real abstol = max(args->tol * norm(lk.first),options.abstol*minnorm);
        l += logquadfunc(integrand_LK_NeighboringTriangles_G3, args, args->th.v1, args->th.v2, NULL, &quadargs, args->tol, abstol);
        if(args->bzSmall) {
        } else {
          quadargs.limitfunc = &NeighboringTriangles_C3overlog_edgelimit;
          abstol = max(args->tol * norm(lk.second),options.abstol);
          k += logquadfunc(integrand_LK_NeighboringTriangles_C3, args, args->th.v1, args->th.v2, NULL, &quadargs, args->tol, abstol);
        }

        //printf("%g %g\n",l,k);
#ifdef DEBUG
        cout << lk << "\n";
#endif
        lk += ComplexPair(Complex(l,0.),Complex(k,0.));
#ifdef DEBUG
        cout << lk << "\n";
#endif
      } else {
      
#ifdef DEBUG  
        printf("Edge Neighbors Non-singular\n");
#endif
        getTriangle(t,&(args->th),1.,btPositive);
        getTriangle(s,&(args->sh),1.,bsPositive);
        
        Vector3 dv = center(trit) - args->sh.v0;
        Real z00 = dot(dv,args->sh.n);
        Real z0 = min(fabs(z00*args->sh.normv1i),fabs(z00*args->sh.normv2i));
        args->bzSmall = z0<=options.inPlaneThresh;

        tquadfunc = &quadgauss105<ComplexPair,ComplexPairValuedFunction>;
        lk = tquadfunc(integrand_LK_NeighboringTriangles, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);

        if(args->tti == NULL) {
          args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,36,36,useTotal);
          ttiCache[key] = args->tti;
        }
        args->tti->start();
        tquadfunc = &quadgauss36<ComplexPair,ComplexPairValuedFunction>;
        args->squadfunc = &quadgauss36<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
        lk += tquadfunc(integrand_LK_NeighboringTriangles_nonsingular, args, args->tti->v1t, args->tti->v2t, NULL, NULL, args->tol, 0.);
        args->tti->end();

        QuadratureFunction<Real,RealValuedFunction>::Type quadfunc = &quadgauss105<Real,RealValuedFunction>;
        Real l = quadfunc(integrand_LK_NeighboringTriangles_G3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        Real k;
        if(args->bzSmall) {
          k = 0.;
        } else {
          k = quadfunc(integrand_LK_NeighboringTriangles_C3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        }

        /*
        getTriangle(t,&(args->th),1-1e-4,btPositive);
        getTriangle(s,&(args->sh),1-1e-4,bsPositive);

        args->th.v0 += 1e-5 * (center(trit) - args->th.v0);
        args->sh.v0 += 1e-5 * (center(tris) - args->sh.v0);

        quadfunc = &quad5553<Real,RealValuedFunction>;
        l = quadfunc(integrand_LK_NeighboringTriangles_G3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        if(args->bzSmall) {
          k = 0.;
        } else {
          k = quadfunc(integrand_LK_NeighboringTriangles_C3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
        }
        
        //printf("%g %g %g %g\n",l,l2,k,k2);
        */

#ifdef DEBUG
        cout << lk << "\n";
#endif
        lk += ComplexPair(Complex(l,0.),Complex(k,0.));
#ifdef DEBUG
        cout << lk << "\n";
#endif
      }
    } else if(shared == 1) {
#ifdef DEBUG
      printf("Corner Neighbors\n");
#endif
      getTriangle(t,&(args->th),1.,btPositive);
      getTriangle(s,&(args->sh),1.,bsPositive);          

      Vector3 dv = center(trit) - args->sh.v0;
      Real z00 = dot(dv,args->sh.n);
      Real z0 = min(fabs(z00*args->sh.normv1i),fabs(z00*args->sh.normv2i));
      args->bzSmall = z0<=options.inPlaneThresh;

      tquadfunc = &quadgauss105<ComplexPair,ComplexPairValuedFunction>;
      lk = tquadfunc(integrand_LK_NeighboringTriangles, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);

#ifdef DEBUG
      cout << lk << "\n";
#endif
      if(args->tti == NULL) {
        args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,36,36,useTotal);
        ttiCache[key] = args->tti;
      }
      args->tti->start();
      tquadfunc = &quadgauss36<ComplexPair,ComplexPairValuedFunction>;
      args->squadfunc = &quadgauss36<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
      lk += tquadfunc(integrand_LK_NeighboringTriangles_nonsingular, args, args->tti->v1t, args->tti->v2t, NULL, NULL, args->tol, 0.);
      args->tti->end();
      
      QuadratureFunction<Real,RealValuedFunction>::Type quadfunc = &quadgauss105<Real,RealValuedFunction>;
      Real l = quadfunc(integrand_LK_NeighboringTriangles_G3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
      Real k;
      if(args->bzSmall) {
        k = 0.;
      } else {
        k = quadfunc(integrand_LK_NeighboringTriangles_C3, args, args->th.v1, args->th.v2, NULL, NULL, args->tol, 0.);
      }
      
#ifdef DEBUG
      cout << lk << "\n";
#endif
      lk += ComplexPair(Complex(l,0.),Complex(k,0.));
#ifdef DEBUG
      cout << lk << "\n";
#endif
    } else {
#ifdef DEBUG
      printf("Far\n");
#endif
      getTriangle(t,&(args->th),1.,btPositive);
      getTriangle(s,&(args->sh),1.,bsPositive);          

      Point pt = center(trit);
      Point ps = center(tris);
      Vector3 v = pt - ps;
      Real d = sqrt(norm(v));
      Real temax = max(max(trit->e12->length,trit->e23->length),trit->e31->length);
      Real semax = max(max(tris->e12->length,tris->e23->length),tris->e31->length);

      if((args->tti && args->tti->nt == 10) || (d * options.farThresh10 > (temax + semax))) {
        if(args->tti == NULL) {
          args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,10,10,useTotal);
          ttiCache[key] = args->tti;
        }
        tquadfunc = &quadgauss10<ComplexPair,ComplexPairValuedFunction>;
        args->squadfunc = &quadgauss10<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
      } else if((args->tti && args->tti->nt == 15) || (d * options.farThresh15 > (temax + semax))) {
        if(args->tti == NULL) {
          args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,15,15,useTotal);
          ttiCache[key] = args->tti;
        }
        tquadfunc = &quadgauss15<ComplexPair,ComplexPairValuedFunction>;
        args->squadfunc = &quadgauss15<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
      } else if((args->tti && args->tti->nt == 21) || (d * options.farThresh21 > (temax + semax))) {
        if(args->tti == NULL) {
          args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,21,21,useTotal);
          ttiCache[key] = args->tti;
        }
        tquadfunc = &quadgauss21<ComplexPair,ComplexPairValuedFunction>;
        args->squadfunc = &quadgauss21<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
      } else {
        if(args->tti == NULL) {
          args->tti = new TriangleTriangleInteraction(args->th.v0,args->th.v1,args->th.v2,args->sh.v0,args->sh.v1,args->sh.v2,36,36,useTotal);
          ttiCache[key] = args->tti;
        }
        tquadfunc = &quadgauss36<ComplexPair,ComplexPairValuedFunction>;
        args->squadfunc = &quadgauss36<ComplexVector3Pair,ComplexVectorPairValuedFunction>;
      }
      args->tti->start();
      lk = tquadfunc(integrand_LK_FarTriangles, args, args->tti->v1t, args->tti->v2t, NULL, NULL, args->tol, 0.);
      args->tti->end();
        
#ifdef DEBUG      
      cout << lk << "\n";
#endif
    }
  }

  if(args->tti->useCount == args->tti->useTotal) {
    ttiCache.erase(key);
    ttiCount.erase(key);
    delete args->tti;
  }

  return lk;
}

void countTriangleTriangleInteractions(RWG *t, RWG *s, bool btPositive, bool bsPositive, TriangleTriangleInteractionCount &ttiCount)
{
  Triangle *trit = btPositive?t->t1:t->t2;
  Triangle *tris = bsPositive?s->t1:s->t2;

  Key key = getKey(trit->index, tris->index);
  TriangleTriangleInteractionCount::iterator ttii = ttiCount.find(key);
  int count;
  if(ttii != ttiCount.end()) {
    count = ttii->second + 1;
  } else {
    count = 1;
  }
  ttiCount[key] = count;
}

void countTriangleTriangleInteractions(RWG *t, RWG *s, TriangleTriangleInteractionCount &ttiCount)
{
  countTriangleTriangleInteractions(t,s,true,true,ttiCount);
  countTriangleTriangleInteractions(t,s,false,true,ttiCount);
  countTriangleTriangleInteractions(t,s,false,false,ttiCount);
  countTriangleTriangleInteractions(t,s,true,false,ttiCount);
}


ComplexPair LK(RWG *t, RWG *s, GreensFunction *g, TriangleTriangleInteractionCache &ttiCache, TriangleTriangleInteractionCount &ttiCount, Options &options)
{
  FunctionArgs args;
  args.green = g;
  args.tol = options.tol;

  ComplexPair lk;
  lk  = LK(t,s,&args,true,true,ttiCache,ttiCount,options);
  lk -= LK(t,s,&args,false,true,ttiCache,ttiCount,options);
  lk += LK(t,s,&args,false,false,ttiCache,ttiCount,options);
  lk -= LK(t,s,&args,true,false,ttiCache,ttiCount,options);

  lk *= t->e->length * s->e->length;

  return lk;
}

ComplexVector3 LK1(RWG *t, int s, FunctionArgs *args, TrianglePointInteractionCache &tpiCache, TrianglePointInteractionCount &tpiCount, bool btPositive)
{
  getTriangle(t,&(args->th),1.,btPositive);
  Triangle *trit = btPositive?t->t1:t->t2;

  Key key = getKey(trit->index, (Key)s);
  int useTotal = tpiCount[key];

  TrianglePointInteractionCache::iterator tpii = tpiCache.find(key);
  if(tpii != tpiCache.end()) {
    args->tpi = tpii->second;
  } else {
    args->tpi = new TrianglePointInteraction(args->th.v0,args->th.v1,args->th.v2,36,useTotal);
    tpiCache[key] = args->tpi;
  }  

  args->tpi->start();
  QuadratureFunction<ComplexVector3,ComplexVectorValuedFunction>::Type quadfunc;
  quadfunc = &quadgauss36<ComplexVector3,ComplexVectorValuedFunction>;
  ComplexVector3 lk = quadfunc(integrand_LK1, args, args->tpi->v1, args->tpi->v2, NULL, NULL, args->tol, 0.);
  args->tpi->end();
  if(args->tpi->useCount == args->tpi->useTotal) {
    tpiCache.erase(key);
    delete args->tpi;
  }
  return lk;
}

ComplexVector3Pair LK2(int t, RWG *s, FunctionArgs *args, TrianglePointInteractionCache &tpiCache, TrianglePointInteractionCount &tpiCount, bool bsPositive)
{
  getTriangle(s,&(args->sh),1.,bsPositive);
  Triangle *tris = bsPositive?s->t1:s->t2;

  Key key = getKey(tris->index, (Key)t);
  int useTotal = tpiCount[key];

  TrianglePointInteractionCache::iterator tpii = tpiCache.find(key);
  if(tpii != tpiCache.end()) {
    args->tpi = tpii->second;
  } else {
    args->tpi = new TrianglePointInteraction(args->sh.v0,args->sh.v1,args->sh.v2,36,useTotal);
    tpiCache[key] = args->tpi;
  }  

  args->tpi->start();
  QuadratureFunction<ComplexVector3Pair,ComplexVectorPairValuedFunction>::Type quadfunc;
  quadfunc = &quadgauss36<ComplexVector3Pair,ComplexVectorPairValuedFunction>;  
  ComplexVector3Pair lk = quadfunc(integrand_LK2, args, args->tpi->v1, args->tpi->v2, NULL, NULL, args->tol, 0.);
  args->tpi->end();
  if(args->tpi->useCount == args->tpi->useTotal) {
    tpiCache.erase(key);
    delete args->tpi;
  }
  return lk;
}

ComplexVector3 LK1(RWG *t, Point *rs, int sindex, GreensFunction *g, TrianglePointInteractionCache &tpiCache, TrianglePointInteractionCount &tpiCount, Options &options)
{
  FunctionArgs args;
  args.green = g;
  args.tol = options.tol;
  args.rs = rs;

  ComplexVector3 lk(Complex(0.,0.),Complex(0.,0.),Complex(0.,0.));
  lk += LK1(t,sindex,&args,tpiCache,tpiCount,true);
  lk -= LK1(t,sindex,&args,tpiCache,tpiCount,false);
  lk *= t->e->length;

  return lk;
}

ComplexVector3Pair LK2(Point *rt, int tindex, RWG *s, GreensFunction *g, TrianglePointInteractionCache &tpiCache,  TrianglePointInteractionCount &tpiCount, Options &options)
{
  FunctionArgs args;
  args.green = g;
  args.tol = options.tol;
  args.rt = rt;

  ComplexVector3Pair lk2;
  lk2 = LK2(tindex,s,&args,tpiCache,tpiCount,true);
  lk2 -= LK2(tindex,s,&args,tpiCache,tpiCount,false);
  lk2 *= s->e->length;

  return lk2;
}

void countTrianglePointInteractions(RWG *t, int s, bool btPositive, TrianglePointInteractionCount &tpiCount)
{
  Triangle *trit = btPositive?t->t1:t->t2;
  Key key = getKey(trit->index, s);
  TrianglePointInteractionCount::iterator tpii = tpiCount.find(key);
  int count;
  if(tpii != tpiCount.end()) {
    count = tpii->second + 1;
  } else {
    count = 1;
  }
  tpiCount[key] = count;
}

void countTrianglePointInteractions(RWG *t, int s, TrianglePointInteractionCount &tpiCount)
{
  countTrianglePointInteractions(t,s,true,tpiCount);
  countTrianglePointInteractions(t,s,false,tpiCount);
}

// get an array of 3d points distributed over the surface described by the surface parameterization
void getSurfacePoints(int n, SurfaceParameterization &u, Point **point, int *points)
{
  *points = SurfaceParameterization :: hexagonalPoints(n);
  *point = new Point[*points];

  int k = 0;
  for(int p = -n+1; p < n; p++) {
    int q0;
    int q1;
    if(p < 0) q0 = -n - p;
    else q0 = -n;
    if(p > 0) q1 = n - p;
    else q1 = n;
    for(int q = q0+1; q < q1; q++) {
      Point2 p2 = SurfaceParameterization :: hexagonalPoint(n,p,q);
      (*point)[k++] = u.inverse(p2);
    }
  }
}

// get all the edge midpoints
void getSurfacePoints(RWGTree *t, Point **point, int  *points)
{
  *points = (int) t->rwgs;
  *point = new Point[*points];
  
  for(int i=0; i<t->rwgs; i++) {
    RWG *rwg = t->rwg[i];
    Point p = *(rwg->e->p1);
    p += *(rwg->e->p2);
    p *= 0.5;
    (*point)[i] = p;
  }
}

void Zdense(RWGTree *t, 
            RWGTree *s,
            ComplexDenseMatrix *dense,
            GreensFunction *gExt,
            Options &options)
{
  TriangleTriangleInteractionCount ttiCount1;
  TriangleTriangleInteractionCount ttiCount2;
  TriangleTriangleInteractionCache ttiCache1;
  TriangleTriangleInteractionCache ttiCache2;
  
  int m = t->rwgs;
  int n = s->rwgs;
  dense->resize(2*m,2*n);

  bool bSeparate = (t->g != s->g);
  GreensFunction *g1 = gExt;
  GreensFunction *g2 = t->g;
  
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      countTriangleTriangleInteractions(t->rwg[i], s->rwg[j], ttiCount1);
      if(!bSeparate)
        countTriangleTriangleInteractions(t->rwg[i], s->rwg[j], ttiCount2);
    }
  }
  ComplexDenseMatrix &Z = *dense;

  Complex I(0.,1.);
  Complex z1eps = 1./(I * g1->omega * g1->eps);
  Complex z1mu = 1./(I * g1->omega * g1->mu);
  Complex z2eps = 1./(I * g2->omega * g2->eps);
  Complex z2mu = 1./(I * g2->omega * g2->mu);

  for(int i=0; i<m; i++) {
    int ii = 2*i;
    for(int j=0; j<n; j++) {
      int jj = 2*j;

      ComplexPair lk1 = LK(t->rwg[i], s->rwg[j], g1, ttiCache1, ttiCount1, options);
      ComplexPair lk2 = bSeparate ? ComplexPair(0.,0.) : LK(t->rwg[i], s->rwg[j], g2, ttiCache2, ttiCount2, options);

      Z(ii,jj) = (z1eps * lk1.first + z2eps * lk2.first);
      Z(ii,jj+1) = (lk1.second + lk2.second);
      Z(ii+1,jj) = -(lk1.second + lk2.second);
      Z(ii+1,jj+1) = (z1mu * lk1.first + z2mu * lk2.first);
    }
    //fprintf(stderr,"%d/%d\n",i+1,m);
  }
}

void Zrk(RWGTree *t,
         RWGTree *s,
         ComplexRk *rk,
         GreensFunction *gExt,
         SurfaceParameterizationCache *cache,
         omp_lock_t *cacheLock,
         Options &options,
         Real gap)
{
  //printf("r %d x %d\n",t->rwgs,s->rwgs);
  bool bSeparate = (t->g != s->g);

  GreensFunction *g1 = gExt;
  GreensFunction *g2 = t->g;

  // choose values for nt, ns, which determine the number of points in the interpolation
  int nt;
  int ns;

  //Vector3 v = t->centroid - s->centroid;
  //Real distance = sqrt(norm(v));
  //Real maxk = max(g1->normk, g2->normk);
  //Real k2 = Square(max(options.gapNormalization / gap, maxk));
  Real k2 = Square(options.gapNormalization / gap);
  
  Real dt = options.interpolationDensity * pow(t->totalArea * k2, options.interpolationPowerLaw);
  Real ds = options.interpolationDensity * pow(s->totalArea * k2, options.interpolationPowerLaw);

  nt = 3 + (int)dt;
  ns = 3 + (int)ds;

  Point *pt;
  int npt;
  Point *ps;
  int nps;

  if(t->rwgs < options.minParameterizableSegmentSize ||
     (Real)SurfaceParameterization::hexagonalPoints(nt) > options.parameterizationThresh * (Real)t->rwgs) {
    getSurfacePoints(t,&pt,&npt);
  } else {
    SurfaceParameterization *ut;
    if(!t->bParameterizable) {
      cout << "RWGTree not parameterizable\n";
      abort();
    }
    if(cache) {
      omp_set_lock(cacheLock);
      SurfaceParameterizationCache::iterator uti = cache->find(t);
      if(uti != cache->end()) {
        ut = uti->second;
      } else {
        ut = t->parameterize();
        (*cache)[t] = ut;
      }
      omp_unset_lock(cacheLock);
    } else {
      ut = t->parameterize();
    }
    //cout << *ut << "\n";
    getSurfacePoints(nt,*ut,&pt,&npt);
    if(!cache) delete ut;    
  }

  if(s->rwgs < options.minParameterizableSegmentSize ||
     (Real)SurfaceParameterization::hexagonalPoints(ns) > options.parameterizationThresh * (Real)s->rwgs) {
    getSurfacePoints(s,&ps,&nps);
  } else {
    SurfaceParameterization *us;
    if(!s->bParameterizable) {
      cout << "RWGTree not parameterizable\n";
      abort();
    }
    if(cache) {
      omp_set_lock(cacheLock);
      SurfaceParameterizationCache::iterator usi = cache->find(s);
      if(usi != cache->end()) {
        us = usi->second;
      } else {
        us = s->parameterize();
        (*cache)[s] = us;
      }
      omp_unset_lock(cacheLock);
    } else {
      us = s->parameterize();
    }
    //cout << *us << "\n";
    getSurfacePoints(ns,*us,&ps,&nps);
    if(!cache) delete us;
  }

  ComplexDenseMatrix G(2*npt,2*nps);
  for(int i=0; i<npt; i++) {
    int ii = 2 * i;
    for(int j=0; j<nps; j++) {
      int jj = 2 * j;
      G(ii,jj) = g1->scalar(pt+i,ps+j);
      G(ii,jj+1) = Complex(0.,0.);
      G(ii+1,jj) = Complex(0.,0.);
      G(ii+1,jj+1) = bSeparate ? Complex(0.,0.) : g2->scalar(pt+i,ps+j);
    }
  }

  /*
  cout << "t:\n";
  for(int i=0; i<npt; i++)
    cout << pt[i] << "\n";

  cout << "s:\n";
  for(int i=0; i<nps; i++)
    cout << ps[i] << "\n";
  */

  vector<int> r;
  vector<int> c;
  aca(G,r,c,options.acaTolerance);

  int m = t->rwgs;
  int n = s->rwgs;
  int rows = (int) r.size();
  int cols = (int) c.size();

  ComplexVector3Matrix Av00(2*m,cols,0);
  ComplexVector3Matrix Av01(2*m,cols,0);
  ComplexVector3Matrix Av10(2*m,cols,0);
  ComplexVector3Matrix Av11(2*m,cols,0);

  ComplexVector3Matrix Bv00(rows,2*n,0);
  ComplexVector3Matrix Bv01(rows,2*n,0);
  ComplexVector3Matrix Bv10(rows,2*n,0);
  ComplexVector3Matrix Bv11(rows,2*n,0);

  Complex I(0.,1.);
  Complex z1eps = 1./(I * g1->omega * g1->eps);
  Complex z1mu = 1./(I * g1->omega * g1->mu);
  Complex z2eps = 1./(I * g2->omega * g2->eps);
  Complex z2mu = 1./(I * g2->omega * g2->mu);

  TrianglePointInteractionCount tpiCount11;
  TrianglePointInteractionCount tpiCount12;
  TrianglePointInteractionCount tpiCount21;
  TrianglePointInteractionCount tpiCount22;
  for(int i=0; i<m; i++) {
    for(int j=0; j<cols; j++) {
      int jj = c[j];
      int sindex = jj>>1;
      if(jj % 2 == 0) {
        countTrianglePointInteractions(t->rwg[i], sindex, tpiCount11);
      } else {
        countTrianglePointInteractions(t->rwg[i], sindex, tpiCount12);
      }
    }
  }
  for(int i=0; i<n; i++) {
    for(int j=0; j<rows; j++) {
      int jj = r[j];
      int tindex = jj>>1;
      if(jj % 2 == 0) {
        countTrianglePointInteractions(s->rwg[i], tindex, tpiCount21);
      } else {
        countTrianglePointInteractions(s->rwg[i], tindex, tpiCount22);
      }
    }
  }

  TrianglePointInteractionCache tpiCache11;
  TrianglePointInteractionCache tpiCache12;
  TrianglePointInteractionCache tpiCache21;
  TrianglePointInteractionCache tpiCache22;
  for(int i=0; i<m; i++) {
    int ii = i<<1;
    for(int j=0; j<cols; j++) {
      int jj = c[j];
      int sindex = jj>>1;
      if(jj % 2 == 0) {
        ComplexVector3 lk11 = LK1(t->rwg[i], ps+sindex, sindex, g1, tpiCache11, tpiCount11, options);
        if(isnan(lk11)) 
          cerr << "ACHTUNG 1!\n";
        Av00(ii,j) = z1eps * lk11;
        Av01(ii,j) = lk11;
        Av10(ii+1,j) = -lk11;
        Av11(ii+1,j) = z1mu * lk11;
      } else {
        ComplexVector3 lk12 = LK1(t->rwg[i], ps+sindex, sindex, g2, tpiCache12, tpiCount12, options);
        if(isnan(lk12)) 
          cerr << "ACHTUNG 2!\n";
        Av00(ii,j) = z2eps * lk12;
        Av01(ii,j) = lk12;
        Av10(ii+1,j) = -lk12;
        Av11(ii+1,j) = z2mu * lk12;
      }
    }
  }
  for(int i=0; i<n; i++) {
    int ii = i<<1;
    for(int j=0; j<rows; j++) {
      int jj = r[j];
      int tindex = jj>>1;
      if(jj % 2 == 0) {
        ComplexVector3Pair lk21 = LK2(pt+tindex, tindex, s->rwg[i], g1, tpiCache21, tpiCount21, options);
        if(isnan(lk21.first) || isnan(lk21.second)) 
          cerr << "ACHTUNG 3!\n";
        Bv00(j,ii) = lk21.first;
        Bv01(j,ii+1) = lk21.second;
        Bv10(j,ii) = lk21.second;
        Bv11(j,ii+1) = lk21.first;
      } else {
        ComplexVector3Pair lk22 = LK2(pt+tindex, tindex, s->rwg[i], g2, tpiCache22, tpiCount22, options);
        if(isnan(lk22.first) || isnan(lk22.second)) 
          cerr << "ACHTUNG 3!\n";
        Bv00(j,ii) = lk22.first;
        Bv01(j,ii+1) = lk22.second;
        Bv10(j,ii) = lk22.second;
        Bv11(j,ii+1) = lk22.first;
      }
    }
  }

  delete [] pt;
  delete [] ps;

  ComplexDenseMatrix S = G(r,c);

  ComplexDenseMatrix A00(2*m,3*cols);
  ComplexDenseMatrix A01(2*m,3*cols);
  ComplexDenseMatrix A10(2*m,3*cols);
  ComplexDenseMatrix A11(2*m,3*cols);

  for(int i=0; i<m; i++) {
    int ii = 2 * i;
    int ii1 = ii + 1;
    for(int j=0; j<cols; j++) {
      int jj = 3 * j;
      ComplexVector3 v;

      v = Av00(ii,j);
      A00(ii,jj) = v.x;
      A00(ii,jj+1) = v.y;
      A00(ii,jj+2) = v.z;
      v = Av00(ii1,j);
      A00(ii1,jj) = v.x;
      A00(ii1,jj+1) = v.y;
      A00(ii1,jj+2) = v.z;

      v = Av01(ii,j);
      A01(ii,jj) = v.x;
      A01(ii,jj+1) = v.y;
      A01(ii,jj+2) = v.z;
      v = Av01(ii1,j);
      A01(ii1,jj) = v.x;
      A01(ii1,jj+1) = v.y;
      A01(ii1,jj+2) = v.z;

      v = Av10(ii,j);
      A10(ii,jj) = v.x;
      A10(ii,jj+1) = v.y;
      A10(ii,jj+2) = v.z;
      v = Av10(ii1,j);
      A10(ii1,jj) = v.x;
      A10(ii1,jj+1) = v.y;
      A10(ii1,jj+2) = v.z;

      v = Av11(ii,j);
      A11(ii,jj) = v.x;
      A11(ii,jj+1) = v.y;
      A11(ii,jj+2) = v.z;
      v = Av11(ii1,j);
      A11(ii1,jj) = v.x;
      A11(ii1,jj+1) = v.y;
      A11(ii1,jj+2) = v.z;
    }
  }

  S.LU();
  solveLU(S,Bv00);
  solveLU(S,Bv01);
  solveLU(S,Bv10);
  solveLU(S,Bv11);

  ComplexDenseMatrix B00(3*rows,2*n);
  ComplexDenseMatrix B01(3*rows,2*n);
  ComplexDenseMatrix B10(3*rows,2*n);
  ComplexDenseMatrix B11(3*rows,2*n);
  for(int i=0; i<rows; i++) {
    int ii = 3 * i;
    for(int j=0; j<n; j++) {
      int jj = 2 * j;
      int jj1 = jj + 1;
      ComplexVector3 v;

      v = Bv00(i,jj);
      B00(ii,jj) = v.x;
      B00(ii+1,jj) = v.y;
      B00(ii+2,jj) = v.z;
      v = Bv00(i,jj1);
      B00(ii,jj1) = v.x;
      B00(ii+1,jj1) = v.y;
      B00(ii+2,jj1) = v.z;

      v = Bv01(i,jj);
      B01(ii,jj) = v.x;
      B01(ii+1,jj) = v.y;
      B01(ii+2,jj) = v.z;
      v = Bv01(i,jj1);
      B01(ii,jj1) = v.x;
      B01(ii+1,jj1) = v.y;
      B01(ii+2,jj1) = v.z;

      v = Bv10(i,jj);
      B10(ii,jj) = v.x;
      B10(ii+1,jj) = v.y;
      B10(ii+2,jj) = v.z;
      v = Bv10(i,jj1);
      B10(ii,jj1) = v.x;
      B10(ii+1,jj1) = v.y;
      B10(ii+2,jj1) = v.z;

      v = Bv11(i,jj);
      B11(ii,jj) = v.x;
      B11(ii+1,jj) = v.y;
      B11(ii+2,jj) = v.z;
      v = Bv11(i,jj1);
      B11(ii,jj1) = v.x;
      B11(ii+1,jj1) = v.y;
      B11(ii+2,jj1) = v.z;
    }
  }         

  ComplexRk rk00(A00,B00,options.rkEps);
  ComplexRk rk01(A01,B01,options.rkEps);
  ComplexRk rk10(A10,B10,options.rkEps);
  ComplexRk rk11(A11,B11,options.rkEps);

  ComplexRk rk0 = rk00 + rk01;
  ComplexRk rk1 = rk10 + rk11;

  *rk = rk0 + rk1;
  //fprintf(stderr,"t: %d/%d/%d/%d, s: %d/%d/%d/%d\n",nt,rk->k,rows,t->rwgs,ns,rk->k,cols,s->rwgs);
  //fprintf(stderr,"rk compresion: %g\n", (Real)((rk->m+rk->n)*rk->k)/(Real)(rk->m*rk->n));
}

Real projectedDistance(RWGTree *t, RWGTree *s, Vector3 &vst)
{
  Real maxd = 0.;
  for(int i=0;i<s->rwgs;i++) {
    RWG *rwg = s->rwg[i];

    Triangle *t;
    Vector3 v;
    Real d;

    t = rwg->t1;

    v = *(t->p1) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p2) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p3) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    t = rwg->t2;

    v = *(t->p1) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p2) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;

    v = *(t->p3) - s->centroid;
    d = dot(v,vst);
    if(d > maxd)
      maxd = d;
  }
  return maxd;
}

Real admissable(RWGTree *t, RWGTree *s, Options &options)
{
  if(t==s) return false;
  if(t->rwgs < options.minAdmissableSegmentSize || s->rwgs < options.minAdmissableSegmentSize) return false;
  if(!t->bParameterizable || !s->bParameterizable) return false;
  Vector3 v = t->centroid - s->centroid;
  Real d = sqrt(norm(v));
  v *= 1. / d;
  Real d1 = projectedDistance(t,s,v);
  v *= -1.;
  Real d2 = projectedDistance(s,t,v);
  //printf("%d %d %g %g %g\n",t->rwgs,s->rwgs,d,d1,d2);
  Real gap = d - (d1 + d2);
  Real r = max(t->radius,s->radius);
  if(gap > r * options.admissableThresh)
    return gap;
  else 
    return -1.;
}

ComplexMatrix * Zh(RWGTree *t,
                   RWGTree *s,
                   GreensFunction *gExt,
                   ImpedanceJobVector &jobs,
                   Options &options)
{
  Real gap = admissable(t,s,options);
  if(gap > 0.) {
    ComplexRk *rk = new ComplexRk;
    ImpedanceJob job(t,s,rk,gap);
    jobs.push_back(job);
    return rk;
  }
  if(!t->sub && !s->sub) {
    ComplexDenseMatrix *dense = new ComplexDenseMatrix;
    ImpedanceJob job(t,s,dense,gap);
    jobs.push_back(job);
    return dense;
  }
  if(!t->sub) {
    ComplexHMatrix *Hp = new ComplexHMatrix(1,s->subs);
    ComplexHMatrix &H = *Hp;
    H.r[0] = 2*t->rwgs;
    for(int j=0; j<s->subs; j++)
      H.c[j] = 2*s->sub[j]->rwgs;
    for(int j=0; j<s->subs; j++) {
      ComplexMatrix *M = Zh(t,s->sub[j],gExt,jobs,options);
      H.set(0,j,M);
    }
    return Hp;
  }
  if(!s->sub) {
    ComplexHMatrix *Hp = new ComplexHMatrix(t->subs,1);
    ComplexHMatrix &H = *Hp;
    for(int i=0; i<t->subs; i++)
      H.r[i] = 2*t->sub[i]->rwgs;
    H.c[0] = 2*s->rwgs;
    for(int i=0; i<t->subs; i++) {
      ComplexMatrix *M = Zh(t->sub[i],s,gExt,jobs,options);
      H.set(i,0,M);
    }
    return Hp;
  }

  ComplexHMatrix *Hp = new ComplexHMatrix(t->subs,s->subs);
  ComplexHMatrix &H = *Hp;
  for(int i=0; i<t->subs; i++)
    H.r[i] = 2*t->sub[i]->rwgs;
  for(int j=0; j<s->subs; j++)
    H.c[j] = 2*s->sub[j]->rwgs;
  for(int i=0; i<t->subs; i++) {
    for(int j=0; j<s->subs; j++) {
      ComplexMatrix *M = Zh(t->sub[i],s->sub[j],gExt,jobs,options);
      H.set(i,j,M);
    }
  }
  return Hp;
}


ComplexMatrix * Z(RWGTree *t,
                  RWGTree *s,
                  GreensFunction *gExt,
                  Progress &progress,
                  Options &options)
{
  ImpedanceJobVector jobs;
  ComplexMatrix *Zhp = Zh(t,s,gExt,jobs,options);

  SurfaceParameterizationCache cache;
  omp_lock_t cacheLock;
  omp_lock_t progressLock;
  omp_init_lock(&cacheLock);
  omp_init_lock(&progressLock);

  int chunk = 1;
  int n = jobs.size();
  int i;
#pragma omp parallel for default(shared) private(i) schedule(static,chunk) 
  for(i=0; i<n; i++) {
    ImpedanceJob &job = jobs[i];
    job.perform(gExt,&cache,&cacheLock,progress,&progressLock,options);
  }

  return Zhp;
}

ImpedanceJob :: ImpedanceJob(RWGTree *t, RWGTree *s, ComplexMatrix *M, Real gap)
{
  this->t = t;
  this->s = s;
  this->M = M;
  this->gap = gap;
}

void ImpedanceJob :: perform(GreensFunction *gExt,
                             SurfaceParameterizationCache *cache,
                             omp_lock_t *cacheLock,
                             Progress &progress,
                             omp_lock_t *progressLock,
                             Options &options)
{
  if(M->type == ComplexMatrix::dense) {
    ComplexDenseMatrix *dense = (ComplexDenseMatrix*)M;
    Zdense(t,s,dense,gExt,options);
    omp_set_lock(progressLock);
    progress.finishedMatrix(*dense);
    omp_unset_lock(progressLock);
  } else {
    ComplexRk *rk = (ComplexRk*)M;
    Zrk(t,s,rk,gExt,cache,cacheLock,options,gap);
    omp_set_lock(progressLock);
    progress.finishedMatrix(*rk);
    omp_unset_lock(progressLock);
  }
}

void LKdense(RWGTree *t, 
             RWGTree *s,
             ComplexDenseMatrix *L_,
             ComplexDenseMatrix *K_,
             GreensFunction *gExt,
             Options &options)
{
  TriangleTriangleInteractionCount ttiCount1;
  TriangleTriangleInteractionCache ttiCache1;
  
  int m = t->rwgs;
  int n = s->rwgs;

  GreensFunction *g1 = s->g;
    
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      countTriangleTriangleInteractions(t->rwg[i], s->rwg[j], ttiCount1);
    }
  }

  ComplexDenseMatrix &L = *L_;
  ComplexDenseMatrix &K = *K_;
  L.resize(m,n);
  K.resize(m,n);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      ComplexPair lk = LK(t->rwg[i], s->rwg[j], g1, ttiCache1, ttiCount1, options);
      L(i,j) = lk.first;
      K(i,j) = lk.second;
    }
  }
}

#undef DEBUG
