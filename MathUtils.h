#ifndef MATHUTILS_H
#define MATHUTILS_H

#include "defs.h"

#define THREERHALVESSQRT3 2.5980762113533159402911695122588
#define SQRT3OVER2 0.86602540378443864676372317075294
#define LOG2 0.69314718055994528622676398299518
#define LOG4 1.3862943611198905724535279659904
#define ONESIXTH 0.16666666666666666666666666666667
#define ONESIXTEENTH 0.0625
#define ONETWENTYFOURTH 0.041666666666666666666666666666667
#define ONETHIRD 0.33333333333333333333333333333333
#define TWOTHIRDS 0.6666666666666666666666666666667
#define ONEOVER48 0.0208333333333333333333333333333
#define FIVESIXTHS 0.833333333333333333333333333333
#define CUBEROOT2 1.2599210498948731906665443602833
#define CUBEROOT2SQUARED 1.5874010519681993613971826562192
#define SQRT3 1.7320508075688772935274463415059
#define SQRT6 2.4494897427831780981972840747059
#define ONEOVERSQRT3 0.57735026918962576450914878050196
#define ONEOVERSQRT6OVER2 0.2041241452319315363705953814133
#define SQRT2 1.4142135623730950488016887242097
#define TWOPI 6.283185307179586476925286766559
#define PI 3.1415926535897932384626433832795
#define PIOVER2 1.5707963267948966192313216916398
#define FOURPI 12.566370614359172953850573533118
#define TWOOVERPI 0.63661977236758138243288840385503
#define ONEOVER32PISQUARED 0.0031662869888230554815677919577865
#define ONEOVER4PI 0.079577471545947667884441881686257
#define ONEOVER16PI 0.019894367886486916971110470421564

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

struct PointAndVector {
  Real *p;
  Real *xi;
};

void mnbrak(Real *ax,
            Real *bx,
            Real *cx,
            Real *fa,
            Real *fb,
            Real *fc,
            Real (*func)(Real, void *args),
            void *args);
void powell(Real *p, 
            Real **xi,
            int n,
            Real ftol,
            Real lintol,
            int *iter,
            Real *fret,
            Real (*func)(Real, void *),
            void *args,
            PointAndVector *pv);
Real brent(Real ax, 
           Real bx,
           Real cx,
           Real (*f)(Real, void *),
           void *args,
           Real tol,
           Real *xmin);
Real dbrent(Real ax,
            Real bx,
            Real cx,
            Real (*f)(Real, void*),
            Real (*df)(Real, void *),
            void *args,
            Real tol,
            Real *xmin);

int searchTable(Real *tab, int n, Real x);
Real cubicHermiteSplineInterpolate(Real *xtab, Real *ytab, int n, int k, Real x);
Real cubicHermiteSplineInterpolate(Real *xtab, Real *ytab, int n, Real x);

inline Real Square(Real x)
{
  return x*x;
}

inline LongReal Square(LongReal x)
{
  return x*x;
}

inline Complex Square(Complex x)
{
  return x*x;
}

inline Real norm(Real x)
{
  return x*x; 
}

Complex parseComplex(char *str);

static const char* errmsg[1] = {"Error Parsing Complex"};

#define MATHUTILS_PARSECOMPLEX_EXCEPTION MathUtilsException(0)

class MathUtilsException: public exception
{
 public:
  MathUtilsException(int type)
  {
    this->type = type;
  }
  virtual const char* what() const throw()
  {
    return errmsg[type];
  }

 protected:
  int type;
  
};


#endif
