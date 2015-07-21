#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include "Vector3.h"
#include "ComplexVector3.h"
#include "Dielectric.h"

class GreensFunction {
 public:
  static Real scale;
  static void setLengthScale(Real s);

  GreensFunction(Dielectric *diel, bool bExternal, bool bNonLocal = false);
  void init(Real feV);
  void scalar_gradient(Vector3 *rt, Vector3 *rs, Complex *g, ComplexVector3 *gg);
  void scalar_gradient_smooth(Vector3 *rt, Vector3 *rs, Complex *g, ComplexVector3 *gg);
  Complex scalar_smooth(Vector3 *rt, Vector3 *rs);
  Complex scalar(Vector3 *rt, Vector3 *rs);
  bool isExternal();

  Dielectric *diel;
  GreensFunction *gLong;
  bool bExternal;
  bool bNonlocal;
  Complex eps0;
  Complex Zeps;
  Complex Zmu;
  Real normk;
  Complex ik;
  Complex k2;
  Complex k;
  Complex ksquared;
  Complex eps;
  Complex mu;
  Complex omega;
};

#endif
