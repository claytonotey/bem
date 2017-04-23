#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include "Vector3.h"
#include "ComplexVector3.h"

class GreensFunction {
 public:
  GreensFunction(bool bExternal);
  GreensFunction(Complex eps, Complex mu, Complex omega, bool bExternal = false);
  void init(Complex eps, Complex mu, Complex omega);
  void scalar_gradient(Vector3 *rt, Vector3 *rs, Complex *g, ComplexVector3 *gg);
  void scalar_gradient_smooth(Vector3 *rt, Vector3 *rs, Complex *g, ComplexVector3 *gg);
  Complex scalar(Vector3 *rt, Vector3 *rs);
  bool isExternal();

  bool bExternal;
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
