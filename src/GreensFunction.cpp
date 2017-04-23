#include "GreensFunction.h"
#include "MathUtils.h"

#define GREENEPSR 1e-9

GreensFunction :: GreensFunction(bool bExternal)
{
  this->bExternal = bExternal;
}

GreensFunction :: GreensFunction(Complex eps, Complex mu, Complex omega, bool bExternal)
{
  this->bExternal = bExternal;
  init(eps,mu,omega);
}

void GreensFunction :: init(Complex eps, Complex mu, Complex omega)
{
  this->eps = eps;
  this->mu = mu;
  this->omega = omega;
  ksquared = eps*mu*omega*omega;
  k = sqrt(ksquared);
  Complex I(0.,1.);
  ik = I * k;
  k2 = -0.5 * ksquared;
  normk = sqrt(norm(k));
  Zeps = 1. / (I * omega * eps);
  Zmu = 1. / (I * omega * mu);
}

Complex GreensFunction :: scalar(Vector3 *rt, Vector3 *rs)
{
  Real x = rt->x - rs->x; 
  Real y = rt->y - rs->y; 
  Real z = rt->z - rs->z; 
  Real R = sqrt( x*x + y*y + z*z );
  return exp(ik*R) / R;
}

void GreensFunction :: scalar_gradient(Vector3 *rt, Vector3 *rs, Complex *g, ComplexVector3 *gg)
{
  Real x = rt->x - rs->x; 
  Real y = rt->y - rs->y; 
  Real z = rt->z - rs->z; 
  Real R2 = x*x + y*y + z*z;
  Real R = sqrt(R2);
  Complex ikR = ik*R;
  *g = exp(ikR) / R;
  Complex q = *g * (1. - ikR) / R2; 
  *gg = ComplexVector3(q*x,q*y,q*z);
}

void GreensFunction :: scalar_gradient_smooth(Vector3 *rt, Vector3 *rs, Complex *g, ComplexVector3 *gg)
{
  Real x = rt->x - rs->x; 
  Real y = rt->y - rs->y; 
  Real z = rt->z - rs->z; 
  Real R2 = x*x + y*y + z*z;
  Real R = sqrt(R2);
  Complex ikR = ik*R;
  Complex expikR = exp(ikR);
  *g = (expikR - 1.0) / R - k2 * R;

  Complex q;
  Real R3 = R2*R;
  q = (expikR * (1.0 - ikR) + k2 * R2 - 1.0) / R3;
  *gg = ComplexVector3(q*x,q*y,q*z);
}

bool GreensFunction :: isExternal()
{
  return bExternal;
}
