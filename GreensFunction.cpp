#include "GreensFunction.h"
#include "MathUtils.h"

#define GREENEPSR 1e-9

Real GreensFunction :: scale;

GreensFunction :: GreensFunction(Dielectric *diel, bool bExternal, bool bNonlocal)
{
  this->bExternal = bExternal;
  this->bNonlocal = bNonlocal;
  this->diel = diel;
  if(!bNonlocal && diel->isNonlocal()) {
    gLong = new GreensFunction(diel, bExternal, true);
  } else {
    gLong = NULL;
  }
}

void GreensFunction :: setLengthScale(Real s) 
{
  GreensFunction :: scale = s;
}

void GreensFunction :: init(Real feV)
{
  if(gLong) gLong->init(feV);
  this->eps = diel->eps(feV);
  this->mu = diel->mu(feV);
  this->eps0 = diel->eps0(feV);
  // unitless
  this->omega = feV * SI::eV / SI::hbar / SI :: c * scale;
  ksquared = diel->dispersion_k2(feV,bNonlocal) * scale * scale;
  k = sqrt(ksquared);
  Complex I(0.,1.);
  ik = I * k;
  k2 = -0.5 * ksquared;
  normk = sqrt(norm(k));
  if(bNonlocal) {
    Zeps =  -I / omega * diel->Zeps(feV);
    Zmu = 0.0f;
  } else {
    Zeps = mu * omega / (I * ksquared);
    Zmu = eps * omega / (I * ksquared);
  }
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

Complex GreensFunction :: scalar_smooth(Vector3 *rt, Vector3 *rs)
{
  Real x = rt->x - rs->x; 
  Real y = rt->y - rs->y; 
  Real z = rt->z - rs->z; 
  Real R2 = x*x + y*y + z*z;
  Real R = sqrt(R2);
  Complex ikR = ik*R;
  Complex expikR = exp(ikR);
  return (expikR - 1.0) / R - k2 * R;
}

bool GreensFunction :: isExternal()
{
  return bExternal;
}
