#ifndef DIELECTRIC_H
#define DIELECTRIC_H

#include "SI.h"
#include "defs.h"
#include "MathUtils.h"

class Dielectric {
 public:
  Dielectric() {}
  virtual ~Dielectric() {};
  virtual Complex eps(Real eV)=0;
  virtual Complex mu(Real eV)=0;
  virtual bool isNonlocal() { return false; }
  
  Complex lorentzianPole(Real Delta_epsilon, Real omega1, Real delta1, Real omega) {
    Real omega12 = Square(omega1);
    return Delta_epsilon * omega12 / (omega12 - omega * (omega + Complex(0,2*delta1)));
  }

  // k in SI m^-1
  virtual Complex dispersion_k2(Real feV, bool bLongitudinal = false) {
    Real omega = feV * SI::eV / SI::hbar;
    return eps(feV)*mu(feV)*Square(omega/ SI::c);
  }

  virtual Complex eps0(Real feV) {
    return eps(feV);
  }
  virtual Complex Zeps(Real feV) {
    abort();
  }

  // old school form 
  Complex lorentzianPole(Real rho, Real feV0, Real gamma, Real eps0, Real feV) {

    Real f = feV/feV0;
    Real f2 = f * f;
    Real gamma2 = gamma * gamma;
    Real d = Square(1. - f2) + gamma2 * f2;
    Real xi = rho * (1. - f2) / d;
    Real sigmaoverf = TWOPI * rho * f * gamma / d;
    Real eps = eps0 + FOURPI * xi;
    d = sqrt(Square(eps) + 4. * Square(sigmaoverf));
    Real n = sqrt(0.5 * (d + eps));
    Real k = sqrt(0.5 * (d - eps));

    return Square(Complex(n,k));
  }
  
};

class ConstantDielectric : public Dielectric {
 public:
  ConstantDielectric(Complex eps, Complex mu) : _eps(eps), _mu(mu) {}
  virtual ~ConstantDielectric() {};
  Complex eps(Real eV) { return _eps; }
  Complex mu(Real eV) { return _mu; }
  Complex _eps;
  Complex _mu;
};

class Vacuum : public Dielectric {
 public:
  Complex eps(Real eV) { return Complex(1.0,0.0); }
  Complex mu(Real eV) { return Complex(1.0,0.0); }
};

#endif
