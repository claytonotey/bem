#include "Al.h"
#include "SI.h"

const Real Al :: epsInf = 2.2281;
const Real Al :: epsLorentz = 1.30601;
const Real Al :: omega0 = 1.50690679e14;
const Real Al :: gamma0 = 3.76734192e13;

Complex Al :: eps(Real feV) {
  Real omega = feV * SI::eV / SI::hbar;
  return epsInf +  epsLorentz * omega0 * omega0 / (omega0 * omega0 - omega * omega - Complex(0,2.0*omega*gamma0));
}

Complex Al :: mu(Real feV) {
  return Complex(1.,0.);
}
