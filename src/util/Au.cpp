#include "Au.h"
#include "SI.h"

const Real Au :: eps0 = 1.0;
const Real Au :: omegap = 1.366592804e16;
const Real Au :: omegac = 4.08407045e13;

Complex Au :: eps(Real feV) {
  Real omega = feV * SI::eV / SI::hbar;
  return eps0 - omegap * omegap / (omega * (omega + Complex(0,omegac)));
}

Complex Au :: mu(Real feV) {
  return Complex(1.,0.);
}
