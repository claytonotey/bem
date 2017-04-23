#include "SiC.h"
#include "SI.h"

SiC :: SiC()
{
  setTemperature(300);
}

void SiC :: setTemperature(Real T) {
  epsinf = 6.7 * exp(5e-5 * (T-300));
  wL = 182.7e12 - 5.463e9 * (T-300);
  wT = 149.5e12 - 4.106e9 * (T-300);
  gamma = 6.6e11 * (1 + 2.0 / (exp(SI::hbar * wT / (2*SI::kb*T)) - 1));
}

Complex SiC :: eps(Real feV) {
  Real w = feV * SI::eV / SI::hbar;
  return epsinf * (wL * wL - w * w - Complex(0,gamma*w)) / (wT * wT - w * w - Complex(0,gamma*w));
}

Complex SiC :: mu(Real feV) {
  return Complex(1.,0.);
}
