#include "SiliconCarbide.h"
#include "SI.h"

const Real SiliconCarbide :: rho = 0.263;
const Real SiliconCarbide :: feV0 = 0.098356643807595;
const Real SiliconCarbide :: gamma = 0.007;
const Real SiliconCarbide :: eps0 = 6.7;

Complex SiliconCarbide :: eps(Real feV) {
  return lorentzianPole(rho,feV0,gamma,eps0,feV);
}

Complex SiliconCarbide :: mu(Real feV) {
  return Complex(1.,0.);
}
