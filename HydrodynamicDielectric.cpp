#include "HydrodynamicDielectric.h"
#include <iostream>

const Real Au :: epsilon_inf = 3.559;
const Real Au :: omega_D = 8.812; //eV
const Real Au :: gamma = 0.0752;
const Real Au :: deltaEpsilon1 = 2.912;
const Real Au :: omega1 = 4.693;
const Real Au :: delta1 = 1.541;
const Real Au :: deltaEpsilon2 = 1.272;
const Real Au :: omega2 = 3.112;
const Real Au :: delta2 = 0.525;
const Real Au :: vFermi = 1.39e6;


bool Au :: isNonlocal() 
{
  return true;
}

Complex Au :: dispersion_k2(Real feV, bool bLongitudinal)
{
  if(bLongitudinal) {
    Real beta2 = 3.0 * Square(vFermi) / 5.0;
    Real omega = feV;
    Complex eps_int = lorentzianPole(deltaEpsilon1,omega1,delta1,omega) + lorentzianPole(deltaEpsilon2,omega2,delta2,omega);

    return (-Square(omega_D) / (epsilon_inf + eps_int) + (omega * (omega + Complex(0,gamma)))) / beta2 * Square(SI::eV / SI::hbar);
  } else {
    return Dielectric :: dispersion_k2(feV);
  }
}

Complex Au :: eps0(Real feV)
{
  Real omega = feV;
  Complex eps0 = epsilon_inf + lorentzianPole(deltaEpsilon1,omega1,delta1,omega) + lorentzianPole(deltaEpsilon2,omega2,delta2,omega);
  return eps0;
}

Complex Au :: Zeps(Real feV) 
{
  Real omega = feV;
  Complex eps0 = epsilon_inf + lorentzianPole(deltaEpsilon1,omega1,delta1,omega) + lorentzianPole(deltaEpsilon2,omega2,delta2,omega);

  return Square(omega_D) / (eps0 * (-eps0 * omega * (omega + Complex(0,gamma)) + Square(omega_D)));
}

Complex Au :: eps(Real feV) 
{
  Real omega = feV;
  Complex eps_int = lorentzianPole(deltaEpsilon1,omega1,delta1,omega) + lorentzianPole(deltaEpsilon2,omega2,delta2,omega);
  Complex eps_intra = -Square(omega_D) / (omega * (omega + Complex(0,gamma)));
  return epsilon_inf + eps_int + eps_intra;
}

Complex Au :: mu(Real feV) {
  return Complex(1.,0.); 
}



