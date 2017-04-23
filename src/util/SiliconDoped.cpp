#include "SiliconDoped.h"
#include "SI.h"

SiliconDoped :: SiliconDoped(Real N) : N(N) {}

Complex SiliconDoped :: eps(Real eV)
{
  Real omega = eV * SI::eV / SI::hbar;
  if(N < 0.0) {
    Real mu1 = .00685;
    Real mumax = .1414;
    Real mu2 = .00561;
    Real Cr = 9.20e22;
    Real Cs = 3.41e26;
    Real alpha = .711;
    Real beta = 1.98;
    Real m = .27 * SI::me;
    Real mu = mu1 + (mumax - mu1) / (1 + pow(fabs(N) / Cr, alpha)) - mu2 / (1 + pow(Cs / fabs(N), beta));
    Real omegap = SI::eV * sqrt(fabs(N) / m / SI::eps0);
    Real gamma = SI::eV / m / mu;
    return (Real)11.7 - omegap * omegap / (omega * (omega + Complex(0.,gamma)));
  } else if(N > 0.0) {
    // ptype
    Real mu1 = .00449;
    Real mumax = .04705;
    Real mu2 = .0029;
    Real Cr = 2.23e23;
    Real Cs = 6.10e26;
    Real alpha = .719;
    Real beta = 2.00;
    Real pc = 9.23e24;
    Real m = .37 * SI::me;
    Real mu = mu1 * exp(-pc / N) + mumax / (1 + pow(N / Cr, alpha)) - mu2 / (1 + pow(Cs / N, beta));
    Real omegap = SI::eV * sqrt(N / m / SI::eps0);
    Real gamma = SI::eV / m / mu;
    return (Real)11.7 - omegap * omegap / (omega * (omega + Complex(0.,gamma)));
  } else {
     return 11.7;
  }
}

Complex SiliconDoped :: mu(Real eV)
{
  return Complex(1.,0.);
}

