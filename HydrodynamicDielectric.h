#ifndef HYDRODYNAMIC_DIELECTRIC_H
#define HYDRODYNAMIC_DIELECTRIC_H

#include "Dielectric.h"

class Au : public Dielectric {
 public:

  static const Real epsilon_inf;
  static const Real omega_D;
  static const Real gamma;
  static const Real deltaEpsilon1;
  static const Real omega1;
  static const Real delta1;
  static const Real deltaEpsilon2;
  static const Real omega2;
  static const Real delta2;
  static const Real vFermi;

  bool isNonlocal();
  Complex dispersion_k2(Real feV, bool bNonLocal);
  Complex eps(Real eV);
  Complex mu(Real eV);
  Complex Zeps(Real eV);
  Complex eps0(Real eV);
};

#endif
