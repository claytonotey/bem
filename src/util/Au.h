#ifndef AU_H
#define AU_H

#include "Dielectric.h"

class Au : public Dielectric {
 public:

  static const Real omegap;
  static const Real omegac;
  static const Real eps0;

  Complex eps(Real eV);
  Complex mu(Real eV);
};

#endif
