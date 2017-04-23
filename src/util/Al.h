#ifndef AL_H
#define AL_H

#include "Dielectric.h"

class Al : public Dielectric {
 public:

  static const Real omega0;
  static const Real gamma0;
  static const Real epsInf;
  static const Real epsLorentz;

  Complex eps(Real eV);
  Complex mu(Real eV);
};

#endif
