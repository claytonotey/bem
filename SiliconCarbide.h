#ifndef SILCONCARBIDE_H
#define SILCONCARBIDE_H

#include "Dielectric.h"

class SiliconCarbide : public Dielectric {
 public:

  static const Real rho;
  static const Real feV0;
  static const Real gamma;
  static const Real eps0;

  Complex eps(Real eV);
  Complex mu(Real eV);
};

#endif
