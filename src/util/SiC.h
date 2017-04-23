#ifndef SILCONCARBIDE_H
#define SILCONCARBIDE_H

#include "Dielectric.h"

class SiC : public Dielectric {
 public:
  SiC();
  Real epsinf, gamma, wL, wT;
  void setTemperature(Real T);
  Complex eps(Real eV);
  Complex mu(Real eV);
};

#endif
