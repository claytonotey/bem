#ifndef SILICONDOPED_H
#define SILICONDOPED_H

#include "Dielectric.h"

class SiliconDoped : public Dielectric {
 public:
  SiliconDoped(Real N);
  Complex eps(Real eV);
  Complex mu(Real eV);
  Real N;
};

#endif
