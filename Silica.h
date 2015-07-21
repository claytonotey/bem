#ifndef SILICA_H
#define SILICA_H

#include "Dielectric.h"

#define SILICA_TABLE_SIZE 85

class Silica : public Dielectric {
 public:
  Complex eps(Real eV);
  Complex mu(Real eV);
  static Real ftab[SILICA_TABLE_SIZE];
  static Real ntab[SILICA_TABLE_SIZE];
  static Real ktab[SILICA_TABLE_SIZE];
};

#endif
