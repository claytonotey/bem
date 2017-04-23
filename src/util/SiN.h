#ifndef SiN_H
#define SiN_H
#include "Dielectric.h"
#define SiN_TABLE_SIZE 1114
class SiN : public Dielectric {
public:
Complex eps(Real eV);
Complex mu(Real eV);
static Real ftab[SiN_TABLE_SIZE];
static Real ntab[SiN_TABLE_SIZE];
static Real ktab[SiN_TABLE_SIZE];
};
#endif