#ifndef FREQUENCYSWEEP_H
#define FREQUENCYSWEEP_H

#include "defs.h"

using namespace std;
#include <vector>

class FrequencySweep {
 public:
  FrequencySweep(char *file);
  vector<Real> f;
};

#endif
