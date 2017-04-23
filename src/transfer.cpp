#include "HeatTransfer.h"

int main(int argc, char **argv)
{
  initReferenceCountingPointerLock();

  if(argc != 10) {
    cerr << "Usage: transfer <x3d file for body 1> <x3d file for body 2> <scale factor applied to x3d files (in meters)> <T_1 (in K)> <T_2 (in K)> <gap x y z>" << "\n";
    exit(-2);
  }
  char *body1x3d = argv[1];
  char *body2x3d = argv[2];
  Real eV = atof(argv[3]);
  Real a = atof(argv[4]);
  Real T1 = atof(argv[5]);
  Real T2 = atof(argv[6]);
  Real gapx = atof(argv[7]);
  Real gapy = atof(argv[8]);
  Real gapz = atof(argv[9]);

  Vector3 gap(gapx, gapy, gapz);
  Options options;
  RWGTree *top;
  RWGTree *s;
  RWGTree *t;

  parseShapes(top, s, t, body1x3d, body2x3d, gap, options);
  CalcHeatTransfer(a, eV, T1, T2, top, s, t, options);
}
