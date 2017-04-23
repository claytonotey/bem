#include "Dielectrics.h"
#include "HeatTransfer.h"
#include "Util.h"
#include "Params.h"
#include <string>
#include <functional>
using namespace std;

//typedef std::unordered_map<pair<string, Params>, string, paramsHasher > RegistryMap;

int main(int argc, char **argv)
{
  

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
  
  //RegistryMap hi;
  //hi[pair<string,Params>("hi",Params())] = string("hi");
  
  //cout << hi[Params()] << "\n";
  Params params;
  params["N"] = ParamFloat(4.0);
  Dielectric *diel = Dielectrics::getDielectric("SiliconDoped",params);
  cout << diel->eps(eV) << "\n";
  //parseShapes(top, s, t, body1x3d, body2x3d, gap, options);
}
