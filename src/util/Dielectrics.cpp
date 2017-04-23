#include "Dielectrics.h"
#include "MathUtils.h"
#include "Silica.h"
#include "SiC.h"
#include "SiliconDoped.h"
#include "SiN.h"
#include "Al.cpp"
#include "Au.cpp"

#include <iostream>
using namespace std;

DielectricRegistry Dielectrics :: registry;

Dielectric *Dielectrics :: getDielectric(const std::string &s)
{
  const DielectricKey key = DielectricKey(s,Params());
  return registry.the[key];
}

Dielectric *Dielectrics :: getDielectric(const std::string &s, const Params &params)
{
  Dielectric *diel;
  const DielectricKey key = DielectricKey(s,params);
  diel = registry.the[key];

  if(!diel) {
    if(s == "SiliconDoped") {
      Real N = 0;
      ParamValue Np = params["N"];
      if(Np.type == ParamValue::ParamFloat) {
        N = Np.f;
      }
      diel = new SiliconDoped(N);
    }
  }

  if(!diel) {
    Complex eps(1.0,0);
    Complex mu(1.0,0);
    ParamValue epsP = params["eps"];
    if(epsP.type == ParamValue::ParamComplex) {
      eps = epsP.z;
    } 
    ParamValue muP = params["mu"];
    if(muP.type == ParamValue::ParamComplex) {
      mu = muP.z;
    }
    diel = new ConstantDielectric(eps,mu);
  }

  registry.the[key] = diel;
  
  return diel;
}

DielectricRegistry :: DielectricRegistry()
{
  the[DielectricKey(string("SiN"),Params())] = new SiN();
  the[DielectricKey(string("SiC"),Params())] = new SiC();
  the[DielectricKey(string("Al"),Params())] = new Al();
  the[DielectricKey(string("Au"),Params())] = new Au();
}

DielectricRegistry :: ~DielectricRegistry()
{
  for(DielectricRegistryMap::iterator i = the.begin(); i != the.end(); ++i) {
    delete i->second;
  }
}

