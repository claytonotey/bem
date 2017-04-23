#ifndef DIELECTRICS_H
#define DIELECTRICS_H

#include "Dielectric.h"
#include <string>
#include <utility>
#include "Params.h"

#include <string>
#include <map>
#include <utility>
#include <functional>

typedef std::pair< std::string, Params> DielectricKey;

struct pair_hash {
  std::size_t operator () (const DielectricKey &p) const 
  {
    std::hash<std::string> hash1;
    std::size_t h1 = hash1(p.first);
    std::size_t h2 = paramsHasher(p.second);
    return h1 ^ (h2<<1);  
  }
};

typedef std::unordered_map<DielectricKey, Dielectric*, pair_hash> DielectricRegistryMap;

class DielectricRegistry {
 public:
  DielectricRegistry();
  ~DielectricRegistry();
  DielectricRegistryMap the;
};


class Dielectrics {
 public:
  ~Dielectrics();
  static Dielectric *getDielectric(const std::string &s, const Params &params);
  static Dielectric *getDielectric(const std::string &s);
  static DielectricRegistry registry;
};


#endif

