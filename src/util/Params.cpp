#include "Params.h"
#include <iostream>

using namespace std;

ParamValue paramNULLInst;

size_t mapHasher(const std::unordered_map<std::string, ParamValue> &m)
{ 
  size_t ret = 0;
  for(std::unordered_map<std::string, ParamValue>::const_iterator i = m.begin(); i != m.end(); ++i) {
    ret ^= (std::hash<std::string>()(i->first) ^ (i->second.hash() << 1));
  }
  return ret;
}

size_t paramsHasher(const Params &m) 
{
  return mapHasher(m.params);
}

ParamValue Params :: parseParam(const std::string &s)
{
  try {
    Complex z = parseComplex(s);
    return ParamComplex(z);
  } catch(std::invalid_argument e) {
    try {
      int i = std::stoi(s);
      return ParamInt(i);
    } catch(std::invalid_argument e) {
      try {
        double f = std::stof(s);
        return ParamFloat(f);
      } catch(std::invalid_argument e) {
        return ParamString(s);
      }
    }
  }
  return paramNULLInst;
}

Params :: Params() {}

ParamValue Params :: operator[](const std::string &s) const {
  ParamMap::type::const_iterator i = params.find(s);
  if(i != params.end()) {
    return (i->second);
  } else {
    return paramNULLInst;
  }
}

ParamValue &Params :: operator[](const std::string &s) {
  return params[s];
}


bool Params::operator==(const Params& other) const {
  return params == other.params;
}

bool Params :: operator!=(const Params& other) const {
  return !(*this == other);
}
