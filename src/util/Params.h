#ifndef PARAMS_H
#define PARAMS_H

#include "MathUtils.h"
#include <string>
#include <unordered_map>


class ParamValue {
public:

  enum ParamType {
    ParamNULL = 0,
    ParamString,
    ParamFloat,
    ParamInt,
    ParamComplex
  };

  virtual size_t hash() const { return 0; }
  virtual operator bool() const { return false; }
  
 ParamValue() : type(0) {}
  
 ParamValue(ParamType type) : type(type) {}
  int type;
  int i;
  double f;
  Complex z;
  std::string s; 
};

class ParamComplex : public ParamValue {
public:
 ParamComplex(const Complex &val) : ParamValue(ParamValue::ParamComplex) { z = val;}
  operator bool() const { return real(z) || imag(z); }
  size_t hash() const { return std::hash<double>()(real(z)) ^ std::hash<double>()(imag(z)); }
};


class ParamInt : public ParamValue {
public:
 ParamInt(int val) : ParamValue(ParamValue::ParamInt) { i = val;}
  operator bool() const { return i; }
  size_t hash() const { return i; }
};

class ParamFloat: public ParamValue {
public:
 ParamFloat(double val) : ParamValue(ParamValue::ParamFloat) { f = val;}
  operator bool() const { return f; }
  size_t hash() const { return f; }
};

class ParamString: public ParamValue {
public:
 ParamString(const std::string &val) : ParamValue(ParamValue::ParamString) { s = val;}
  operator bool() const { return !(s.empty()); }
  size_t hash() const { return std::hash<std::string>()(s); }
};

struct ParamMap
{
  typedef std::unordered_map<std::string, ParamValue> type;
};

class Params {
public:
  Params();
  static ParamValue parseParam(const std::string &name);
  ParamMap::type params;

  bool operator==(const Params& other) const;
  bool operator!=(const Params& other) const;
  ParamValue operator[](const std::string &s) const;
  ParamValue &operator[](const std::string &s);
};


size_t mapHasher(const std::unordered_map<std::string, ParamValue> &m);
size_t paramsHasher(const Params &m);
extern ParamValue paramNULLInst;

#endif
