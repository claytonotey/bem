#ifndef DEFS_H
#define DEFS_H

#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <limits>
#include <float.h>
using namespace std;

inline double real(double x)
{
  return x;
}

template<class R>
inline R norm(const R &x)
{
  return x * x;
}

template<class T>
inline T& min(T &x, T &y) {
  return (x<y)?x:y;
}

template<class T>
inline T& max(T &x, T &y) {
  return (x>y)?x:y;
}

/* CONS(a,b) should return ab, the concatenation
   of its arguments */
#if  __STDC__ || __APOGEE__ || WIN32
#define CONS(a,b) a##b
#else
#define CONS(a,b) a/**/b
#endif

#define FORTRAN_RET_T int
#define FORTRAN(name,NAME) CONS(name,_)

typedef unsigned long Key;
typedef long double LongReal;
typedef complex<LongReal> LongComplex;
typedef complex<double> ComplexDouble;
typedef complex<float> ComplexFloat;

//#define REAL_FLOAT

#ifdef REAL_FLOAT
typedef float Real;
typedef ComplexFloat Complex;
#define MIN_EXP FLT_MIN_EXP
#define MAX_EXP FLT_MAX_EXP
#else
typedef double Real;
typedef ComplexDouble Complex;
#define MIN_EXP DBL_MIN_EXP
#define MAX_EXP DBL_MAX_EXP
#endif

typedef long double LongReal;
typedef complex<long double> LongComplex;

#define Infinity (numeric_limits<Real>::max())
#define NegativeInfinity (numeric_limits<Real>::min())

inline bool isnan(const Complex &z)
{
  return isnan(real(z)) || isnan(imag(z)) || isinf(real(z)) || isinf(imag(z));
}

//using namespace __gnu_cxx;
//#define GCC_BEFORE_UNORDERED 1
#if WIN32 || GCC_BEFORE_UNORDERED
#include <ext/hash_map>
#include <ext/hash_set>
#define UnorderedMap hash_map
#define UnorderedSet hash_set
namespace __gnu_cxx
{
template<class T>
struct hash<T*>
{
  size_t operator()(T *s) const { return (size_t)s; }
};
}

#else
#include <unordered_map>
#include <unordered_set>
#define UnorderedMap std::unordered_map
#define UnorderedSet std::unordered_set
#endif

#endif
