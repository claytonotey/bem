#ifndef COMPLEXVECTOR3PAIR_H
#define COMPLEXVECTOR3PAIR_H

#include "defs.h"
#include "ComplexVector3.h"

class ComplexVector3Pair {
public:
  ComplexVector3 first, second;

  ComplexVector3Pair() {}

  ComplexVector3Pair(const ComplexVector3 &first, const ComplexVector3 &second) {
    this->first = first;
    this->second = second;
  }
  
  inline ComplexVector3Pair operator-() {
    ComplexVector3Pair p;
    p.first = -first;
    p.second = -second;
    return p;
  }

  inline void operator+=(const ComplexVector3Pair &p) {
    first += p.first;
    second += p.second;
  }
  
  inline void operator-=(const ComplexVector3Pair &p) {
    first -= p.first;
    second -= p.second;
  }
  
  inline void operator*=(Real x) {
    first *= x;
    second *= x;
  }

  inline void operator*=(const Complex &x) {
    first *= x;
    second *= x;
  }
  
  inline friend ComplexVector3Pair operator*(Real x, const ComplexVector3Pair &p) {
    ComplexVector3Pair q(x*p.first,x*p.second);
    return q;
  }

  inline friend Real norm(const ComplexVector3Pair &p) {
    return norm(p.first) + norm(p.second);
  }
  
  inline friend ostream& operator<<(ostream& os, const ComplexVector3Pair &r) {
    os << r.first << ";" << r.second;
    return os;
  }
  
};

#endif
