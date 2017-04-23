#ifndef COMPLEXPAIR_H
#define COMPLEXPAIR_H

class ComplexPair {
 public:
  Complex first, second;

  ComplexPair() {}

  ComplexPair(const Complex &first, const Complex &second) {
    this->first = first;
    this->second = second;
  }

  inline ComplexPair operator-() {
    ComplexPair p;
    p.first = -first;
    p.second = -second;
    return p;
  }
  
  inline void operator+=(const ComplexPair &p) {
    first += p.first;
    second += p.second;
  }
  
  inline void operator-=(const ComplexPair &p) {
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
  
  inline friend ComplexPair operator*(Real x, const ComplexPair &p) {
    ComplexPair q(x*p.first,x*p.second);
    return q;
  }
  
  inline friend Real norm(const ComplexPair &p) {
    return norm(p.first) + norm(p.second);
  }
  
  inline friend ostream& operator<<(ostream& os, const ComplexPair &r) {
    os << r.first << ";" << r.second;
    return os;
  }
  
};

#endif
