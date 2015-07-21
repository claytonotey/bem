#ifndef ARRAY_H
#define ARRAY_H

#include <string.h>
using namespace std;
#include <iostream>

// a wrapper for a pointer with a destructor
template<typename T>
class Array {
 public:

  Array() { 
    p = NULL; 
  }

  Array(int n) {
    this->p = new T[n];
    this->n = n;
  }

  ~Array() {
    if(p) delete [] p;
  }

  void resize(int n) {
    if(p) delete [] p;
    p = new T[n];
    this->n = n;
  }

  inline void zero() {
    memset(p,0,n*sizeof(T));
  }

  inline int operator-(T *p) {
    return this->p - p;
  }

  inline T& operator[](int i) {
    return p[i];
  }

  inline T& operator[](int i) const {
    return p[i];
  }

  inline operator T*() const {
    return p;
  }

  inline T *operator+(int i) {
    return p + i;
  }

  friend ostream& operator<<(ostream& os, Array &v) {
    for(int i = 0; i < v.n; i++)
      os << v[i] << "\n";
    return os;
  }

  int n;
  T *p;
};

#endif
