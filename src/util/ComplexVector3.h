#ifndef COMPLEXVECTOR3_H
#define COMPLEXVECTOR3_H

#include "Vector3.h"

#include <iostream>
using namespace std;

class ComplexVector3 {
 public:
  Complex x,y,z;

  ComplexVector3() {}

  ComplexVector3(const Complex &x, const Complex &y, const Complex &z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }
  
  inline void operator*=(const Real c) {
    x *= c;
    y *= c;
    z *= c;
  }
  
  inline void operator*=(const Complex &c) {
    x *= c;
    y *= c;
    z *= c;
  }
  
  inline void operator+=(const ComplexVector3 &v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }

  inline void operator-=(const ComplexVector3 &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }
  
  inline ComplexVector3 operator-() {
    ComplexVector3 result(-x,-y,-z);
    return result;
  }

  inline friend ComplexVector3 operator/(const Real c, const ComplexVector3 &v) {
    ComplexVector3 result(c/v.x,c/v.y,c/v.z);
    return result;
  }

  inline friend ComplexVector3 operator*(const ComplexVector3 &v, const Real c) {
    ComplexVector3 result(c*v.x,c*v.y,c*v.z);
    return result;
  }
  
  inline friend ComplexVector3 operator*(const Real c, const ComplexVector3 &v) {
    ComplexVector3 result(c*v.x,c*v.y,c*v.z);
    return result;
  }
  
  inline friend ComplexVector3 operator*(const ComplexVector3 &v, const Complex &c) {
    ComplexVector3 result(c*v.x,c*v.y,c*v.z);
    return result;
  }

  inline friend ComplexVector3 operator*(const Complex &c, const ComplexVector3 &v) {
    ComplexVector3 result(c*v.x,c*v.y,c*v.z);
    return result;
  }
  
  inline friend ComplexVector3 operator+(const ComplexVector3 &v1, const ComplexVector3 &v2) {
    ComplexVector3 result(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
    return result;
  }

  inline friend ComplexVector3 operator+(const ComplexVector3 &v1, const Vector3 &v2) {
    ComplexVector3 result(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
    return result;
  }

  inline friend ComplexVector3 operator+(const Vector3 &v1, const ComplexVector3 &v2) {
    ComplexVector3 result(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
    return result;
  }

  inline friend ComplexVector3 operator-(const ComplexVector3 &v1, const ComplexVector3 &v2) {
    ComplexVector3 result(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
    return result;
  }
  
  // cross product
  inline friend ComplexVector3 operator*(const ComplexVector3 &v1, const ComplexVector3 &v2) {
    ComplexVector3 result(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x);
    return result;
  }
  
  inline friend ComplexVector3 operator*(const Vector3 &v1, const ComplexVector3 &v2) {
    ComplexVector3 result(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x);
    return result;
  }
  
  inline friend ComplexVector3 operator*(const ComplexVector3 &v1, const Vector3 &v2) {
    ComplexVector3 result(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x);
    return result;
  }
  
  inline friend ComplexVector3 conj(const ComplexVector3 &v) {
    ComplexVector3 result(conj(v.x),conj(v.y),conj(v.z));
    return result;
  }

  inline void zero() {
    x = y = z = Complex(0,0);
  }

  inline friend bool isnan(const ComplexVector3 &v) {
    return isnan(v.x) || isnan(v.y) || isnan(v.z);
  }
  
  inline friend Complex dot(const Vector3 &v1, const ComplexVector3 &v2) {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  }
  
  inline friend Complex dot(const ComplexVector3 &v1, const ComplexVector3 &v2) {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  }
  
  inline friend Vector3 real(const ComplexVector3 &v) {
    Vector3 result(real(v.x),real(v.y),real(v.z));
    return result;
  }

  inline friend Real norm(const ComplexVector3 &v) {
    return norm(v.x) + norm(v.y) + norm(v.z);
  }
  
  inline friend ostream& operator<<(ostream& os, const ComplexVector3 &v) {
    os << "[" << v.x << " " << v.y << " " << v.z << "]" ;
    return os;
  }

};

inline ComplexVector3 operator*(const Vector3 &v, const Complex &c) {
  ComplexVector3 result(c*v.x,c*v.y,c*v.z);
  return result;
}
  
inline ComplexVector3 operator*(const Complex &c, const Vector3 &v) {
  ComplexVector3 result(c*v.x,c*v.y,c*v.z);
  return result;
}

#endif

