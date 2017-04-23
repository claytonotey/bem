#ifndef VECTOR3_H
#define VECTOR3_H

#include "Util.h"
#include <iostream>
#include <exception>

class Vector3Exception: public std::exception
{
 public:
  Vector3Exception() {};
  virtual const char* what() const throw()
  {
    return "Error parsing Vector3";
  }
};

class Vector3 {
 public:
  Real x,y,z;

  Vector3() {}

  Vector3(Real x, Real y, Real z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  Vector3(char *str) {
    char *s = str;
    
    while(*s == ' ') s++; 
    if(*s == 0) throw Vector3Exception();
    x = atof(s);
    while(*s != ' ') s++; 
    
    while(*s == ' ') s++; 
    if(*s == 0) throw Vector3Exception();
    y = atof(s);
    while(*s != ' ') s++; 
    
    while(*s == ' ') s++; 
    if(*s == 0) throw Vector3Exception();
    z = atof(s);
  }

  inline void operator+=(const Vector3 &v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }
  
  inline void operator-=(const Vector3 &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }

  inline void operator*=(const Real c) {
    x *= c;
    y *= c;
    z *= c;
  }

  inline friend Vector3 operator*(const Vector3 &v, const Real c) {
    Vector3 result(c*v.x,c*v.y,c*v.z);
    return result;
  }

  inline friend Vector3 operator*(const Real &c, const Vector3 &v) {
    Vector3 result(c*v.x,c*v.y,c*v.z);
    return result;
  }
  
  inline friend Vector3 operator+(const Vector3 &v1, const Vector3 &v2) {
    Vector3 result(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
    return result;
  }
  
  inline friend Vector3 operator-(const Vector3 &v1, const Vector3 &v2) {
    Vector3 result(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
    return result;
  }
  
  // cross product
  inline friend Vector3 operator*(const Vector3 &v1, const Vector3 &v2) {
    Vector3 result(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x);
    return result;
  }
  
  inline friend Real norm(const Vector3 &v) {
    return v.x*v.x + v.y*v.y + v.z*v.z;
  }
  
  inline friend Real dot(const Vector3 &v1, const Vector3 &v2) {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  }
  
  inline Vector3 scale(const Vector3 &v) {
    Vector3 result(x*v.x,y*v.y,z*v.z);
    return result;
  }

  inline friend bool isnan(const Vector3 &v) {
    return isnan(v.x) || isnan(v.y) || isnan(v.z);
  }
  
  inline void normalize() {
    Real normi = 1. / sqrt(x*x+y*y+z*z);
    x *= normi;
    y *= normi;
    z *= normi;
  }

  inline friend Vector3 expandCoords(Vector3 &vx, Vector3 &zz, Real x, Real y, Real z) {
    Vector3 r;
    Vector3 vy = zz * vx;
    r.x = x*vx.x + y*vy.x + z*zz.x;
    r.y = x*vx.y + y*vy.y + z*zz.y;
    r.z = x*vx.z + y*vy.z + z*zz.z;
    return r;
  }

  inline friend istream& operator >> (istream& is, Vector3 &v)
  {
    return is >> v.x >> v.y >> v.z;
  }

  inline friend ostream& operator<<(ostream& os, const Vector3 &v) {
    os << "[" << v.x << " " << v.y << " " << v.z << "]" ;
    return os;
  }
};


typedef Vector3 Point;

#endif
