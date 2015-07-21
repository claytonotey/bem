#ifndef FACEINDEX_H
#define FACEINDEX_H

#include <exception>
#include <iostream>
using namespace std;

class FaceIndex {
 public:
  FaceIndex(const char *str);

  int i[4];
};

ostream& operator<<(ostream& os, const FaceIndex &f);

class FaceIndexException: public exception
{
  virtual const char* what() const throw()
  {
    return "Error parsing FaceIndex";
  }
};

#endif
