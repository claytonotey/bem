#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <iostream>
#include <fstream>

template<class T>
void pushFileLines(const char *filename, vector<T> &v)
{  
  ifstream file(filename, ios_base::in)  ;
  if(!file) {
    abort();
  }
  T x;
  while(!file.eof()) {
    file >> x;
    if(!file.fail()) v.push_back(x);
  }
  file.close();
}

#endif
