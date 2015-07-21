#include "FaceIndex.h"
#include <stdlib.h>

FaceIndex :: FaceIndex(const char *str)
{
  char *s = (char*)str;
  for(int k=0;k<4;k++) {
    while(*s == ' ') s++; 
    if(*s == 0) throw FaceIndexException();
    i[k] = atoi(s);
    while(*s != ' ') s++; 
  }
  if(i[3] == 1) {
    int tmp = i[1];
    i[1] = i[2];
    i[2] = tmp;
  }
}

ostream& operator<<(ostream& os, const FaceIndex &f)
{
  os << "(" << f.i[0] << "," << f.i[1] << "," << f.i[2] << "," << f.i[3] << ")" ;
  return os;
}
