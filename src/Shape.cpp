#include "Util.h"
#include "Shape.h"

Shape :: Shape()
{

}

void Shape :: transform(Transform &tform)
{
  for(int i=0; i<coord.size(); i++) {
    coord[i] = tform.transform(coord[i]);
  }
}

ostream& operator<<(ostream& os, const Shape &s)
{

  int fisize = s.faceIndex.size();
  os << "Face Index:\n";
  for(int i=0; i<fisize; i++) {
    os << s.faceIndex[i] << "\n";
  }

  int csize = s.coord.size();
  os << "\nCoordinates:\n";
  for(int i=0; i<csize; i++) {
    os << s.coord[i] << "\n";
  }
  
  os << "\n";

  return os;
}
