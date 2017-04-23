#include "VolumeMesh.h"

ostream& operator<<(ostream& os, const VolumeMesh &m)
{
  os << "v = [\n";
  for(int i=0; i<m.tetrahedra.n; i++) {
    os << m.tetrahedra[i] << "\n";
  }
  os << "]\n;";
  return os;
}
