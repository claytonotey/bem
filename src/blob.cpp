#include "ComplexVector3.h"
#include "HydrodynamicDielectric.h"
#include "Vector3.h"
#include "defs.h"
#include "MathUtils.h"
#include "Impedance.h"
#include "GreensFunction.h"
#include "RWG.h"
#include "ShapeParser.h"
#include "ShapeMesh.h"
#include "VolumeMesh.h"
#include <time.h>
#include "TriangleQuadrature.h"
#include "TetrahedronQuadrature.h"
#include "AnalyticIntegrals.h"
#include "Options.h"
#include "MatrixOperations.h"
#include "Source.h"
#include "Output.h"
#include "Flux.h"
#include "FrequencySweep.h"
#include "SI.h"
#include <xmmintrin.h>

int main(int argc, char **argv)
{
  Au diel;
  for(int i=0; i<100; i++)  {
    Real feV = 0.5 + 5.5 * i/(float)100;
    Complex eps = diel.eps(feV);
    Real n = real(sqrt(eps));
    Real k = imag(sqrt(eps));
    cout << feV << "\t" << n << "\t" << k << "\n";
  }
}
