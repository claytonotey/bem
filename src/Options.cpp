#include "Options.h"

Options :: Options()
{
  gapNormalization = 1.0;
  tol = 1e-7;
  abstol = 1e-16;
  admissableThresh = 1.5;
  interpolationDensity = 8.0;
  interpolationPowerLaw = 0.5;
  parameterizationThresh = 4.0;
  acaTolerance = 1e-8;
  rkEps = 3e-4;
  inPlaneThresh = 1e-4;
  sharpThresh = 2.0;
  concaveThresh = 0.1;
  relativeCostAngleConcave = 1.0;
  minSegmentSize = 32;
  minAdmissableSegmentSize = 24;
  minParameterizableSegmentSize = 32;
  relativeMinSegmentRadius = 0.0;
  relativeMaxEdgeLength = 0.1;
  farThresh6 = 12.0;
  farThresh10 = 8.0;
  farThresh15 = 4.0;
  farThresh21 = 2.0;
}
