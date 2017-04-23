#ifndef OPTIONS_H
#define OPTIONS_H

#include "Util.h"

class Options {
 public:
  Options();

  Real gapNormalization;
  Real tol;
  Real abstol;
  Real admissableThresh;
  Real interpolationDensity;
  Real interpolationPowerLaw;
  Real parameterizationThresh;
  Real acaTolerance;
  Real rkEps;
  Real inPlaneThresh;
  Real sharpThresh;
  Real concaveThresh;
  Real relativeCostAngleConcave;
  int minAdmissableSegmentSize;
  int minParameterizableSegmentSize;
  int minSegmentSize;
  Real relativeMinSegmentRadius;
  Real relativeMaxEdgeLength;
  Real farThresh6;
  Real farThresh10;
  Real farThresh15;
  Real farThresh21;
};

#endif
