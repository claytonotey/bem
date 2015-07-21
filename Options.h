#ifndef OPTIONS_H
#define OPTIONS_H

#include "defs.h"

class Options {
 public:

  static Real gapNormalization;
  static Real tol;
  static Real abstol;
  static Real admissableThresh;
  static Real interpolationDensity;
  static Real interpolationPowerLaw;
  static Real parameterizationThresh;
  static Real acaTolerance;
  static Real rkEps;
  static Real inPlaneThresh;
  static Real sharpThresh;
  static Real concaveThresh;
  static Real relativeCostAngleConcave;
  static int minAdmissableSegmentSize;
  static int minParameterizableSegmentSize;
  static int minSegmentSize;
  static Real relativeMinSegmentRadius;
  static Real relativeMaxEdgeLength;
  static Real farThresh6;
  static Real farThresh10;
  static Real farThresh15;
  static Real farThresh21;
};

#endif
