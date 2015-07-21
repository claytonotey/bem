#ifndef PROGRESS_H
#define PROGRESS_H

#include "Matrix.h"
#include "DenseMatrix.h"
#include "Rk.h"

class Progress {
 public:
  Progress();

  // time
  time_t start;

  // completion
  int lastPercent;

  void startRHSConstruction(int m);
  int rhsSize;
  int rhsDone;
  void finishedRHS(ComplexDenseMatrix &M);

  // matrix construction
  long matrixSize;
  long matrixDone;
  long dense;
  long rk;
  void startMatrixConstruction(int m, int n);
  void startedMatrix(ComplexMatrix &M);
  void finishedMatrix(ComplexDenseMatrix &M);
  void finishedMatrix(ComplexRk &M);
  void finishedMatrix(ComplexMatrix &M);
  
  // flux
  int numFluxVolumeElements;
  int numFluxVolumeElementsDone;
  void startFluxCalculation(int n);
  void finishedFluxVolumeElement();

};

#endif
