#include "Progress.h"
#include <time.h>

Progress :: Progress()
{
  start = time(NULL);
}

void Progress :: startRHSConstruction(int m)
{
  rhsSize = m;
  rhsDone = 0;
  lastPercent = -1;
}

void Progress :: finishedRHS(ComplexDenseMatrix &M)
{
  rhsDone += M.getNumRows();
  int percent = (int)(1000. * ((Real)rhsDone / (Real)rhsSize));

  if(lastPercent != percent) {
    long elapsed = (long)(time(NULL)-start);

    fprintf(stderr,"\rRHS Construction: completion: %3.1f%%, elapsed: %lds    ", percent>1000?100.:percent/10., elapsed);
    lastPercent = percent;
    if(percent == 1000)
      fprintf(stderr,"\n");
    fflush(stderr);
  }
}

void Progress :: startMatrixConstruction(int m, int n)
{
  matrixSize = m * n;
  matrixDone = 0;
  dense = 0;
  rk = 0;
  lastPercent = -1;
}

void Progress :: startedMatrix(ComplexMatrix &M)
{
}

void Progress :: finishedMatrix(ComplexDenseMatrix &M)
{
  dense += M.getNumRows() * M.getNumCols();
  finishedMatrix((ComplexMatrix&)M);
}

void Progress :: finishedMatrix(ComplexRk &M)
{
  rk += M.getNumRows() * M.getNumCols();
  finishedMatrix((ComplexMatrix&)M);
}

void Progress :: finishedMatrix(ComplexMatrix &M)
{
  matrixDone += M.getNumRows() * M.getNumCols();
  int percent = (int)(1000. * ((Real)matrixDone / (Real)matrixSize));

  if(lastPercent != percent) {
    int densepercent = (int)(1000. * ((Real)dense / (Real)matrixDone));
    int rkpercent = (int)(1000. * ((Real)rk / (Real)matrixDone));
    long elapsed = (long)(time(NULL)-start);

    fprintf(stderr,"\rMatrix Construction: %3.1f%% dense, %3.1f%% rk, completion: %3.1f%%, mem: %ldM, elapsed: %lds      ", densepercent>1000?100.:densepercent/10., rkpercent>1000?100.:rkpercent/10., percent>1000?100.:percent/10., referenceCountedMemUsage>>20, elapsed);
    lastPercent = percent;
    if(percent == 1000)
      fprintf(stderr,"\n");
    fflush(stderr);
  }
}

void Progress :: startFluxCalculation(int n)
{
  numFluxVolumeElements = n;
  numFluxVolumeElementsDone = 0;
  lastPercent = -1;
}

void Progress :: finishedFluxVolumeElement()
{
  numFluxVolumeElementsDone++;

  long elapsed = (long)(time(NULL)-start);
  
  fprintf(stderr,"\rFlux Calculation: completion: %d/%d elapsed: %lds      ", numFluxVolumeElementsDone, numFluxVolumeElements, elapsed);
  if(numFluxVolumeElementsDone == numFluxVolumeElements)
    fprintf(stderr,"\n");
  fflush(stderr);
}
