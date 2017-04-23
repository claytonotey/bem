#include "ReferenceCountingPointer.h"

#if HAVE_OMP
omp_lock_t referenceCountingPointerLock;
#endif
unsigned long referenceCountedMemUsage = 0;

bool bReferenceCountingPointerLockUnitialized = true;

void initReferenceCountingPointerLock()
{
  OmpInit
  bReferenceCountingPointerLockUnitialized = false;
}

void destroyReferenceCountingPointerLock()
{
  OmpDestroy
  bReferenceCountingPointerLockUnitialized = true;
}
