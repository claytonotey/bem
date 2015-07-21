#ifdef _DEBUG
#undef _DEBUG
#include <omp.h>
#define _DEBUG
#else
#include <omp.h>
#endif

unsigned long referenceCountedMemUsage = 0;
omp_lock_t referenceCountingPointerLock;
bool bReferenceCountingPointerLockUnitialized = true;

void initReferenceCountingPointerLock()
{
  omp_init_lock(&referenceCountingPointerLock);
  bReferenceCountingPointerLockUnitialized = false;
}

void destroyReferenceCountingPointerLock()
{
  omp_destroy_lock(&referenceCountingPointerLock);
  bReferenceCountingPointerLockUnitialized = true;
}
