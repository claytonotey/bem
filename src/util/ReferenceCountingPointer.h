#ifndef REFERENCECOUNTINGPOINTER_H
#define REFERENCECOUNTINGPOINTER_H

#include <iostream>
#include <stdlib.h>


#if HAVE_OMP
#include <omp.h>
#define OmpInit   omp_init_lock(&referenceCountingPointerLock);
#define OmpDestroy omp_destroy_lock(&referenceCountingPointerLock);
#define OmpLock   omp_init_lock(&referenceCountingPointerLock);
#define OmpUnlock   omp_init_lock(&referenceCountingPointerLock);
extern omp_lock_t referenceCountingPointerLock;
#else
#define OmpInit
#define OmpDestroy
#define OmpLock 
#define OmpUnlock 
#endif

extern unsigned long referenceCountedMemUsage;
extern bool bReferenceCountingPointerLockUnitialized;

void initReferenceCountingPointerLock();
void destroyReferenceCountingPointerLock();

template<class T>
class ReferenceCountingPointer {

 public:

  ReferenceCountingPointer() 
  {
    p = NULL;
    pp = NULL;
    count = NULL;
    size = 0;
  }

  ReferenceCountingPointer(int size)
  { 
    count = new int;
    (*count) = 1;
    this->size = size;
    if(size) {
      p = (T*)malloc(size*sizeof(T));
      pp = p;
    } else {
      p = NULL;
      pp = NULL;
    }
    referenceCountedMemUsage += size * sizeof(T);
  }
  
  ReferenceCountingPointer(const ReferenceCountingPointer<T> &q)
  {
    if(bReferenceCountingPointerLockUnitialized) {
      std::cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }
      
    OmpLock
    p = q.p;
    pp = q.pp;
    count = q.count;
    size = q.size;
    if(p) (*count)++;                               
    OmpUnlock
  }

  ReferenceCountingPointer(const ReferenceCountingPointer<T> &q, int offset) 
  {
    if(bReferenceCountingPointerLockUnitialized) {
      std::cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }

    OmpLock
    p = q.p;
    pp = q.pp + offset;
    count = q.count;
    size = q.size;
    if(p) (*count)++;
    OmpUnlock
  }

  ~ReferenceCountingPointer()
  {
    if(bReferenceCountingPointerLockUnitialized) {
      std::cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }

    OmpLock
    if(count != NULL) {
      (*count)--;
      if((*count) == 0) destroy();
    }
    OmpUnlock
  }

  inline ReferenceCountingPointer<T>& operator=(const ReferenceCountingPointer<T> &q)
  {
    if(this == &q)
      return *this;

    if(bReferenceCountingPointerLockUnitialized) {
      std::cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }

    OmpLock
    if(p && p == q.p) {
      pp = q.pp;
    } else {
      if(count != NULL) {
        (*count)--;
        if((*count) == 0) destroy();
      }
      p = q.p;
      pp = q.pp;
      count = q.count;
      size = q.size;
      if(count) (*count)++;
    }
    OmpUnlock
    return *this;
  }
  
  inline const T *operator+(int i) const
  {
    return pp + i;
  }

  inline T *operator+(int i)
  {
    return pp + i;
  }

  inline T& operator[](int i) const
  {
    return pp[i];
  }

  inline T& operator[](int i)
  {
    return pp[i];
  }

  inline operator bool() const
  {
    return (p!=NULL);
  }

  inline operator T*() const
  {
    return pp;
  }
  
  protected:
  
  inline void destroy()
  {
    if(p) {
      free(p);
      p = NULL;
      pp = NULL;
      referenceCountedMemUsage -= size * sizeof(T);
    }
    if(count) {
      delete count;
      count = NULL;
    }
  }
  
  T *p;
  T *pp;
  int *count;
  int size;

};

#endif
