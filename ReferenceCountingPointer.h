#ifndef REFERENCECOUNTINGPOINTER_H
#define REFERENCECOUNTINGPOINTER_H

#include <stdlib.h>
#include <iostream>
using namespace std;

#ifdef _DEBUG
#undef _DEBUG
#include <omp.h>
#define _DEBUG
#else
#include <omp.h>
#endif

extern omp_lock_t referenceCountingPointerLock;
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
      cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }
      
    omp_set_lock(&referenceCountingPointerLock);
    p = q.p;
    pp = q.pp;
    count = q.count;
    size = q.size;
    if(p) (*count)++;
    omp_unset_lock(&referenceCountingPointerLock);
  }

  ReferenceCountingPointer(const ReferenceCountingPointer<T> &q, int offset) 
  {
    if(bReferenceCountingPointerLockUnitialized) {
      cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }

    omp_set_lock(&referenceCountingPointerLock);
    p = q.p;
    pp = q.pp + offset;
    count = q.count;
    size = q.size;
    if(p) (*count)++;
    omp_unset_lock(&referenceCountingPointerLock);
  }

  ~ReferenceCountingPointer()
  {
    if(bReferenceCountingPointerLockUnitialized) {
      cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }

    omp_set_lock(&referenceCountingPointerLock);
    if(count != NULL) {
      (*count)--;
      if((*count) == 0) destroy();
    }
    omp_unset_lock(&referenceCountingPointerLock);
  }

  inline ReferenceCountingPointer<T>& operator=(const ReferenceCountingPointer<T> &q)
  {
    if(this == &q)
      return *this;

    if(bReferenceCountingPointerLockUnitialized) {
      cerr << "ReferenceCountingPointer is not initialized: call initReferenceCountingPointerLock() before use\n";
      abort();
    }

    omp_set_lock(&referenceCountingPointerLock);
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
    omp_unset_lock(&referenceCountingPointerLock);
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
