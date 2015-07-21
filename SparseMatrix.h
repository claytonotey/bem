#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "SparseBLAS.h"

using namespace SparseBLAS;
typedef blas_sparse_matrix SparseMatrixDef;
#define SparseMatrixVectorMultiply(A,x,xstride,y,ystride) BLAS_xusmv<T>(blas_no_trans,1.,A,x,xstride,y,ystride)
#define SparseMatrixBegin(m,n) BLAS_xuscr_begin<T>(m,n);
#define SparseMatrixInsert(A,a,i,j) BLAS_xuscr_insert_entry<T>(A,a,i,j)
#define SparseMatrixEnd(A) BLAS_xuscr_end(A)
#define SparseMatrixDestroy(A) BLAS_usds(A)

template<class T>
class SparseMatrix {

 public:

  SparseMatrix(int m, int n) {
    def = SparseMatrixBegin(m,n);
  }

  ~SparseMatrix() {
    SparseMatrixDestroy(def);
  }

  inline void insert(int i, int j, const T &a) {
    SparseMatrixInsert(def,a,i,j);
  }
  
  inline void end() {
    SparseMatrixEnd(def);
  }

  inline void MatrixVectorMultiply(const T *x, int xstride, T *y, int ystride) {
    SparseMatrixVectorMultiply(def,x,xstride,y,ystride);
  }

  SparseMatrixDef def;
  
};

#endif
