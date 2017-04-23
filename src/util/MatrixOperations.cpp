#include "Util.h"
#include "Matrix.h"
#include "DenseMatrix.h"

#define DEBUG 0

extern "C" {
  FORTRAN_RET_T 
  FORTRAN(zlaswp,ZLASWP)(const int& nrhs,
                         ComplexDouble *B,
                         const int& ldb,
                         const int& k1,
                         const int& k2,
                         int *ipiv,
                         const int& incx);

  FORTRAN_RET_T 
  FORTRAN(ztrsm,ZTRSM)(const int& side,
                       const int& uplo,
                       const int& transa,
                       const int& diag,
                       const int& m,
                       const int& n,
                       const ComplexDouble& alpha,
                       ComplexDouble *A,
                       const int& lda,
                       ComplexDouble *B,
                       const int& ldb);
};

template<>
  void ComplexDenseMatrix :: forwardSubstL(const ComplexDenseMatrix &A)
{
  FORTRAN(zlaswp,ZLASWP)(n,p,lda,1,m,A.piv,1);
  FORTRAN(ztrsm,ZTRSM)('L','L','N','U',m,n,ComplexDouble(1.,0.),A.p,A.lda,p,lda);
}

template<>
  void ComplexDenseMatrix :: forwardSubstU(const ComplexDenseMatrix &A)
{
  FORTRAN(ztrsm,ZTRSM)('R','U','N','N',m,n,ComplexDouble(1.,0.),A.p,A.lda,p,lda);
}

template<>
void ComplexDenseMatrix :: backwardSubstU(const ComplexDenseMatrix &A)
{
  FORTRAN(ztrsm,ZTRSM)('L','U','N','N',m,n,ComplexDouble(1.,0.),A.p,A.lda,p,lda);
}


extern "C"
{
  FORTRAN_RET_T 
  FORTRAN(zgemm,ZGEMM)(const int& transa,
                       const int& transb,
                       const int& m,
                       const int& n,
                       const int& k,
                       const ComplexDouble& alpha,
                       ComplexDouble* A,
                       const int& lda,
                       ComplexDouble* B,
                       const int& ldb,
                       const ComplexDouble& beta,
                       ComplexDouble* C,
                       const int& ldc);

};

template<>
  void ComplexDenseMatrix :: addproduct(const ComplexDenseMatrix &A, const ComplexDenseMatrix &B) 
{
  FORTRAN(zgemm,ZGEMM)('N','N',A.m,B.n,A.n,ComplexDouble(1.,0.),A.p,A.lda,B.p,B.lda,ComplexDouble(1.,0.),p,lda);
}

template<>
  void ComplexDenseMatrix :: subtractproduct(const ComplexDenseMatrix &A, const ComplexDenseMatrix &B) 
{
  FORTRAN(zgemm,ZGEMM)('N','N',A.m,B.n,A.n,ComplexDouble(-1.,0.),A.p,A.lda,B.p,B.lda,ComplexDouble(1.,0.),p,lda);
}

#undef DEBUG
