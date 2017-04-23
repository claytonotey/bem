#include "Util.h"
#include "DenseMatrix.h"
#include "MatrixOperations.h"

extern "C" {
  FORTRAN_RET_T
  FORTRAN(zgeev, ZGEEV) (const char&, //JOBVL
                          const char&, //JOBVR
                          const int &, //N
                          ComplexDouble*, //A
                          const int&, //LDA
                          ComplexDouble*, //W
                          ComplexDouble*, //VL
                          const int&, //LDVL
                          ComplexDouble*, //VR
                          const int&, //LDVR
                          ComplexDouble*, //WORK
                          const int&, //LWORK
                          double*, //RWORK
                          int &); //INFO
}

void eig(ComplexDenseMatrix &A, Array<Complex> &w)
{
  w.resize(A.n);
  Array<double> rwork(2*A.n);
  Array<ComplexDouble> work(1);
  int worksize = -1;
  int info;
  FORTRAN(zgeev,ZGEEV)('N','N',A.n,A.p,A.lda,w.p,NULL,1,NULL,1,work.p,worksize,rwork.p,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgeev,ZGEEV)('N','N',A.n,A.p,A.lda,w.p,NULL,1,NULL,1,work.p,worksize,rwork.p,info);
  if(info)
    abort();
}
    
extern "C" {
  FORTRAN_RET_T
  FORTRAN(zgetrf, ZGETRF) (const int&, // M
                            const int&, // N
                            ComplexDouble*, // A
                            const int&, // LDA
                            int *, //IPIV
                            int &); //INFO

  FORTRAN_RET_T
  FORTRAN(zgetri, ZGETRI) (const int&, // M
                            ComplexDouble*, // A
                            const int&, // LDA
                            int *, //IPIV
                            ComplexDouble *, // WORK
                            const int&, // LWORK
                            int &); //INFO

}

void lu(DenseMatrix<double,ComplexDouble> &A, int *ipiv)
{
  int n = A.n;
  int info;
  FORTRAN(zgetrf,ZGETRF)(n,n,A.p,A.lda,ipiv,info);
}

void inverse(DenseMatrix<double,ComplexDouble> &A)
{
  int n = A.n;
  Array<ComplexDouble> work(1);
  Array<int> ipiv(n);
  int worksize = -1;
  int info;
  FORTRAN(zgetrf,ZGETRF)(n,n,A.p,A.lda,ipiv.p,info);
  FORTRAN(zgetri,ZGETRI)(n,A.p,A.lda,ipiv.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgetri,ZGETRI)(n,A.p,A.lda,ipiv.p,work.p,worksize,info);
}

extern "C" {
  FORTRAN_RET_T
  FORTRAN(zgesvd, ZGESVD) (const int&, // JOBU
                            const int&, // JOBVT
                            const int&, // M
                            const int&, // N
                            ComplexDouble*, // A
                            const int&, // LDA
                            double*, // S
                            ComplexDouble*, // U
                            const int&, // LDU
                            ComplexDouble*, // VT
                            const int&, // LDVT
                            ComplexDouble*, // WORK
                            const int&, // LWORK
                            double*, // RWORK
                            int&); //INFO

  FORTRAN_RET_T
  FORTRAN(zgesdd, ZGESDD) (const int&, // JOBZ
                           const int&, // M
                           const int&, // N
                           ComplexDouble*, // A
                           const int&, // LDA
                           double*, // S
                           ComplexDouble*, // U
                           const int&, // LDU
                           ComplexDouble*, // VT
                           const int&, // LDVT
                           ComplexDouble*, // WORK
                           const int&, // LWORK
                           double*, // RWORK
                           int*, //IWORK
                           int&); //INFO
}

void svd(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &U, Array<double> &S, DenseMatrix<double,ComplexDouble> &VT)
{
  int m = A.m;
  int n = A.n;
  int k = min(m,n);
  Array<ComplexDouble> work(1);
  Array<double> rwork(5 * k);
  int worksize = -1;
  int info;
  char jobU = 'S';
  char jobVT = 'S';
  if(U.n != k)
    U.resize(m,k);
  if(VT.m != k)
    VT.resize(k,n);
  if(S.n != k)
    S.resize(k);
  FORTRAN(zgesvd,ZGESVD)(jobU,jobVT,m,n,A.p,A.lda,S.p,U.p,U.lda,VT.p,VT.lda,work.p,worksize,rwork.p,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgesvd,ZGESVD)(jobU,jobVT,m,n,A.p,A.lda,S.p,U.p,U.lda,VT.p,VT.lda,work.p,worksize,rwork.p,info);
}

void svd2(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &U, Array<double> &S, DenseMatrix<double,ComplexDouble> &VT)
{
  int m = A.m;
  int n = A.n;
  int k = min(m,n);
  Array<ComplexDouble> work(1);
  Array<double> rwork(5*k*k + 7*k);
  Array<int> iwork(8*k);
  int worksize = -1;
  int info;
  if(U.n != k)
    U.resize(m,k);
  if(VT.m != k)
    VT.resize(k,n);
  if(S.n != k)
    S.resize(k);
  FORTRAN(zgesdd,ZGESDD)('S',m,n,A.p,A.lda,S.p,U.p,U.lda,VT.p,VT.lda,work.p,worksize,rwork.p,iwork.p,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgesdd,ZGESDD)('S',m,n,A.p,A.lda,S.p,U.p,U.lda,VT.p,VT.lda,work.p,worksize,rwork.p,iwork.p,info);
}

Real condition(ComplexDenseMatrix &A)
{
  Array<double> S;
  DenseMatrix<double,ComplexDouble> U;
  DenseMatrix<double,ComplexDouble> VT;
  svd(A,U,S,VT);
  Real minS = S[0];
  Real maxS = minS;
  for(int i=0; i<S.n; i++) {
    Real s = S[i];
    if(s < minS)
      minS = s;
    else if(s > maxS)
      maxS = s;
  }

  return (maxS / minS);
}

extern "C"
{
  FORTRAN_RET_T
  FORTRAN(zgeqrf, ZGEQRF) (const int&,
                            const int&,
                            ComplexDouble*,
                            const int&,
                            ComplexDouble*,
                            ComplexDouble*,
                            const int&,
                            int&);

  FORTRAN_RET_T
  FORTRAN(zungqr, ZUNGQR) (const int&,
                            const int&,
                            const int&,
                            ComplexDouble*,
                            const int&,
                            ComplexDouble*,
                            ComplexDouble*,
                            const int&,
                            int&);

  FORTRAN_RET_T
  FORTRAN(zunmqr, ZUNMQR) (const int& side,
                           const int& trans,
                           const int& m,
                           const int& n,
                           const int& k,
                           ComplexDouble *A,
                           const int& lda,
                           ComplexDouble *tau,
                           ComplexDouble *C,
                           const int& ldc,
                           ComplexDouble *work,
                           const int& lwork,
                           int& info);
}

void qr(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &Q, DenseMatrix<double,ComplexDouble> &R)
{
  int m = A.m;
  int n = A.n;
  int k = min(m,n);
  Array<ComplexDouble> tau(k);
  Array<ComplexDouble> work(1);
  int worksize = -1;
  int info;
  R.copyFrom(A,true);
  FORTRAN(zgeqrf,ZGEQRF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgeqrf,ZGEQRF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);

  // form Q from elementary reflectors
  Q.copyFrom(R,true);
  worksize = -1;
  FORTRAN(zungqr,ZUNGQR)(m,k,k,Q.p,Q.lda,tau.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zungqr,ZUNGQR)(m,k,k,Q.p,Q.lda,tau.p,work.p,worksize,info);
  if(k != n) Q.resize(m,k);

  // set lower triangular elements of R to zero
  for(int i=0; i<k; i++) {
    for(int j=0; j<i; j++) {
      R(i,j) = 0.;
    }
  }
  if(k != m) R.resize(k,n);
}

void qr(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau, DenseMatrix<double,ComplexDouble> &R)
{
  int m = A.m;
  int n = A.n;
  int k = min(m,n);
  Array<ComplexDouble> work(1);
  tau.resize(k);
  int worksize = -1;
  int info;
  R.copyFrom(A,true);
  FORTRAN(zgeqrf,ZGEQRF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgeqrf,ZGEQRF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);

  // Q contains elementary reflectors
  Q.copyFrom(R,true);
  if(k != n) Q.resize(m,k);

  // set lower triangular elements of R to zero
  for(int i=0; i<k; i++) {
    for(int j=0; j<i; j++) {
      R(i,j) = 0.;
    }
  }
  if(k != m) R.resize(k,n);
}

void lmq(DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau, DenseMatrix<double,ComplexDouble> &C)
{
  int m = C.m;
  int n = C.n;
  int k = tau.n;
  Array<ComplexDouble> work(1);
  int worksize = -1;
  int info;
  FORTRAN(zunmqr,ZUNMQR)('L','N',m,n,k,Q.p,Q.lda,tau.p,C.p,C.lda,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zunmqr,ZUNMQR)('L','N',m,n,k,Q.p,Q.lda,tau.p,C.p,C.lda,work.p,worksize,info);
}

extern "C"
{
  FORTRAN_RET_T
  FORTRAN(zgerqf, ZGERQF) (const int&,
                            const int&,
                            ComplexDouble*,
                            const int&,
                            ComplexDouble*,
                            ComplexDouble*,
                            const int&,
                            int&);

  FORTRAN_RET_T
  FORTRAN(zungrq, ZUNGRQ) (const int&,
                            const int&,
                            const int&,
                            ComplexDouble*,
                            const int&,
                            ComplexDouble*,
                            ComplexDouble*,
                            const int&,
                            int&);


  FORTRAN_RET_T
  FORTRAN(zunmrq, ZUNMRQ) (const int& side,
                           const int& trans,
                           const int& m,
                           const int& n,
                           const int& k,
                           ComplexDouble *A,
                           const int& lda,
                           ComplexDouble *tau,
                           ComplexDouble *C,
                           const int& ldc,
                           ComplexDouble *work,
                           const int& lwork,
                           int& info);
}

void rq(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &R, DenseMatrix<double,ComplexDouble> &Q)
{
  int m = A.m;
  int n = A.n;
  int k = min(m,n);
  Array<ComplexDouble> tau(k);
  Array<ComplexDouble> work(1);
  int worksize = -1;
  int info;
  R.copyFrom(A,true);
  FORTRAN(zgerqf,ZGERQF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgerqf,ZGERQF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);

  // form Q from elementary reflectors
  Q.copyFrom(R,true);
  worksize = -1;
  if(k != m) Q.resize(m-n,0,k,n);
  FORTRAN(zungrq,ZUNGRQ)(k,n,k,Q.p,Q.lda,tau.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zungrq,ZUNGRQ)(k,n,k,Q.p,Q.lda,tau.p,work.p,worksize,info);

  if(m <= n) {
    R.resize(0,n-m,m,m);
    for(int i=0; i<m; i++) {
      for(int j=0; j<i; j++) {
        R(i,j) = 0.;
      }
    }
  } else {
    for(int i=m-n; i<m; i++) {
      int inm = i+n-m;
      for(int j=0; j<inm; j++) {
        R(i,j) = 0.;
      }
    }
  }
}

void rq(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &R, DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau)
{
  int m = A.m;
  int n = A.n;
  int k = min(m,n);
  tau.resize(k);
  Array<ComplexDouble> work(1);
  int worksize = -1;
  int info;
  R.copyFrom(A,true);
  FORTRAN(zgerqf,ZGERQF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zgerqf,ZGERQF)(m,n,R.p,R.lda,tau.p,work.p,worksize,info);

  // Q contains elementary reflectors
  Q.copyFrom(R,true);
  if(k != m) Q.resize(m-n,0,k,n);

  if(m <= n) {
    R.resize(0,n-m,m,m);
    for(int i=0; i<m; i++) {
      for(int j=0; j<i; j++) {
        R(i,j) = 0.;
      }
    }
  } else {
    for(int i=m-n; i<m; i++) {
      int inm = i+n-m;
      for(int j=0; j<inm; j++) {
        R(i,j) = 0.;
      }
    }
  }
}

void rmq(DenseMatrix<double,ComplexDouble> &C, DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau)
{
  int m = C.m;
  int n = C.n;
  int k = tau.n;
  Array<ComplexDouble> work(1);
  int worksize = -1;
  int info;
  FORTRAN(zunmrq,ZUNMRQ)('R','N',m,n,k,Q.p,Q.lda,tau.p,C.p,C.lda,work.p,worksize,info);
  worksize = (int) real(work[0]);
  work.resize(worksize);
  FORTRAN(zunmrq,ZUNMRQ)('R','N',m,n,k,Q.p,Q.lda,tau.p,C.p,C.lda,work.p,worksize,info);
}

template<>
ComplexDenseMatrix eye(int n) 
{
  ComplexDenseMatrix I(n,n,0);
  for(int i=0; i<n; i++)
    I(i,i) = Complex(1.,0);
  return I;  
}

template<>
void ComplexDenseMatrix :: LU()
{
  piv = ReferenceCountingPointer<int>(m);
  lu(*this,piv);
}


extern "C"
{
  FORTRAN_RET_T
  FORTRAN(zgemm, ZGEMM) (const int& transa,
                         const int& transb,
                         const int& m,
                         const int& n,
                         const int& k,
                         const ComplexDouble &alpha,
                         ComplexDouble* a,
                         const int& lda,
                         ComplexDouble* b,
                         const int& ldb,
                         const ComplexDouble &beta,
                         ComplexDouble* c,
                         const int& ldc);
}

void multiply_A_Bdagger(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &B, DenseMatrix<double,ComplexDouble> &C)
{
  int m = A.m;
  int n = B.m;
  int k = A.n;
  C.resize(A.m,B.m);
  FORTRAN(zgemm,ZGEMM)('N','C',m,n,k,1.0,A.p,A.lda,B.p,B.lda,0.0,C.p,C.lda);
}

void multiply_Adagger_B(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &B, DenseMatrix<double,ComplexDouble> &C)
{
  int m = A.n;
  int n = B.n;
  int k = A.m;
  C.resize(A.n,B.n);
  FORTRAN(zgemm,ZGEMM)('C','N',m,n,k,1.0,A.p,A.lda,B.p,B.lda,0.0,C.p,C.lda);
}

void symm(DenseMatrix<double,ComplexDouble> &A)
{
  for(int i=0; i<A.m; i++) {
    for(int j=0; j<i; j++) {
      A(i,j) = A(i,j) + conj(A(j,i));
      A(j,i) = conj(A(i,j));
    }
    A(i,i) = real(A(i,i)) + real(A(i,i));
  }
}

extern "C"
{
  FORTRAN_RET_T
  FORTRAN(zgetrs, ZGETRS) (const int& trans,
                           const int& n,
                           const int& nb,
                           ComplexDouble* a,
                           const int& lda,
                           int* ipiv,
                           ComplexDouble* b,
                           const int& ldb,
                           int &info);
}

template<>
void solveLU(ComplexDenseMatrix &A, ComplexDenseMatrix &B)
{
  int info;
  FORTRAN(zgetrs,ZGETRS)('N',A.n,B.n,A.p,A.lda,A.piv,B.p,B.lda,info);
}

template<>
void ComplexDenseMatrix :: output(ostream &os)
{
  os << "[\n";
  for(int i = 0; i < m; i++) {
    os << "[ ";
    for(int j = 0; j < n; j++) {
      Complex &z = (*this)(i,j);
      os << real(z) << "+" << imag(z) << "i ";
    }
    os << "]\n";
  }
  os << "]";
}
