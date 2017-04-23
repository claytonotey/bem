#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

#include "MatrixBase.h"
#include "DenseMatrix.h"
#include "Rk.h"
#include "HMatrix.h"

#define DEBUG 0

// DenseMatrix LHS
template<class R, class T>
void DenseMatrix<R,T> :: copyFrom(const Rk<R,T> &A, bool bDeepCopy)
{
  init(A.getNumRows(),A.getNumCols());
  zero();
  addproduct(A.A,A.B);
}

template<class R, class T>
void DenseMatrix<R,T> :: copyFrom(const HMatrix<R,T> &A, bool bDeepCopy)
{
  int m = A.getNumRows();
  int n = A.getNumCols();
  init(m,n);
  zero();
  add(A);
}

template<class R, class T>
  void DenseMatrix<R,T> :: add(const DenseMatrix<R,T> &B) 
{
  if(m != B.m) abort();
  if(n != B.n) abort();

  for(int j=0; j<n; j++) {
    for(int i=0; i<m; i++) {
      (*this)(i,j) += B(i,j);
    }
  }
}

template<class R, class T>
  void DenseMatrix<R,T> :: add(const Rk<R,T> &A)
{
  addproduct(A.A,A.B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: add(const HMatrix<R,T> &A)
{
  HMatrix<R,T> C(*this,&(A.r),&(A.c));
  C.add(A);
}

template<class R, class T>
void DenseMatrix<R,T> :: subtract(const DenseMatrix<R,T> &B) 
{
  if(m != B.m) abort();
  if(n != B.n) abort();

  for(int j=0; j<n; j++) 
    for(int i=0; i<m; i++)
      (*this)(i,j) -= B(i,j);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtract(const Rk<R,T> &A)
{
  subtractproduct(A.A,A.B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtract(const HMatrix<R,T> &A)
{
  HMatrix<R,T> C(*this,&(A.r),&(A.c));
  C.subtract(A);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B) 
{
  T *c = p;
  T *a = A.p;
  T *b = B.p;
  for(int i=0; i<m; i++) {
    T *ai = a + i;
    for(int j=0; j<n; j++) {
      T *cij = c + (i+lda*j);
      T *bj = b + B.lda*j;
      for(int k=0; k<A.n; k++) {
        *cij += ai[A.lda*k] * bj[k];
      }
    }
  }
}

template<>
void ComplexDenseMatrix :: addproduct(const ComplexDenseMatrix &A, const ComplexDenseMatrix &B);

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B) 
{
  T *c = p;
  T *a = A.p;
  T *b = B.p;
  for(int i=0; i<m; i++) {
    T *ai = a + i;
    for(int j=0; j<n; j++) {
      T *cij = c + (i+lda*j);
      T *bj = b + B.lda*j;
      for(int k=0; k<A.n; k++) {
        *cij -= ai[A.lda*k] * bj[k];
      }
    }
  }
}

template<>
void ComplexDenseMatrix :: subtractproduct(const ComplexDenseMatrix &A, const ComplexDenseMatrix &B);

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)
{ 
  if(DEBUG) cerr << "+= ddr\n";
  DenseMatrix<R,T> C = A * B.A;
  addproduct(C,B.B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= ddr\n";
  DenseMatrix<R,T> C = A * B.A;
  subtractproduct(C,B.B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= ddh\n";
  HMatrix<R,T> C(A,NULL,&(B.r));
  addproduct(C,B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= ddh\n";
  HMatrix<R,T> C(A,NULL,&(B.r));
  subtractproduct(C,B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= drd\n";
  DenseMatrix<R,T> C = A.B * B;
  addproduct(A.A,C);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= drd\n";
  DenseMatrix<R,T> C = A.B * B;
  subtractproduct(A.A,C);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const Rk<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= drr\n";
  DenseMatrix<R,T> C = A.B * B.A;
  DenseMatrix<R,T> D = A.A * C;
  addproduct(D,B.B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const Rk<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= drr\n";
  DenseMatrix<R,T> C = A.B * B.A;
  DenseMatrix<R,T> D = A.A * C;
  subtractproduct(D,B.B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= drh\n";
  HMatrix<R,T> C(A,NULL,&(B.r));
  addproduct(C,B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= drh\n";
  HMatrix<R,T> C(A,NULL,&(B.r));
  subtractproduct(C,B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= dhd\n";
  HMatrix<R,T> C(B,&(A.c),NULL);
  addproduct(A,C);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= dhd " << getNumRows() << " " << getNumCols() << "\n";
  if(DEBUG) cerr << "-= dhd " << A.getNumRows() << " " << A.getNumCols() << "\n";
  if(DEBUG) cerr << "-= dhd " << B.getNumRows() << " " << B.getNumCols() << "\n";
  HMatrix<R,T> C(B,&(A.c),NULL);
  subtractproduct(A,C);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= dhr\n";
  HMatrix<R,T> C(B,&(A.c),NULL);
  addproduct(A,C);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= dhr\n";
  HMatrix<R,T> C(B,&(A.c),NULL);
  subtractproduct(A,C);
}

template<class R, class T>
  void DenseMatrix<R,T> :: addproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  HMatrix<R,T> C(*this,&(A.r),&(B.c));
  C.addproduct(A,B);
}

template<class R, class T>
  void DenseMatrix<R,T> :: subtractproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  HMatrix<R,T> C(*this,&(A.r),&(B.c));
  C.subtractproduct(A,B);
}

template<class R, class T>
void DenseMatrix<R,T> :: forwardSubstL(const DenseMatrix<R,T> &A)
{
  fprintf(stderr,"Undefined matrix operation: forwardSubstL\n");
  abort();
}

template<>
void ComplexDenseMatrix :: forwardSubstL(const ComplexDenseMatrix &A);

template<class R, class T>
void DenseMatrix<R,T> :: forwardSubstU(const DenseMatrix<R,T> &A)
{
  fprintf(stderr,"Undefined matrix operation: forwardSubstU\n");
  abort();
}

template<>
void ComplexDenseMatrix :: forwardSubstU(const ComplexDenseMatrix &A);

template<class R, class T>
void DenseMatrix<R,T> :: backwardSubstU(const DenseMatrix<R,T> &A)
{
  fprintf(stderr,"Undefined matrix operation: backwardSubstU\n");
  abort();
}

template<>
void ComplexDenseMatrix :: backwardSubstU(const ComplexDenseMatrix &A);

template<class R, class T>
  void DenseMatrix<R,T> :: forwardSubstL(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve L*X=B for singular (Rk) L\n");
  abort();
}

template<class R, class T>
  void DenseMatrix<R,T> :: forwardSubstL(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsl dh " << m << " " << n << " " << A.getNumRows() << " " << A.getNumCols() << "\n";
  HMatrix<R,T> B(*this,&(A.c),NULL);
  B.forwardSubstL(A);
}

template<class R, class T>
  void DenseMatrix<R,T> :: forwardSubstU(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve X*U=B for singular (Rk) U\n");
  abort();
}

template<class R, class T>
  void DenseMatrix<R,T> :: forwardSubstU(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsu dh\n";
  HMatrix<R,T> B(*this,NULL,&(A.r));
  B.forwardSubstU(A);
}

template<class R, class T>
  void DenseMatrix<R,T> :: backwardSubstU(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve U*X=B for singular (Rk) U\n");
  abort();  
}

template<class R, class T>
  void DenseMatrix<R,T> :: backwardSubstU(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "bsu dh\n";
  HMatrix<R,T> B(*this,&(A.c),NULL);
  B.backwardSubstU(A);
}

// Rk LHS
template<class R, class T>
  void Rk<R,T> :: truncate()
{
  fprintf(stderr,"Attempting to truncate undefined Rk matrix\n");
  abort();
}

template<class R, class T>
  void Rk<R,T> :: copyFrom(const HMatrix<R,T> &A, bool bDeepCopy)
{
  abort();
}

template<class R, class T>
  void Rk<R,T> :: copyFrom(const DenseMatrix<R,T> &A, bool bDeepCopy)
{
  DenseMatrix<R,T> B;
  B.copyFrom(A,true);
  init(B,eye<R,T>(A.n),eps);
}

template<class R, class T>
  void Rk<R,T> :: add(const DenseMatrix<R,T> &A)
{
  Rk<R,T> rk(A,eye<R,T>(A.n),eps);
  add(rk);
}

template<class R, class T>
  void Rk<R,T> :: add(const Rk<R,T> &rk2)
{
  Rk<R,T> rk1(*this);
  int m = rk1.m;
  int n = rk1.n;
  if(m != rk2.m || n != rk2.n)
    abort();
  int k = rk1.k + rk2.k;
  Real eps = min(rk1.eps, rk2.eps);

  init(m,n,k,eps);
  
  T *cp = A.p;
  T *ap = rk1.A.p;
  T *bp = rk2.A.p;

  for(int j=0; j<rk1.k; j++) {
    memcpy(cp,ap,m*sizeof(T));
    cp += m;
    ap += rk1.A.lda;
  }
  for(int j=0; j<rk2.k; j++) {
    memcpy(cp,bp,m*sizeof(T));
    cp += m;
    bp += rk2.A.lda;
  }
  
  cp = B.p;
  ap = rk1.B.p;
  bp = rk2.B.p;
  
  for(int j=0; j<n; j++) {
    memcpy(cp,ap,rk1.k*sizeof(T));
    cp += rk1.k;
    ap += rk1.B.lda;
    memcpy(cp,bp,rk2.k*sizeof(T));
    cp += rk2.k;
    bp += rk2.B.lda;
  }
  
  truncate();
}

template<class R, class T>
  void Rk<R,T> :: add(const HMatrix<R,T> &A)
{
  abort();
}

template<class R, class T>
  void Rk<R,T> :: subtract(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "-= rd\n";
  cout << m << " " << n << " " << k << " " << A.m << " " << A.n << "\n";
  Rk<R,T> rk(A,eye<R,T>(A.n),RK_DEFAULTEPS);
  subtract(rk);
  cout << k << "\n";
}

template<class R, class T>
  void Rk<R,T> :: subtract(const Rk<R,T> &rk2)
{
  Rk<R,T> rk1(*this);
  int m = rk1.m;
  int n = rk1.n;
  if(m != rk2.m || n != rk2.n)
    abort();
  int k = rk1.k + rk2.k;
  Real eps = min(rk1.eps, rk2.eps);

  init(m,n,k,eps);

  T *cp = A.p;
  T *ap = rk1.A.p;
  T *bp = rk2.A.p;

  for(int j=0; j<rk1.k; j++) {
    memcpy(cp,ap,m*sizeof(T));
    cp += m;
    ap += rk1.A.lda;
  }
  for(int j=0; j<rk2.k; j++) {
    memcpy(cp,bp,m*sizeof(T));
    cp += m;
    bp += rk2.A.lda;
  }
  
  cp = B.p;
  ap = rk1.B.p;
  
  for(int j=0; j<n; j++) {
    memcpy(cp,ap,rk1.k*sizeof(T));
    for(int i=0; i<rk2.k; i++) {
      B(rk1.k+i,j) = -rk2.B(i,j);
    }
    cp += k;
    ap += rk1.B.lda;
  }
  
  truncate();
}  

template<class R, class T>
  void Rk<R,T> :: subtract(const HMatrix<R,T> &A)
{
  abort();
}

template<class R, class T>
void Rk<R,T> :: addproduct(const Rk<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= rrr\n";
  Real eps = min(A.eps, B.eps);
  DenseMatrix<R,T> CA1 = A.B * B.A;
  DenseMatrix<R,T> CA2 = A.A * CA1;
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(CA2,B.B,eps);
    add(rk);
  } else {
    DenseMatrix<R,T> C;
    C.copyFrom(B.B,true);
    init(CA2,C,eps);
  }
}

template<class R, class T>
void Rk<R,T> :: subtractproduct(const Rk<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= rrr\n";
  Real eps = min(A.eps, B.eps);
  DenseMatrix<R,T> CA1 = A.B * B.A;
  DenseMatrix<R,T> CA2 = A.A * CA1;
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(CA2,B.B,eps);
    subtract(rk);
  } else {
    DenseMatrix<R,T> C = -B.B;
    init(CA2,C,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= rdd\n";
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(A,B,eps);
    add(rk);
  } else {
    DenseMatrix<R,T> C;
    DenseMatrix<R,T> D;
    C.copyFrom(A,true);
    D.copyFrom(B,true);
    init(C,D,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= rdd\n";
  if(this->A.p && this->B.p) {
    cout << A.m << " " << A.n << " " << B.m << " " << B.n << "\n";
    DenseMatrix<R,T> C = A * B;
    cout << norm(C) << "\n";
    cout << this->A.m << " " << this->A.n << " " << this->B.m << " " << this->B.n << "\n";
    DenseMatrix<R,T> D = this->A * this->B;
    cout << norm(D) << "\n";
    Rk<R,T> rk(A,B,eps);
    subtract(rk);
  } else {
    DenseMatrix<R,T> C = -A;
    DenseMatrix<R,T> D;
    D.copyFrom(B,true);
    init(C,D,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= rdr\n";
  DenseMatrix<R,T> C = A * B.A;
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(C,B.B,eps);
    add(rk);
  } else {
    DenseMatrix<R,T> D;
    D.copyFrom(B.B,true);
    init(C,D,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= rdr\n";
  DenseMatrix<R,T> C = A * B.A;
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(C,B.B,eps);
    subtract(rk);
  } else {
    DenseMatrix<R,T> D = -B.B;
    init(C,D,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= rdh\n";
  HMatrix<R,T> C(A,NULL,&(B.r));
  addproduct(C,B);
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= rdh\n";
  HMatrix<R,T> C(A,NULL,&(B.r));
  subtractproduct(C,B);
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= rrd\n";
  DenseMatrix<R,T> C = A.B * B;
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(A.A,C,eps);
    add(rk);
  } else {
    DenseMatrix<R,T> D;
    D.copyFrom(A.A,true);
    init(D,C,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= rrd\n";
  DenseMatrix<R,T> C = A.B * B;
  if(this->A.p && this->B.p) {
    Rk<R,T> rk(A.A,C,eps);
    subtract(rk);
  } else {
    DenseMatrix<R,T> D = -A.A;
    init(D,C,eps);
  }
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= rrh\n";
  DenseMatrix<R,T> C = A.B * B;
  addproduct(A.A,C);
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= rrh\n";
  DenseMatrix<R,T> C = A.B * B;
  subtractproduct(A.A,C);
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= rhd\n";
  HMatrix<R,T> C(B,&(A.c),NULL);
  addproduct(A,C);
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= rhd\n";
  HMatrix<R,T> C(B,&(A.c),NULL);
  subtractproduct(A,C);
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= rhr\n";
  DenseMatrix<R,T> C = A * B.A;
  addproduct(C,B.B);
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= rhr\n";
  DenseMatrix<R,T> C = A * B.A;
  subtractproduct(C,B.B);
}

template<class R, class T>
  void Rk<R,T> :: addproduct(const HMatrix<R,T> &C, const HMatrix<R,T> &D)
{
  if(DEBUG) cerr << "+= rhh\n";
  int m = C.getNumRows();
  int n = D.getNumCols();

  int rr = 0;
  int cc = 0;
  for(int i=0; i<C.m; i++) {
    cc = 0;
    for(int j=0; j<D.n; j++) {
      Rk<R,T> rk2;
      for(int k=0; k<C.n; k++) {
        rk2 += C(i,k) * D(k,j);
      }

      Rk<R,T> rk1(*this);
      
      int k = rk1.k + rk2.k;
      init(m,n,k,eps);

      A.zero();
      B.zero();

      T *cp = A.p;
      T *ap = rk1.A.p;
      T *bp = rk2.A.p;
      
      for(int q=0; q<rk1.k; q++) {
        memcpy(cp,ap,m*sizeof(T));
        cp += m;
        ap += rk1.A.lda;
      }
      cp += rr;
      for(int q=0; q<rk2.k; q++) {
        memcpy(cp,bp,rk2.m*sizeof(T));
        cp += m;
        bp += rk2.A.lda;
      }
  
      cp = B.p;
      ap = rk1.B.p;
      bp = rk2.B.p;

      if(rk1.k) {
        for(int q=0; q<n; q++) {
          memcpy(cp,ap,rk1.k*sizeof(T));
          cp += k;
          ap += rk1.B.lda;
        }
      }
      cp = B.p + (rk1.k + k*cc);
      for(int q=0; q<rk2.n; q++) {
        memcpy(cp,bp,rk2.k*sizeof(T));
        cp += k;
        bp += rk2.B.lda;
      }
    
      truncate();
    
      cc += D.c[j];
    }
    rr += C.r[i];
  }
}

template<class R, class T>
  void Rk<R,T> :: subtractproduct(const HMatrix<R,T> &C, const HMatrix<R,T> &D)
{
  if(DEBUG) cerr << "-= rhh\n";
  int m = C.getNumRows();
  int n = D.getNumCols();

  int rr = 0;
  int cc = 0;
  for(int i=0; i<C.m; i++) {
    cc = 0;
    for(int j=0; j<D.n; j++) {
      Rk<R,T> rk2;
      for(int k=0; k<C.n; k++) {
        rk2 -= C(i,k) * D(k,j);
      }

      Rk<R,T> rk1(*this);
      
      int k = rk1.k + rk2.k;
      init(m,n,k,eps);

      A.zero();
      B.zero();

      T *cp = A.p;
      T *ap = rk1.A.p;
      T *bp = rk2.A.p;
      
      for(int q=0; q<rk1.k; q++) {
        memcpy(cp,ap,m*sizeof(T));
        cp += m;
        ap += rk1.A.lda;
      }
      cp += rr;
      for(int q=0; q<rk2.k; q++) {
        memcpy(cp,bp,rk2.m*sizeof(T));
        cp += m;
        bp += rk2.A.lda;
      }
  
      cp = B.p;
      ap = rk1.B.p;
      bp = rk2.B.p;
      
      if(rk1.k) {
        for(int q=0; q<n; q++) {
          memcpy(cp,ap,rk1.k*sizeof(T));
          cp += k;
          ap += rk1.B.lda;
        }
      }
      cp = B.p + (rk1.k + k*cc);
      for(int q=0; q<rk2.n; q++) {
        memcpy(cp,bp,rk2.k*sizeof(T));
        cp += k;
        bp += rk2.B.lda;
      }
    
      truncate();
    
      cc += D.c[j];
    }
    rr += C.r[i];
  }
}

// L*X = B
template<class R, class T>
  void Rk<R,T> :: forwardSubstL(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsl rd\n";
  this->A.forwardSubstL(A);
}

template<class R, class T>
  void Rk<R,T> :: forwardSubstL(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve L*X=B for singular (Rk) L\n");
  abort();
}

template<class R, class T>
  void Rk<R,T> :: forwardSubstL(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsl rh\n";
  this->A.forwardSubstL(A);
}

// X*U = B
template<class R, class T>
  void Rk<R,T> :: forwardSubstU(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsu rd\n";
  this->B.forwardSubstU(A);
}

template<class R, class T>
  void Rk<R,T> :: forwardSubstU(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve X*U=B for singular (Rk) U\n");
  abort();
}

template<class R, class T>
  void Rk<R,T> :: forwardSubstU(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsu rh\n";
  this->B.forwardSubstU(A);
}

// U*X = B
template<class R, class T>
  void Rk<R,T> :: backwardSubstU(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "bsu rd\n";
  this->A.backwardSubstU(A);
}

template<class R, class T>
  void Rk<R,T> :: backwardSubstU(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve U*X=B for singular (Rk) U\n");
  abort();
}

template<class R, class T>
  void Rk<R,T> :: backwardSubstU(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "bsu rh\n";
  this->A.backwardSubstU(A);
}

// HMatrix LHS
template<class R, class T>
void HMatrix<R,T> :: copyFrom(const Rk<R,T> &A, bool bDeepCopy)
{
  fprintf(stderr,"Attempting to copy unformatted Rk matrix to H-matrix format\n");
  abort();
}

template<class R, class T>
void HMatrix<R,T> :: copyFrom(const DenseMatrix<R,T> &A, bool bDeepCopy)
{
  fprintf(stderr,"Attempting to copy unformatted dense matrix to H-matrix format\n");
  abort();
}

// binary ops
template<class R, class T>
  void HMatrix<R,T> :: add(const DenseMatrix<R,T> &B)
{
  HMatrix<R,T> C(B,&(this->r),&(this->c));
  add(C);
}

template<class R, class T>
  void HMatrix<R,T> :: add(const Rk<R,T> &B)
{
  HMatrix<R,T> C(B,&(this->r),&(this->c));
  add(C);
}

template<class R, class T>
  void HMatrix<R,T> :: add(const HMatrix<R,T> &B)
{
  int m = B.m;
  int n = B.n;
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      (*this)(i,j).add(B(i,j));
    }
  }
}

template<class R, class T>
  void HMatrix<R,T> :: subtract(const DenseMatrix<R,T> &B)
{
  HMatrix<R,T> C(B,&(this->r),&(this->c));
  subtract(C);
}

template<class R, class T>
  void HMatrix<R,T> :: subtract(const Rk<R,T> &B)
{
  HMatrix<R,T> C(B,&(this->r),&(this->c));
  subtract(C);
}

template<class R, class T>
  void HMatrix<R,T> :: subtract(const HMatrix<R,T> &B)
{
  int m = B.m;
  int n = B.n;
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      (*this)(i,j).subtract(B(i,j));
    }
  }
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= hdd\n";
  DenseMatrix<R,T> C = A * B;
  add(C);
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= hdd\n";
  DenseMatrix<R,T> C = A * B;
  subtract(C);
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= hdr\n";
  DenseMatrix<R,T> C = A * B.A;
  addproduct(C,B.B);      
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= hdr\n";
  DenseMatrix<R,T> C = A * B.A;
  subtractproduct(C,B.B);      
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= hdh\n";
  HMatrix<R,T> C(A,&(this->r),&(B.r));
  addproduct(C,B);      
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= hdh\n";
  HMatrix<R,T> C(A,&(this->r),&(B.r));
  subtractproduct(C,B);      
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= hrd\n";
  DenseMatrix<R,T> C = A.B * B;
  addproduct(A.A,C);      
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= hrd\n";
  DenseMatrix<R,T> C = A.B * B;
  subtractproduct(A.A,C);      
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const Rk<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= hrr\n";
  DenseMatrix<R,T> C = A.B * B.A;
  DenseMatrix<R,T> D = A.A * C;
  addproduct(D,B.B);      
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const Rk<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= hrr\n";
  DenseMatrix<R,T> C = A.B * B.A;
  DenseMatrix<R,T> D = A.A * C;
  subtractproduct(D,B.B);      
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= hrh\n";
  DenseMatrix<R,T> C = A.B * B;
  addproduct(A.A,C);
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= hrh\n";
  DenseMatrix<R,T> C = A.B * B;
  subtractproduct(A.A,C);
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= hhd\n";
  HMatrix<R,T> C(B,&(A.c),&(this->c));
  addproduct(A,C);      
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= hhd\n";
  HMatrix<R,T> C(B,&(A.c),&(this->c));
  subtractproduct(A,C);      
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "+= hhr\n";
  DenseMatrix<R,T> C = A * B.A;
  addproduct(C,B.B);
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)
{
  if(DEBUG) cerr << "-= hhr\n";
  DenseMatrix<R,T> C = A * B.A;
  subtractproduct(C,B.B);
}

template<class R, class T>
  void HMatrix<R,T> :: addproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "+= hhh\n";
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      Matrix<R,T> &C = (*this)(i,j);
      for(int k=0; k<A.n; k++) {
        C += A(i,k) * B(k,j);
      }
    }
  }   
}

template<class R, class T>
  void HMatrix<R,T> :: subtractproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B)
{
  if(DEBUG) cerr << "-= hhh\n";
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      Matrix<R,T> &C = (*this)(i,j);
      for(int k=0; k<A.n; k++) {
        C -= A(i,k) * B(k,j);
      }
    }
  }   
}

template<class R, class T>
  void HMatrix<R,T> :: forwardSubstL(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsl hd " << m << " " << n;
  if(m != 1)
    abort();
  for(int j=0; j<n; j++) {
    Matrix<R,T> &C = (*this)(0,j);
    C.forwardSubstL(A);
  }
}

template<class R, class T>
  void HMatrix<R,T> :: forwardSubstL(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve L*X=B for singular (Rk) L\n");
  abort();  
}

template<class R, class T>
  void HMatrix<R,T> :: forwardSubstL(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsl hh " << getNumRows() << " " << getNumCols() << " " << A.getNumRows() << " " << A.getNumCols() << "\n";
  HMatrix<R,T> &B = *this;

  for(int j=0; j<n; j++) {
    for(int i=0; i<m; i++) {
      Matrix<R,T> &C = B(i,j);
      for(int k=0; k<i; k++) {
        C -= A(i,k) * B(k,j);
      }
      C.forwardSubstL(A(i,i));
    }
  }
}

template<class R, class T>
  void HMatrix<R,T> :: forwardSubstU(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsu hd " << m << " " << n;
  if(n != 1)
    abort();
  for(int i=0; i<m; i++) {
    Matrix<R,T> &C = (*this)(i,0);
    C.forwardSubstU(A);
  }
}

template<class R, class T>
  void HMatrix<R,T> :: forwardSubstU(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve X*U=B for singular (Rk) U\n");
  abort();
}

template<class R, class T>
  void HMatrix<R,T> :: forwardSubstU(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "fsu hh\n";
  HMatrix<R,T> &B = *this;
  
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      Matrix<R,T> &C = B(i,j);
      for(int k=0; k<j; k++) {
        C -= B(i,k) * A(k,j);
      }
      C.forwardSubstU(A(j,j));
    }
  }           
}

template<class R, class T>
  void HMatrix<R,T> :: backwardSubstU(const DenseMatrix<R,T> &A)
{
  if(DEBUG) cerr << "bsu hd " << m << " " << n;
  if(m != 1)
    abort();
  for(int j=0; j<n; j++) {
    Matrix<R,T> &C = (*this)(0,j);
    C.backwardSubstU(A);
  }
}

template<class R, class T>
  void HMatrix<R,T> :: backwardSubstU(const Rk<R,T> &A)
{
  fprintf(stderr, "Attempting to solve U*X=B for singular (Rk) U\n");
  abort();  
}

template<class R, class T>
  void HMatrix<R,T> :: backwardSubstU(const HMatrix<R,T> &A)
{
  if(DEBUG) cerr << "bsu hh\n";
  HMatrix<R,T> &B = *this;

  for(int j=0; j<n; j++) {
    for(int i=m-1; i>=0; i--) {
      Matrix<R,T> &C = B(i,j);
      for(int k=i+1; k<m; k++) {
        C -= A(i,k) * B(k,j);
      }
      C.backwardSubstU(A(i,i));
    }
  }
}

#undef DEBUG

#endif
