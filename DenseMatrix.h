#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include "MatrixBase.h"
#include "ReferenceCountingPointer.h"
#include <iostream>
#include <fstream>
#include "Array.h"
#include <vector>
using namespace std;

template<class R, class T>
class Rk;

template<class R, class T>
class HMatrix;

template<class R, class T>
class DenseMatrix : public Matrix<R,T> {
 public:

  enum colon { allrows = -999, allcols };

  DenseMatrix() : Matrix<R,T>(Matrix<R,T>::dense), m(0), n(0), lda(0), p(), piv() {}
  
  DenseMatrix(int m_, int n_, char s) : Matrix<R,T>(Matrix<R,T>::dense), m(m_), n(n_), lda(m_), p(m_*n_), piv()
  {
    memset(p,s,m*n*sizeof(T));
  }

  DenseMatrix(int m_, int n_) : Matrix<R,T>(Matrix<R,T>::dense), m(m_), n(n_), lda(m_), p(m_*n_), piv() {
  }

  DenseMatrix(const Matrix<R,T> &A, bool bDeepCopy = false) : Matrix<R,T>(Matrix<R,T>::dense)
  {
    Matrix<R,T>::copyFrom(A,bDeepCopy);
  }

  void copyFrom(const DenseMatrix &A, bool bDeepCopy = false)
  {
    m = A.m;
    n = A.n;
    if(bDeepCopy) {
      lda = m;
      p = ReferenceCountingPointer<T>(m*n);
      for(int j=0; j<n; j++) {
        memcpy(p+j*lda,A.p+j*A.lda,m*sizeof(T));
      }
      if(A.piv) {
        piv = ReferenceCountingPointer<int>(m);
        memcpy(piv,A.piv,m*sizeof(int));
      }
    } else {
      lda = A.lda;
      p = A.p;
      if(A.piv) {
        piv = A.piv;
      }
    }
  }
  void copyFrom(const Rk<R,T> &A, bool bDeepCopy = false);
  void copyFrom(const HMatrix<R,T> &A, bool bDeepCopy = false);

  DenseMatrix(const DenseMatrix<R,T> &A, int i, int j, int m_, int n_) : Matrix<R,T>(Matrix<R,T>::dense), m(m_), n(n_), lda(A.lda), p(A.p,i+j*A.lda) {}

  void init(int m_, int n_)
  {
    m = m_;
    n = n_;
    lda = m_;
    p = ReferenceCountingPointer<T>(m*n);
  }

  ~DenseMatrix() {}

  int getNumRows() const { return m; }
  int getNumCols() const { return n; }

  void zero() 
  {
    if(p) {
      for(int j=0; j<n; j++) {
        memset(p+j*lda,0,m*sizeof(T));
      }
    }
  }

  void resize(int is, int js, int m, int n) 
  {
    if(p) {
      p = ReferenceCountingPointer<T>(p,is+js*lda);
    } else {
      p = ReferenceCountingPointer<T>(m*n);
      lda = m;
    }
    this->m = m;
    this->n = n;
  }

  void resize(int m, int n) 
  {
    resize(0,0,m,n);
  }

  DenseMatrix operator()(vector<int> &r, vector<int> c) const
  {
    int mm = (int) r.size();
    int nn = (int) c.size();
    DenseMatrix C(mm,nn);
    
    for(int i=0; i<mm; i++) {
      int ri = r[i];
      for(int j=0; j<nn; j++) {
        C(i,j) = p[ri + lda * c[j]];
      }
    }

    if(mm==nn && piv) {
      C.piv = ReferenceCountingPointer<int>(mm);
      for(int i=0; i<mm; i++) {
        int ri = r[i];
        if(ri != c[i])
          abort();
        C.piv[i] = piv[ri];
      }
    }

    return C;
  }

  DenseMatrix operator()(int r, vector<int> c) const
  {
    int mm = (r == allrows) ? m : 1;
    int nn = (int) c.size();
    DenseMatrix C(mm,nn);
    
    for(int i=0; i<mm; i++) {
      int ii = (r == allrows) ? i : r;
      for(int j=0; j<nn; j++) {
        C(i,j) = p[ii + lda * c[j]];
      }
    }

    if(mm==nn && piv.p) {
      C.piv = ReferenceCountingPointer<int>(mm);
      for(int i=0; i<mm; i++) {
        int ri = c[i];
        C.piv[i] = piv[ri];
      }
    }

    return C;
  }

  DenseMatrix operator()(vector<int> r, int c) const
  {
    int mm = (int) r.size();
    int nn = (c == allcols) ? n : 1;
    DenseMatrix C(mm,nn);
    
    for(int i=0; i<mm; i++) {
      for(int j=0; j<nn; j++) {
        int jj = (c == allcols) ? j : c;      
        C(i,j) = p[r[i] + lda * jj];
      }
    }

    if(mm==nn && piv.p) {
      C.piv = ReferenceCountingPointer<int>(mm);
      for(int i=0; i<mm; i++) {
        int ri = r[i];
        C.piv[i] = piv[ri];
      }
    }

    return C;
  }

  void negate() {
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)     
        (*this)(i,j) = -(*this)(i,j);
  }

  DenseMatrix operator-() const
  {
    DenseMatrix A(m,n);
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++)     
        A(i,j) = -(*this)(i,j);
    return A;
  }

  inline T& operator()(int i, int j)
  {
    return p[i + lda * j];
  }

  inline T& operator()(int i, int j) const
  {
    return p[i + lda * j];
  }

  inline void writeColumn(T *dest, int j) {
    memcpy(dest, p + lda * j, m * sizeof(T));
  }

  inline void readColumn(T *src, int j) {
    memcpy(p + lda * j, src, m * sizeof(T));
  }

  friend void aca(DenseMatrix<R,T> &A, vector<int> &r, vector<int> &c, R eps2)
  {
    int nu = 0;
    R epsrel2 = 0.;
    R epsabs2 = 0.;
    R epsabstot2 = 0.;
    int m = A.m;
    int n = A.n;
    DenseMatrix<R,T> M(A,true); 
    Array<T> a(m);
    Array<T> b(n);
    int k = min(m,n);
    bool bCont;

    do {
      R maxM2 = -1.;
      int is;
      int js;
      for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
          R M2 = norm(M(i,j));
          if(M2 > maxM2) {
            maxM2 = M2;
            is = i;
            js = j;
          }
        }
      }
      
      T delta = M(is,js);
      if(norm(delta) == 0.)
        break;
      
      R norma = 0.;
      for(int i=0; i<m; i++) {
        T ai = M(i,js);
        a[i] = ai;
        norma += norm(ai);
      }
      
      R normb = 0.;
      for(int j=0; j<n; j++) {
        T bj = M(is,j) / delta;
        b[j] = bj;
        normb += norm(bj);
      }
      
      for(int i=0; i<m; i++)
        for(int j=0; j<n; j++)
          M(i,j) -= a[i] * b[j];
      
      epsabs2 = norma * normb;
      epsabstot2 = epsabstot2 + epsabs2;
      epsrel2 = epsabs2 / epsabstot2;
      
      nu++;

      bCont = (nu < k && epsrel2 > eps2);
      
      if(bCont) {
        r.push_back(is);
        c.push_back(js);
      }     
    } while(bCont);
    
    //printf("aca %d/%d\n",(int)r.size(),k);
  }
  
  friend R norm(DenseMatrix &A)
  {
    R n = 0.;
    for(int i=0; i<A.m; i++)
      for(int j=0; j<A.n; j++)
        n += norm(A(i,j));
    return n;
  }

  // binary ops
  void add(const DenseMatrix<R,T> &A);
  void add(const Rk<R,T> &A);
  void add(const HMatrix<R,T> &A);
  
  void subtract(const DenseMatrix<R,T> &A);
  void subtract(const Rk<R,T> &A);
  void subtract(const HMatrix<R,T> &A);
  
  void addproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B);
  void addproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B);
  void addproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B);
  void addproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B);
  void addproduct(const Rk<R,T> &A, const Rk<R,T> &B);
  void addproduct(const Rk<R,T> &A, const HMatrix<R,T> &B);
  void addproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B);
  void addproduct(const HMatrix<R,T> &A, const Rk<R,T> &B);
  void addproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B);

  void subtractproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B);
  void subtractproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B);
  void subtractproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B);
  void subtractproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B);
  void subtractproduct(const Rk<R,T> &A, const Rk<R,T> &B);
  void subtractproduct(const Rk<R,T> &A, const HMatrix<R,T> &B);
  void subtractproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B);
  void subtractproduct(const HMatrix<R,T> &A, const Rk<R,T> &B);
  void subtractproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B);

  void forwardSubstL(const DenseMatrix<R,T> &A);
  void forwardSubstL(const Rk<R,T> &A);
  void forwardSubstL(const HMatrix<R,T> &A);

  void forwardSubstU(const DenseMatrix<R,T> &A);
  void forwardSubstU(const Rk<R,T> &A);
  void forwardSubstU(const HMatrix<R,T> &A);

  void backwardSubstU(const DenseMatrix<R,T> &A);
  void backwardSubstU(const Rk<R,T> &A);
  void backwardSubstU(const HMatrix<R,T> &A);

  // binary op constructors
  DenseMatrix(const MatrixAdd<R,T> &op) : Matrix<R,T>(Matrix<R,T>::dense)
  {
    Matrix<R,T>::copyFrom(op.A,true);
    Matrix<R,T>::add(op.B);
  }

  DenseMatrix(const MatrixSubtract<R,T> &op) : Matrix<R,T>(Matrix<R,T>::dense)
  {
    Matrix<R,T>::copyFrom(op.A,true);
    Matrix<R,T>::subtract(op.B);
  }

  DenseMatrix(const MatrixMultiply<R,T,T> &op) : Matrix<R,T>(Matrix<R,T>::dense)
  {
    init(op.A.getNumRows(),op.B.getNumCols());
    zero();
    Matrix<R,T>::addproduct(op.A,op.B);
  }

  void output(ostream &os)
  {
    os << "[\n";
    for(int i = 0; i < m; i++) {
      os << "[ ";
      for(int j = 0; j < n; j++) {
        os << (*this)(i,j) << " ";
      }
      os << "]\n";
    }
    os << "]";
  }

  void writeBinary(fstream &fs)
  {
    int s = sizeof(T);
    fs.write((char*)&s,sizeof(int));
    fs.write((char*)&m,sizeof(int));
    fs.write((char*)&n,sizeof(int));
    for(int j=0; j<n; j++)
      fs.write((char*)(p+j*lda),m*s);
  }

  void readBinary(fstream &fs)
  {
    int s;
    fs.read((char*)&s,sizeof(int));
    fs.read((char*)&m,sizeof(int));
    fs.read((char*)&n,sizeof(int));
    lda = m;
    p = ReferenceCountingPointer<T>(m*n);
    for(int j=0; j<n; j++)
      fs.read((char*)(p+j*lda),m*s);
  }

  void LU() {
    fprintf(stderr,"Undefined operation: LU decomposition\n");
    abort();
  }

  template<class T1, class T2>
    void operator=(const MatrixMultiply<R,T1,T2> &op)
  {
    init(op.A.getNumRows(),op.B.getNumCols());
    zero();
    Matrix<R,T>::addproduct(op.A,op.B);
  }

 public:
  int m;
  int n;
  int lda;
  ReferenceCountingPointer<T> p;
  ReferenceCountingPointer<int> piv;
};

template<class R, class T1, class T2>
void solveLU(DenseMatrix<R,T1> &A, DenseMatrix<R,T2> &B)
{
  int m = B.m;
  int n = B.n;

  if(A.m != A.n || A.n != m)
    abort();

  int *piv = A.piv;
  
  // swap rows of B
  for(int i=0; i<m; i++) {
    int pi = piv[i]-1;
    for(int k=0; k<n; k++) {
      swap(B(i,k),B(pi,k));
    }
  }    

  // Solve LY=B
  for(int j=0; j<n; j++) {
    for(int i=0; i<m; i++) {
      T2 &c = B(i,j);
      for(int k=0; k<i; k++) {
        c -= A(i,k) * B(k,j);
      }
      B(i,j) = c;
    }
  }

  // Solve UX=Y
  for(int j=0; j<n; j++) {
    for(int i=m-1; i>=0; i--) {
      T2 &c = B(i,j);
      for(int k=i+1; k<m; k++) {
        c -= A(i,k) * B(k,j);
      }
      c *= (R)1. / A(i,i);
      B(i,j) = c;
    }
  }
}

template<class R, class T1, class T2>
void solve(DenseMatrix<R,T1> &A, DenseMatrix<R,T2> &B)
{
  int m = B.m;

  if(A.m != A.n || A.n != m)
    abort();

  A.LU();
  solveLU(A,B);
}

template<class R, class T>
DenseMatrix<R,T> diag(Array<T> &d) {
  return diag(d,d.n,d.n);
}

template<class R, class T>
DenseMatrix<R,T> diag(Array<T> &d, int m, int n) 
{
  DenseMatrix<R,T> D(m,n,0);
  for(int i=0; i<d.n; i++)
    D(i,i) = d[i];
  
  return D;
}

template<class R, class T>
  DenseMatrix<R,T> transpose(DenseMatrix<R,T> &A)
{
  DenseMatrix<R,T> B(A.n,A.m);
  for(int i=0; i<A.m; i++)
    for(int j=0; j<A.n; j++)
      B(j,i) = A(i,j);
  
  return B;
}


template<class R, class T>
  T trace(DenseMatrix<R,T> &A)
{
  T t = 0.0;
  for(int i=0; i<A.m; i++)
    t += A(i,i);
  
  return t;
}

typedef DenseMatrix<Real,Real> RealDenseMatrix;
typedef DenseMatrix<Real,Complex> ComplexDenseMatrix;

template<class R, class T>
  DenseMatrix<R,T> eye(int n) {
  fprintf(stderr,"Undefined matrix operation: eye\n");
  abort();
}


template<>
ComplexDenseMatrix eye(int n);

template<>
void ComplexDenseMatrix :: LU();

template<>
void ComplexDenseMatrix :: output(ostream &os);

template<>
void solveLU(ComplexDenseMatrix &A, ComplexDenseMatrix &B);

Real condition(ComplexDenseMatrix &A);
void eig(ComplexDenseMatrix &A, Array<Complex> &w);
void lu(DenseMatrix<double,ComplexDouble> &A, int *ipiv);
void inverse(DenseMatrix<double,ComplexDouble> &A);
void svd(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &U, Array<double> &S, DenseMatrix<double,ComplexDouble> &V);
void svd2(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &U, Array<double> &S, DenseMatrix<double,ComplexDouble> &V);
void qr(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &Q, DenseMatrix<double,ComplexDouble> &R);
void rq(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &R, DenseMatrix<double,ComplexDouble> &Q);
void qr(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau, DenseMatrix<double,ComplexDouble> &R);
void lmq(DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau, DenseMatrix<double,ComplexDouble> &C);
void rq(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &R, DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau);
void rmq(DenseMatrix<double,ComplexDouble> &C, DenseMatrix<double,ComplexDouble> &Q, Array<ComplexDouble> &tau);
void multiply_A_Bdagger(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &B, DenseMatrix<double,ComplexDouble> &C);
void multiply_Adagger_B(DenseMatrix<double,ComplexDouble> &A, DenseMatrix<double,ComplexDouble> &B, DenseMatrix<double,ComplexDouble> &C);
void symm(DenseMatrix<double,ComplexDouble> &A);

#endif
