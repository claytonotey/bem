#ifndef HMATRIX_H
#define HMATRIX_H

#include "ReferenceCountingPointer.h"
#include "MatrixBase.h"
#include <iostream>

#include <vector>
using namespace std;

template<class R, class T>
class DenseMatrix;

template<class R, class T>
class Rk;

#define DEBUG 0

template<class R, class T>
class HMatrix : public Matrix<R,T> {
 public:

  int getNumRows() const { 
    int mm = 0;
    for(int i=0; i<m; i++)
      mm += r[i];
    return mm;
  }

  int getNumCols() const {     
    int nn = 0;
    for(int j=0; j<n; j++)
      nn += c[j];
    return nn;
  }

  void init(int m, int n) {
    abort();
    this->m = 1;
    this->n = 1;
    sub = ReferenceCountingPointer< Matrix<R,T>* >(1);
    r.resize(1);
    c.resize(1);
    r[0] = m;
    c[0] = n;
  }

  Matrix<R,T> &operator()(int i, int j) { return *(sub[i+m*j]); }
  Matrix<R,T> &operator()(int i, int j) const { return *(sub[i+m*j]); }
  void set(int i, int j, Matrix<R,T>* A) { sub[i+m*j] = A; }

  HMatrix() : Matrix<R,T>(Matrix<R,T>::h), m(0), n(0), sub(), r(0), c(0) {}
  HMatrix(int m_, int n_) : Matrix<R,T>(Matrix<R,T>::h), m(m_), n(n_), sub(m_*n_), r(m_), c(n_) {}
  HMatrix(const DenseMatrix<R,T> &A, 
          const Array<int> *r_,
          const Array<int> *c_) 
    : Matrix<R,T>(Matrix<R,T>::h),
      m(r_?r_->n:1), n(c_?c_->n:1),
      sub((r_?r_->n:1)*(c_?c_->n:1)), 
      r(r_?r_->n:1), c(c_?c_->n:1) {
    
    if(r_) {
      for(int i=0; i<m; i++) {
        r[i] = (*r_)[i];
      }
    } else {
      r[0] = A.m;
    }
    
    if(c_) {
      for(int j=0; j<n; j++) {
        c[j] = (*c_)[j];
      }
    } else {
      c[0] = A.n;
    }
    
    int rr = 0;
    int cc = 0;
    for(int i=0; i<m; i++) {
      cc = 0;
      for(int j=0; j<n; j++) {
        DenseMatrix<R,T> *B = new DenseMatrix<R,T>(A,rr,cc,r[i],c[j]);
        set(i,j,B);
        cc += c[j];
      }
      rr += r[i];
    }
  }

  HMatrix(const Rk<R,T> &A, 
          const Array<int> *r_,
          const Array<int> *c_) 
    : Matrix<R,T>(Matrix<R,T>::h),
      m(r_?r_->n:1), n(c_?c_->n:1),
      sub((r_?r_->n:1)*(c_?c_->n:1)), 
      r(r_?r_->n:1), c(c_?c_->n:1) {
    
    if(r_) {
      for(int i=0; i<m; i++)
        r[i] = (*r_)[i];
    } else {
      r[0] = A.m;
    }
    
    if(c_) {
      for(int j=0; j<n; j++) {
        c[j] = (*c_)[j];
      }
    } else {
      c[0] = A.n;
    }
    
    int rr = 0;
    int cc = 0;
    for(int i=0; i<m; i++) {
      cc = 0;
      for(int j=0; j<n; j++) {
        DenseMatrix<R,T> CA(A.A,rr,0,r[i],A.k);
        DenseMatrix<R,T> CB(A.B,0,cc,A.k,c[j]);
        Rk<R,T> *B = new Rk<R,T>(CA,CB,A.eps);
        set(i,j,B);
        cc += c[j];
      }
      rr += r[i];
    }
  }

  ~HMatrix()
  {
    deleteSubs();
  }

  void deleteSubs()
  {
    int size = m * n;
    for(int i=0; i<size; i++)
      delete sub[i];
  }

  void zero() 
  {
    int size = m * n;
    for(int i=0; i<size; i++)
      sub[i]->zero();
  }

  void copyFrom(const DenseMatrix<R,T> &A, bool bDeepCopy = false);
  void copyFrom(const Rk<R,T> &A, bool bDeepCopy = false);
  void copyFrom(const HMatrix<R,T> &A, bool bDeepCopy = false)
  {
    m = A.m;
    n = A.n;    
    if(bDeepCopy) {
      abort();
    } else {
      sub = A.sub;
    }
    r.resize(m);
    c.resize(n);
    for(int i=0; i<m; i++)
      r[i] = A.r[i];
    for(int j=0; j<n; j++)
      c[j] = A.c[j];
  }

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
    os << "]\n";
  }

  void writeBinary(fstream &fs) {
    fs.write((char*)&m,sizeof(int));
    fs.write((char*)&n,sizeof(int));
    for(int i=0; i<m; i++) {
      int x = r[i];
      fs.write((char*)&x,sizeof(int));
    }
    for(int j=0; j<n; j++) {
      int x = c[j];
      fs.write((char*)&x,sizeof(int));
    }
    for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
        Matrix<R,T> *A = sub[i+m*j];
        int type = A->type;
        fs.write((char*)&type,sizeof(int));
        A->writeBinary(fs);
      }
    }
  }

  void readBinary(fstream &fs) {
    fs.read((char*)&m,sizeof(int));
    fs.read((char*)&n,sizeof(int));
    sub = ReferenceCountingPointer< Matrix<R,T>* >(m*n);
    r.resize(m);
    for(int i=0; i<m; i++) {
      int x;
      fs.read((char*)&x,sizeof(int));
      r[i] = x;
    }
    c.resize(n);
    for(int j=0; j<n; j++) {
      int x;
      fs.read((char*)&x,sizeof(int));
      c[j] = x;
    }
    for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
        int type;
        fs.read((char*)&type,sizeof(int));
        Matrix<R,T> *A;
        if(type==Matrix<R,T>::dense)
          A = new DenseMatrix<R,T>;
        else if(type==Matrix<R,T>::rk)
          A = new Rk<R,T>;
        else if(type==Matrix<R,T>::h)
          A = new HMatrix<R,T>;
        else
          abort();

        A->readBinary(fs);
        sub[i+m*j] = A;
      }
    }
  }

  void LU()
  {
    if(DEBUG) cerr << "lu h " << getNumRows() << " " << getNumCols() << "\n";

    HMatrix<R,T> &B = *this;
    for(int i=0; i<m; i++) {
      
      // diagonal
      Matrix<R,T> &C = B(i,i);
      for(int k=0; k<i; k++) {
        C -= B(i,k) * B(k,i);
      }
      C.LU();
      
      // lower
      for(int j=i+1; j<m; j++) {
        Matrix<R,T> &E = B(j,i);
        for(int k=0; k<i; k++) {
          E -= B(j,k) * B(k,i);
        }
        E.forwardSubstU(C);
      }
      
      // upper
      for(int j=i+1; j<n; j++) {
        Matrix<R,T> &E = B(i,j);
        for(int k=0; k<i; k++) {
          E -= B(i,k) * B(k,j);
        }
        E.forwardSubstL(C);
      }
    }
  }

  template<class T1, class T2>
    void operator=(const MatrixMultiply<R,T1,T2> &op)
  {
    zero();
    addproduct(op.A,op.B);
  }

  int m;
  int n;
  ReferenceCountingPointer< Matrix<R,T>* > sub;
  Array<int> r;
  Array<int> c;

};


template<class R, class T>
void solveLU(const HMatrix<R,T> &A, HMatrix<R,T> &B)
{
  int m = B.m;
  int n = B.n;

  if(A.m != A.n || A.n != m)
    abort();
  
  if(DEBUG) cerr << "slu hh " << B.getNumRows() << " " << B.getNumCols() << " " << A.getNumRows() << " " << A.getNumCols() << "\n";

  // Solve LY=B
  for(int j=0; j<n; j++) {
    for(int i=0; i<m; i++) {
      Matrix<R,T> &C = B(i,j);
      for(int k=0; k<i; k++) {
        C -= A(i,k) * B(k,j);
      }
      C.forwardSubstL(A(i,i));
    }
  }

  // Solve UX=Y
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

typedef HMatrix<Real,Complex> ComplexHMatrix;

#undef DEBUG

#endif
