#ifndef RK_H
#define RK_H

#include "MatrixBase.h"
#include "DenseMatrix.h"

#define RK_DEFAULTEPS 0.0

template<class R, class T>
class HMatrix;

template<class R, class T>
class Rk : public Matrix<R,T> {
 public:
  
  Rk() : Matrix<R,T>(Matrix<R,T>::rk), A(), B(), eps(RK_DEFAULTEPS), m(0), n(0), k(0) {}

  Rk(const DenseMatrix<R,T> &A_, const DenseMatrix<R,T> &B_, R eps_, bool bDeepCopy = false) 
    : Matrix<R,T>(Matrix<R,T>::rk), A(A_,bDeepCopy), B(B_,bDeepCopy), eps(eps_), m(A.m), n(B.n), k(A.n) {}

  Rk(int m_, int n_, int k_, R eps_) : Matrix<R,T>(Matrix<R,T>::rk), A(m_,k_), B(k_,n_), eps(eps_), m(m_), n(n_), k(k_) {}

  Rk(const Matrix<R,T> &A, bool bDeepCopy = false) : Matrix<R,T>(Matrix<R,T>::rk)
  {
    Matrix<R,T>::copyFrom(A,bDeepCopy);
  }

  void init(const DenseMatrix<R,T> &A_, const DenseMatrix<R,T> &B_, R eps_, bool bDeepCopy = false)
  { 
    A.copyFrom(A_,bDeepCopy);
    B.copyFrom(B_,bDeepCopy);
    eps = eps_;
    m = A.m;
    n = B.n;
    k = A.n;
  }

  void init(int m_, int n_, int k_, R eps_) {
    A.init(m_,k_);
    B.init(k_,n_);
    eps = eps_;
    m = m_;
    n = n_;
    k = k_;
  }

  void init(int m, int n) {
    init(m,n,0,eps);
  }

  ~Rk() {}

  int getNumRows() const { return m; }
  int getNumCols() const { return n; }

  void zero() {
    if(A.p && B.p) {
      A.zero();
      B.zero();
    }
  }

  void copyFrom(const DenseMatrix<R,T> &A, bool bDeepCopy = false);
  void copyFrom(const Rk<R,T> &A, bool bDeepCopy = false)
  {
    init(A.A,A.B,A.eps,bDeepCopy);
  }
  void copyFrom(const HMatrix<R,T> &A, bool bDeepCopy = false);
  
  void truncate();

  DenseMatrix<R,T> toFull()
  {
    DenseMatrix<R,T> C = A * B;
    return C;
  }

  void output(ostream &os)
  {
    os << "[A = \n" << A << "\n";
    os << "B = \n" << B << "]\n";
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


  // binary op constructors
  Rk(const MatrixAdd<R,T> &op) : Matrix<R,T>(Matrix<R,T>::rk)
  {
    Matrix<R,T>::copyFrom(op.A,true);
    Matrix<R,T>::add(op.B);
  }

  Rk(const MatrixSubtract<R,T> &op) : Matrix<R,T>(Matrix<R,T>::rk)
  {
    Matrix<R,T>::copyFrom(op.A,true);
    Matrix<R,T>::subtract(op.B);
  }

  Rk(const MatrixMultiply<R,T,T> &op) : Matrix<R,T>(Matrix<R,T>::rk)
  {
    zero();
    Matrix<R,T>::addproduct(op.A,op.B);
  }

  void LU() {
    fprintf(stderr,"Attempting to perform LU decomposition of singular (Rk) matrix\n");
    abort();
  }

  void writeBinary(fstream &fs) {
    fs.write((char*)&eps,sizeof(R));
    fs.write((char*)&m,sizeof(int));
    fs.write((char*)&n,sizeof(int));
    fs.write((char*)&k,sizeof(int));
    A.writeBinary(fs);
    B.writeBinary(fs);
  }

  void readBinary(fstream &fs) {
    fs.read((char*)&eps,sizeof(R));
    fs.read((char*)&m,sizeof(int));
    fs.read((char*)&n,sizeof(int));
    fs.read((char*)&k,sizeof(int));
    A.readBinary(fs);
    B.readBinary(fs);
  }

  template<class T1, class T2>
    void operator=(const MatrixMultiply<R,T1,T2> &op)
  {
    zero();
    addproduct(op.A,op.B);
  }

  DenseMatrix<R,T> A;
  DenseMatrix<R,T> B;
  R eps;
  int m;
  int n;
  int k;
};

typedef Rk<Real, Complex> ComplexRk;

template<>
void ComplexRk :: truncate();

#endif
