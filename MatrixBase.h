#ifndef MATRIX_H
#define MATRIX_H

#include "defs.h"

#include <ostream>
#include <fstream>
using namespace std;

// Binary Matrix Operations
template<class R, class T>
class Matrix;

template<class R, class T1, class T2>
class BinaryMatrixOperation {
 public:
  BinaryMatrixOperation(const Matrix<R,T1> &A_, const Matrix<R,T2> &B_) : A(A_), B(B_) {}
  const Matrix<R,T1> &A;
  const Matrix<R,T2> &B;
};

template<class R, class T>
  class MatrixAdd : public BinaryMatrixOperation<R,T,T> {
 public:   
 MatrixAdd(const Matrix<R,T> &A_, const Matrix<R,T> &B_) : BinaryMatrixOperation<R,T,T>(A_,B_) {}
};

template<class R, class T>
  class MatrixSubtract : public BinaryMatrixOperation<R,T,T> {
 public:   
 MatrixSubtract(const Matrix<R,T> &A_, const Matrix<R,T> &B_) : BinaryMatrixOperation<R,T,T>(A_,B_) {}
};

template<class R, class T1, class T2>
  class MatrixMultiply : public BinaryMatrixOperation<R,T1,T2> {
 public:
 MatrixMultiply(const Matrix<R,T1> &A_, const Matrix<R,T2> &B_) : BinaryMatrixOperation<R,T1,T2>(A_,B_) {}
};

template<class R, class T>
inline MatrixAdd<R,T> operator+(const Matrix<R,T> &A, const Matrix<R,T> &B)
{
  return MatrixAdd<R,T>(A,B);
}

template<class R, class T>
inline MatrixSubtract<R,T> operator-(const Matrix<R,T> &A, const Matrix<R,T> &B)
{
  return MatrixSubtract<R,T>(A,B);
}

template<class R, class T1, class T2>
inline MatrixMultiply<R,T1,T2> operator*(const Matrix<R,T1> &A, const Matrix<R,T2> &B)
{
  return MatrixMultiply<R,T1,T2>(A,B);
}

// Matrix Base
template<class R, class T>
class DenseMatrix;

template<class R, class T>
class Rk;

template<class R, class T>
class HMatrix;

template<class R, class T>
class Matrix {
 public:
  enum Type { dense, rk, h };
  int type;
  Matrix(int type_) : type(type_) {}
  virtual ~Matrix() {}
  virtual void zero()=0;
  virtual void init(int m, int n)=0;
  virtual void writeBinary(fstream &fs)=0;
  virtual void readBinary(fstream &fs)=0;
  virtual void output(ostream &os)=0;
  friend ostream& operator<<(ostream& os, Matrix &A) {
    A.output(os);
    return os;
  }
  virtual int getNumRows() const=0;
  virtual int getNumCols() const=0;
  // operations
  virtual void copyFrom(const DenseMatrix<R,T> &A, bool bDeepCopy = false)=0;
  virtual void copyFrom(const Rk<R,T> &A, bool bDeepCopy = false)=0;
  virtual void copyFrom(const HMatrix<R,T> &A, bool bDeepCopy = false)=0;
  inline void copyFrom(const Matrix &A, bool bDeepCopy = false)
  {
    switch(A.type) {
    case Matrix<R,T>::dense:
      copyFrom((const DenseMatrix<R,T>&)A,bDeepCopy);
      break;
    case Matrix<R,T>::rk:
      copyFrom((const Rk<R,T>&)A,bDeepCopy);
      break;
    case Matrix<R,T>::h:
      copyFrom((const HMatrix<R,T>&)A,bDeepCopy);
      break;
    default:
      abort();
    }
  }

  inline void operator=(const Matrix &A)
  {
    copyFrom(A);
  }

  // add
  virtual void add(const DenseMatrix<R,T> &A)=0;
  virtual void add(const Rk<R,T> &A)=0;
  virtual void add(const HMatrix<R,T> &A)=0;
  inline void add(const Matrix<R,T> &A)
  {
    switch(A.type) {
    case Matrix<R,T>::dense:
      add((const DenseMatrix<R,T>&)A);
      break;
    case Matrix<R,T>::rk:
      add((const Rk<R,T>&)A);
      break;
    case Matrix<R,T>::h:
      add((const HMatrix<R,T>&)A);
      break;
    default:
      abort();
    }
  }
  inline void operator+=(const Matrix &A)
  {
    add(A);
  }
  inline void operator=(MatrixAdd<R,T> op)
  {
    this->copyFrom(op.A,true);
    add(op.B);
  }

  // subtract
  virtual void subtract(const DenseMatrix<R,T> &A)=0;
  virtual void subtract(const Rk<R,T> &A)=0;
  virtual void subtract(const HMatrix<R,T> &A)=0;
  inline void subtract(const Matrix<R,T> &A)
  {
    switch(A.type) {
    case Matrix<R,T>::dense:
      subtract((const DenseMatrix<R,T>&)A);
      break;
    case Matrix<R,T>::rk:
      subtract((const Rk<R,T>&)A);
      break;
    case Matrix<R,T>::h:
      subtract((const HMatrix<R,T>&)A);
      break;
    default:
      abort();
    }
  }
  inline void operator-=(const Matrix &A)
  {
    subtract(A);
  }
  inline void operator=(MatrixSubtract<R,T> op)
  {
    this->copyFrom(op.A,true);
    subtract(op.B);
  }

  // addproduct
  virtual void addproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B)=0;
  virtual void addproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)=0;
  virtual void addproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)=0;
  virtual void addproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)=0;
  virtual void addproduct(const Rk<R,T> &A, const Rk<R,T> &B)=0;
  virtual void addproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)=0;
  virtual void addproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)=0;
  virtual void addproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)=0;
  virtual void addproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B)=0;
  inline void addproduct(const Matrix<R,T> &A, const Matrix<R,T> &B)
  {
    switch(A.type) {

    case Matrix<R,T>::dense:
      switch(B.type) {
      case Matrix<R,T>::dense:
        addproduct((const DenseMatrix<R,T>&)A,(const DenseMatrix<R,T>&)B); break;
      case Matrix<R,T>::rk:
        addproduct((const DenseMatrix<R,T>&)A,(const Rk<R,T>&)B); break;
      case Matrix<R,T>::h:
        addproduct((const DenseMatrix<R,T>&)A,(const HMatrix<R,T>&)B); break;
      default:
        abort();
      }
      break;

    case Matrix<R,T>::rk:
      switch(B.type) {
      case Matrix<R,T>::dense:
        addproduct((const Rk<R,T>&)A,(const DenseMatrix<R,T>&)B); break;
      case Matrix<R,T>::rk:
        addproduct((const Rk<R,T>&)A,(const Rk<R,T>&)B); break;
      case Matrix<R,T>::h:
        addproduct((const Rk<R,T>&)A,(const HMatrix<R,T>&)B); break;
      default:
        abort();
      }
      break;

    case Matrix<R,T>::h:
      switch(B.type) {
      case Matrix<R,T>::dense:
        addproduct((const HMatrix<R,T>&)A,(const DenseMatrix<R,T>&)B); break;
      case Matrix<R,T>::rk:
        addproduct((const HMatrix<R,T>&)A,(const Rk<R,T>&)B); break;
      case Matrix<R,T>::h:
        addproduct((const HMatrix<R,T>&)A,(const HMatrix<R,T>&)B); break;
      default:
        abort();
      }
      break;
      
    default:
      abort();
    }
  }


  // subtractproduct
  virtual void subtractproduct(const DenseMatrix<R,T> &A, const DenseMatrix<R,T> &B)=0;
  virtual void subtractproduct(const DenseMatrix<R,T> &A, const Rk<R,T> &B)=0;
  virtual void subtractproduct(const DenseMatrix<R,T> &A, const HMatrix<R,T> &B)=0;
  virtual void subtractproduct(const Rk<R,T> &A, const DenseMatrix<R,T> &B)=0;
  virtual void subtractproduct(const Rk<R,T> &A, const Rk<R,T> &B)=0;
  virtual void subtractproduct(const Rk<R,T> &A, const HMatrix<R,T> &B)=0;
  virtual void subtractproduct(const HMatrix<R,T> &A, const DenseMatrix<R,T> &B)=0;
  virtual void subtractproduct(const HMatrix<R,T> &A, const Rk<R,T> &B)=0;
  virtual void subtractproduct(const HMatrix<R,T> &A, const HMatrix<R,T> &B)=0;
  inline void subtractproduct(const Matrix<R,T> &A, const Matrix<R,T> &B)
  {
    switch(A.type) {

    case Matrix<R,T>::dense:
      switch(B.type) {
      case Matrix<R,T>::dense:
        subtractproduct((const DenseMatrix<R,T>&)A,(const DenseMatrix<R,T>&)B); break;
      case Matrix<R,T>::rk:
        subtractproduct((const DenseMatrix<R,T>&)A,(const Rk<R,T>&)B); break;
      case Matrix<R,T>::h:
        subtractproduct((const DenseMatrix<R,T>&)A,(const HMatrix<R,T>&)B); break;
      default:
        abort();
      }
      break;

    case Matrix<R,T>::rk:
      switch(B.type) {
      case Matrix<R,T>::dense:
        subtractproduct((const Rk<R,T>&)A,(const DenseMatrix<R,T>&)B); break;
      case Matrix<R,T>::rk:
        subtractproduct((const Rk<R,T>&)A,(const Rk<R,T>&)B); break;
      case Matrix<R,T>::h:
        subtractproduct((const Rk<R,T>&)A,(const HMatrix<R,T>&)B); break;
      default:
        abort();
      }
      break;

    case Matrix<R,T>::h:
      switch(B.type) {
      case Matrix<R,T>::dense:
        subtractproduct((const HMatrix<R,T>&)A,(const DenseMatrix<R,T>&)B); break;
      case Matrix<R,T>::rk:
        subtractproduct((const HMatrix<R,T>&)A,(const Rk<R,T>&)B); break;
      case Matrix<R,T>::h:
        subtractproduct((const HMatrix<R,T>&)A,(const HMatrix<R,T>&)B); break;
      default:
        abort();
      }
      break;
      
    default:
      abort();
    }
  }

  
  virtual void forwardSubstL(const DenseMatrix<R,T> &A)=0;
  virtual void forwardSubstL(const Rk<R,T> &A)=0;
  virtual void forwardSubstL(const HMatrix<R,T> &A)=0;
  inline void forwardSubstL(const Matrix<R,T> &A) {
    switch(A.type) {
    case Matrix<R,T>::dense:
      forwardSubstL((const DenseMatrix<R,T>&)A); break;
    case Matrix<R,T>::rk:
      forwardSubstL((const Rk<R,T>&)A); break;
    case Matrix<R,T>::h:
      forwardSubstL((const HMatrix<R,T>&)A); break;
    default:
      abort();
    }
  }
  
  virtual void forwardSubstU(const DenseMatrix<R,T> &A)=0;
  virtual void forwardSubstU(const Rk<R,T> &A)=0;
  virtual void forwardSubstU(const HMatrix<R,T> &A)=0;
  inline void forwardSubstU(const Matrix<R,T> &A) {
    switch(A.type) {
    case Matrix<R,T>::dense:
      forwardSubstU((const DenseMatrix<R,T>&)A); break;
    case Matrix<R,T>::rk:
      forwardSubstU((const Rk<R,T>&)A); break;
    case Matrix<R,T>::h:
      forwardSubstU((const HMatrix<R,T>&)A); break;
    default:
      abort();
    }
  }

  virtual void backwardSubstU(const DenseMatrix<R,T> &A)=0;
  virtual void backwardSubstU(const Rk<R,T> &A)=0;
  virtual void backwardSubstU(const HMatrix<R,T> &A)=0;
  inline void backwardSubstU(const Matrix<R,T> &A) {
    switch(A.type) {
    case Matrix<R,T>::dense:
      backwardSubstU((const DenseMatrix<R,T>&)A); break;
    case Matrix<R,T>::rk:
      backwardSubstU((const Rk<R,T>&)A); break;
    case Matrix<R,T>::h:
      backwardSubstU((const HMatrix<R,T>&)A); break;
    default:
      abort();
    }
  }

  template<class T1, class T2>
  inline void operator+=(const MatrixMultiply<R,T1,T2> &op)
  {
    addproduct(op.A,op.B);
  }

  template<class T1, class T2>
  inline void operator-=(const MatrixMultiply<R,T1,T2> &op)
  {
    subtractproduct(op.A,op.B);
  }

  virtual void LU()=0;
};

typedef Matrix<Real,Complex> ComplexMatrix;

#endif
