#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

template<class T>
class MatrixVectorMultiply {
 public:
  virtual void multiply(Array<T> &x, Array<T> &y) = 0;
  virtual ~MatrixVectorMultiply() {}
};

template<class T1, class T2>
T1 dot(Array<T2> &v, Array<T2> &w)
{
  T1 r = 0.;
  for(int i=0; i<v.n; i++) r += dot(v[i],w[i]);
  return r;
}

template<class T>
class Preconditioner {
 public:
  virtual void solve(Array<T> &x, Array<T> &b) = 0;
  virtual ~Preconditioner() {}
};

template<class T1, class T2>
class DiagonalPreconditioner : public Preconditioner<T2> {
 protected:
  Array<T1> di;

 public:
 DiagonalPreconditioner(Array<T1> &diag) : di(diag.n) {
    int n = diag.n;
    for(int i=0;i<n;i++) {
      di[i] = 1. / diag[i];
    }
  }
  
  void solve(Array<T2> &x, Array<T2> &b) {
    int n = x.n;
    for(int i=0; i<n; i++) {
      x[i] = di[i] * b[i];
    }
  }
};

template<class R, class C, class T>
Real conjugateGradientSquared(MatrixVectorMultiply<T> &A,
                              Preconditioner<T> &M,
                              Array<T> &x, 
                              Array<T> &b,
                              Real eps,
                              int maxiters)
{
  int n = b.n;
  Array<T> r(n);
  Array<T> rt(n);
  Array<T> u(n);
  Array<T> p(n);
  Array<T> q(n);
  Array<T> w1(n);
  Array<T> w2(n);
  Array<T> w3(n);

  R abstol = eps * real(dot<C,T>(b,b));
  A.multiply(x,w1);
  for(int i=0; i<n; i++) {
    r[i] = b[i] - w1[i];
    rt[i] = r[i];
  }
  C rhoPrev = 0.;
  for(int iters = 0; iters < maxiters; iters++) {
    C rho = dot<C,T>(rt,r);
    if(rho == 0.)
      return -1.;
    if(!iters) {
      for(int i=0; i<n; i++) {
        u[i] = r[i];
        p[i] = r[i];
      }
    } else {
      C beta = rho / rhoPrev;
      for(int i=0; i<n; i++) {
        T betaq = beta * q[i];
        u[i] = r[i] + betaq;
        T betap = beta * p[i];
        T qbetap = q[i] + betap;
        T betaqbetap = beta * qbetap;
        p[i] = u[i] + betaqbetap; 
      }
    }
    M.solve(w1,p);
    A.multiply(w1,w2);
    C alpha = rho / dot<C,T>(rt,w2);
    for(int i=0; i<n; i++) {
      T alphaw2 = alpha * w2[i];
      q[i] = u[i] - alphaw2;
      w1[i] = u[i] + q[i];
    }
    M.solve(w3,w1);
    for(int i=0; i<n; i++) {
      T alphaw3 = alpha * w3[i];
      x[i] += alphaw3;
    }
    A.multiply(x,w1);
    for(int i=0; i<n; i++) {
      w2[i] = b[i] - w1[i];
    }
    R res = real(dot<C,T>(w2,w2));
    if(res < abstol) 
      return res;
    A.multiply(w3,w1);
    for(int i=0; i<n; i++) {
      T alphaw1 = alpha * w1[i];
      r[i] -= alphaw1;
    }
    rhoPrev = rho;
  }
  return -1.;
}

#endif
