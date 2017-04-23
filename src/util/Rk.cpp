#include "Rk.h"
#include "MatrixOperations.h"

template<>
void ComplexRk :: truncate()
{
  ComplexDenseMatrix QA;
  ComplexDenseMatrix RA;
  ComplexDenseMatrix QB;
  ComplexDenseMatrix RB;
  Array<Complex> tauA;
  Array<Complex> tauB;
  
  qr(A,QA,tauA,RA);
  rq(B,RB,QB,tauB);
  
  int kold = k;
  ComplexDenseMatrix RAB = RA * RB;
  k = min(RAB.m,RAB.n);
  ComplexDenseMatrix U(m,k,0);
  ComplexDenseMatrix VT1(k,n,0);
  ComplexDenseMatrix VT(VT1,0,n-RAB.n,k,RAB.n);

  Array<Real> S(k);

  svd(RAB,U,S,VT);
  
  Real Smax = S[0];
  int knew = 0;
  int kmax = k;
  for(int i=0; i<kmax; i++) {
    if(S[i] >= eps * Smax) {
      knew++;
    }
  }
  
  //cout << "svd " << m << " " << n << " " << kold << " " << knew << "\n";

  U.resize(m,knew);
  VT1.resize(knew,n);

  for(int i = 0; i < U.m; i++)
    for(int j = 0; j < knew; j++)
      U(i,j) *= S[j];

  lmq(QA,tauA,U);
  rmq(VT1,QB,tauB);

  A.copyFrom(U,true);
  B.copyFrom(VT1,true);

  k = knew;
}
