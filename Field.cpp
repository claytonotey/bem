#include "Field.h"
#include "DenseMatrix.h"
#include "HMatrix.h"
#include "MatrixOperations.h"

ComplexDenseMatrix * Bdense(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options)
{
  GreensFunction *g1 = t->g;
  GreensFunction *g2 = gExt;

  ComplexDenseMatrix *Bp = new ComplexDenseMatrix(2*t->rwgs,1,0);
  ComplexDenseMatrix &B = *Bp;

  for(int i=0; i<t->rwgs; i++) {
    int ii = 2 * i;
    RWG *rwg = t->rwg[i];
    for(SourceSet::iterator j=sources.begin(); j != sources.end(); ++j) {
      Source *s = *j;
      if(s->isInRegion(g1) || s->isInRegion(g2)) {
        ComplexPair EH = s->EH(rwg,options);
        if(s->isInRegion(g1)) {
          if(g1->isExternal()) {
            B(ii,0) += EH.first;
            B(ii+1,0) += EH.second;
          } else {
            B(ii,0) -= EH.first;
            B(ii+1,0) -= EH.second;
          }
        } else if(s->isInRegion(g2)) {
          if(g2->isExternal()) {
            B(ii,0) += EH.first;
            B(ii+1,0) += EH.second;
          } else {
            B(ii,0) -= EH.first;
            B(ii+1,0) -= EH.second;
          }
        }
      }
    }
  }
  // progress.finishedRHS(B);
  return Bp;
}

ComplexMatrix * BH(RWGTree *t, GreensFunction *gExt, SourceSet &sources, Progress &progress, Options &options)
{
  if(!t->sub) {
    return Bdense(t,gExt,sources,progress,options);
  }
  ComplexHMatrix *Hp = new ComplexHMatrix(t->subs,1);
  ComplexHMatrix &H = *Hp;
  for(int i=0; i<t->subs; i++) {
    H.r[i] = 2*t->sub[i]->rwgs;
  }
  H.c[0] = 1;
  for(int i=0; i<t->subs; i++) {
    ComplexMatrix *M = BH(t->sub[i],gExt,sources,progress,options);
    H.set(i,0,M);
  }
  return Hp;
}
