#include "Util.h"
#include "MathUtils.h"
#include <stdexcept> 

using namespace std;
#include <iostream>
#include <string>

Real LinearInterpolate(Real *xtab, Real *ytab, int n, Real x)
{
  int k = searchTable(xtab,n,x);
  return LinearInterpolate(xtab, ytab, n, k, x);
}


Real LinearInterpolate(Real *xtab, Real *ytab, int n, int k, Real x)
{
  Real h = xtab[k+1] - xtab[k];
  Real slope = (ytab[k+1] - ytab[k]) / h ;
  return ytab[k] + slope * ( x - xtab[k] );
}

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
Real brent(Real ax,
           Real bx,
           Real cx,
           Real (*f)(Real, void *), 
           void *args,
           Real tol,
           Real *xmin)
{
  int iter;
  Real a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  Real e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x,args);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
        d=p/q;
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u,args);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
        SHFT(fv,fw,fx,fu)
        } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }
  abort();
  *xmin=x;
  return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT


#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
Real dbrent(Real ax,
            Real bx,
            Real cx,
            Real (*f)(Real, void*),
            Real (*df)(Real, void *), 
            void *args,
            Real tol,
            Real *xmin)
{
	int iter,ok1,ok2;
	Real a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
	Real fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x,args);
	dw=dv=dx=(*df)(x,args);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			if (dv != dx) d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*f)(u,args);
		} else {
			u=x+SIGN(tol1,d);
			fu=(*f)(u,args);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du=(*df)(u,args);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	abort();
	return 0.0;
}
#undef ITMAX
#undef ZEPS
#undef MOV3

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(Real *ax,
            Real *bx,
            Real *cx,
            Real *fa,
            Real *fb,
            Real *fc,
            Real (*func)(Real, void *args),
            void *args)
{
  Real ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax,args);
  *fb=(*func)(*bx,args);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx,args);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(max(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u,args);
      if (fu < *fc) {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u,args);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u,args);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,(*func)(u,args))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u,args);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u,args);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

void linmin(Real *p, 
            Real *xi,
            int n,
            Real *fret,
            Real (*func)(Real, void *), 
            void *args,
            PointAndVector *pv,
            Real tol)
{
  Real xx,xmin,fx,fb,fa,bx,ax;
  for(int j=0;j<n;j++) {
    pv->p[j] = p[j];
    pv->xi[j] = xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,func,args);
  *fret=brent(ax,xx,bx,func,args,tol,&xmin);
  for(int j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
}

#define ITMAX 200
#define SQR(x) ((x)*(x))
void powell(Real *p, 
            Real **xi,
            int n,
            Real ftol,
            Real lintol,
            int *iter,
            Real *fret,
            Real (*func)(Real, void *args),
            void *args,
            PointAndVector *pv)
{
  int i,ibig,j;
  Real del,fp,fptt,t,*pt,*ptt,*xit;

  pt = new Real[n];
  ptt = new Real[n];
  xit = new Real[n];

  for(j=0;j<n;j++) {
    pv->p[j] = p[j];
    pv->xi[j] = xi[j][0];
  }
  *fret=(*func)(0.,args);
  for(j=0;j<n;j++) pt[j]=p[j];
  for(*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func,args,pv,lintol);
      if (fabs(fptt-(*fret)) > del) {
        del=fabs(fptt-(*fret));
        ibig=i;
      }
    }
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
      delete [] pt;
      delete [] ptt;
      delete [] xit;
      return;
    }
    if (*iter == ITMAX) abort();
    for (j=0;j<n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
      pv->p[j] = ptt[j];
    }
    fptt=(*func)(0.,args);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t < 0.0) {
        linmin(p,xit,n,fret,func,args,pv,lintol);
        for (j=0;j<n;j++) {
          xi[j][ibig]=xi[j][n-1];
          xi[j][n-1]=xit[j];
        }
      }
    }
  }
}
#undef ITMAX

int searchTable(Real *tab, int n, Real x)
{
  if(n == 1)
    return 0;
  if(n == 2) {
    if(x == tab[1])
      return 1;
    else
      return 0;
  }
  int k = n/2;
  if(x >= tab[k])
    return k + searchTable(tab+k,n-k,x);
  else
    return searchTable(tab,k+1,x);

}

Real cubicHermiteSplineInterpolate(Real *xtab, Real *ytab, int n, Real x)
{
  int k = searchTable(xtab,n,x);
  return cubicHermiteSplineInterpolate(xtab, ytab, n, k, x);
}

Real cubicHermiteSplineInterpolate(Real *xtab, Real *ytab, int n, int k, Real x)
{
  Real m0, m1;
  Real h = xtab[k+1] - xtab[k];

  if(k == 0)
    m0 = (ytab[1] - ytab[0]) / h;
  else
    m0 = 0.5 * ((ytab[k] - ytab[k-1]) / (xtab[k] - xtab[k-1]) + (ytab[k+1] - ytab[k]) / h);

  if(k == n-2)
    m1 = (ytab[n-1] - ytab[n-2]) / (xtab[n-1] - xtab[n-2]);
  else
    m1 = 0.5 * ((ytab[k+1] - ytab[k]) / h + (ytab[k+2] - ytab[k+1]) / (xtab[k+2] - xtab[k+1]));

  Real t = (x - xtab[k]) / h;
  Real t1 = t - 1.;
  Real t12 = t1 * t1;
  Real twot = t + t;
  Real h00 = t12 * (1. + twot);
  Real h10 = t12 * t;
  Real t2 = t * t;
  Real h01 = t2 * (3 - twot);
  Real h11 = t2 * t1;
  
  return h00 * ytab[k] + h10 * h * m0 + h01 * ytab[k+1] + h11 * h * m1;
}

Complex parseComplex(const string &s)
{
  return parseComplex(s.c_str());
}

Complex parseComplex(char *str)
{
  char *s = str;

  while(*s == ' ') s++; 
  if(*s == 0) throw invalid_argument("Invalid Complex argument");
  Real r = atof(s);
  while(*s != '+' && *s != '-') s++; 
  while(*s == '+') s++; 
  if(*s == 0) throw invalid_argument("Invalid Complex argument");
  Real i = atof(s);

  Complex result(r,i);
  return result;
}
