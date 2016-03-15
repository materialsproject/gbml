/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  Routines for computing weight diagrams.
 *     wdiag(lf,des,lx,deg,ty,exp)
 *  Must locfit() first, unless ker==WPARM and has par. comp.
 *  
 */

#include "local.h"

static double *wd;
extern double robscale;
void nnresproj(lfd,sp,des,u,m,p)
lfdata *lfd;
smpar *sp;
design *des;
double *u;
int m, p;
{ int i, j;
  double link[LLEN];
  setzero(des->f1,p);
  for (j=0; j<m; j++)
  { stdlinks(link,lfd,sp,(int)des->ind[j],des->th[j],robscale);
    for (i=0; i<p; i++) des->f1[i] += link[ZDDLL]*d_xij(des,j,i)*u[j];
  }
  jacob_solve(&des->xtwx,des->f1);
  for (i=0; i<m; i++)
    u[i] -= innerprod(des->f1,d_xi(des,i),p)*des->w[i];
}

void wdexpand(l,n,ind,m)
double *l;
Sint *ind;
int n, m;
{ int i, j, t;
  double z;
  for (j=m; j<n; j++) { l[j] = 0.0; ind[j] = -1; }
  j = m-1;
  while (j>=0)
  { if (ind[j]==j) j--;
    else
    { i = ind[j];
      z = l[j]; l[j] = l[i]; l[i] = z;
      t = ind[j]; ind[j] = ind[i]; ind[i] = t;
      if (ind[j]==-1) j--;
    }
  }

/*  for (i=n-1; i>=0; i--)
  { l[i] = ((j>=0) && (ind[j]==i)) ? l[j--] : 0.0; } */
}

int wdiagp(lfd,sp,des,lx,pc,dv,deg,ty,exp)
lfdata *lfd;
smpar *sp;
design *des;
paramcomp *pc;
deriv *dv;
double *lx;
int deg, ty, exp;
{ int i, j, p, nd;
  double *l1;

  p = des->p;

  fitfun(lfd,sp,des->xev,pc->xbar,des->f1,dv);
  if (exp)
  { jacob_solve(&pc->xtwx,des->f1);
    for (i=0; i<lfd->n; i++)
      lx[i] = innerprod(des->f1,d_xi(des,i),p);
    return(lfd->n);
  }
  jacob_hsolve(&pc->xtwx,des->f1);
  for (i=0; i<p; i++) lx[i] = des->f1[i];

  nd = dv->nd;
  dv->nd = nd+1;
  if (deg>=1)
    for (i=0; i<lfd->d; i++)
    { dv->deriv[nd] = i;
      l1 = &lx[(i+1)*p];
      fitfun(lfd,sp,des->xev,pc->xbar,l1,dv);
      jacob_hsolve(&pc->xtwx,l1);
    }

  dv->nd = nd+2;
  if (deg>=2)
    for (i=0; i<lfd->d; i++)
    { dv->deriv[nd] = i;
      for (j=0; j<lfd->d; j++)
      { dv->deriv[nd+1] = j;
        l1 = &lx[(i*lfd->d+j+lfd->d+1)*p];
        fitfun(lfd,sp,des->xev,pc->xbar,l1,dv);
        jacob_hsolve(&pc->xtwx,l1);
    } }
  dv->nd = nd;
  return(p);
}

int wdiag(lfd,sp,des,lx,dv,deg,ty,exp)
lfdata *lfd;
smpar *sp;
design *des;
deriv *dv;
double *lx;
int deg, ty, exp;
/* deg=0: l(x) only.
   deg=1: l(x), l'(x)
   deg=2: l(x), l'(x), l''(x)
   ty = 1: e1 (X^T WVX)^{-1} X^T W        -- hat matrix
   ty = 2: e1 (X^T WVX)^{-1} X^T WV^{1/2} -- scb's
*/
{ double w, *X, *lxd=NULL, *lxdd=NULL, wdd, wdw, *ulx, link[LLEN], h;
  double dfx[MXDIM], hs[MXDIM];
  int i, ii, j, k, l, m, d, p, nd;

  h = des->h;
  nd = dv->nd;
  wd = des->wd;
  d = lfd->d; p = des->p; X = d_x(des);
  ulx = des->res;
  m = des->n;
  for (i=0; i<d; i++) hs[i] = h*lfd->sca[i];
  if (deg>0)
  { lxd = &lx[m];
    setzero(lxd,m*d);
    if (deg>1)
    { lxdd = &lxd[d*m];
      setzero(lxdd,m*d*d);
  } }

  if (nd>0) fitfun(lfd,sp,des->xev,des->xev,des->f1,dv); /* c(0) */
    else unitvec(des->f1,0,p);
  jacob_solve(&des->xtwx,des->f1);   /* c(0) (X^TWX)^{-1} */
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    lx[i] = innerprod(des->f1,&X[i*p],p); /* c(0)(XTWX)^{-1}X^T */
    if (deg>0)
    { wd[i] = Wd(des->di[ii]/h,ker(sp));
      for (j=0; j<d; j++)
      { dfx[j] = datum(lfd,j,ii)-des->xev[j];
        lxd[j*m+i] = lx[i]*des->w[i]*weightd(dfx[j],lfd->sca[j],
          d,ker(sp),kt(sp),h,lfd->sty[j],des->di[ii]);
             /* c(0) (XTWX)^{-1}XTW' */
      }
      if (deg>1)
      { wdd = Wdd(des->di[ii]/h,ker(sp));
        for (j=0; j<d; j++)
          for (k=0; k<d; k++)
          { w = (des->di[ii]==0) ? 0 : h/des->di[ii];
            w = wdd * (des->xev[k]-datum(lfd,k,ii)) * (des->xev[j]-datum(lfd,j,ii))
                  * w*w / (hs[k]*hs[k]*hs[j]*hs[j]);
            if (j==k) w += wd[i]/(hs[j]*hs[j]);
            lxdd[(j*d+k)*m+i] = lx[i]*w;
              /* c(0)(XTWX)^{-1}XTW'' */
          }
      }
    }
    lx[i] *= des->w[i];
  }

  dv->nd = nd+1;
  if (deg==2)
  { for (i=0; i<d; i++)
    { dv->deriv[nd] = i;
      fitfun(lfd,sp,des->xev,des->xev,des->f1,dv);
      for (k=0; k<m; k++)
      { stdlinks(link,lfd,sp,(int)des->ind[k],des->th[k],robscale);
        for (j=0; j<p; j++)
          des->f1[j] -= link[ZDDLL]*lxd[i*m+k]*X[k*p+j];
        /* c'(x)-c(x)(XTWX)^{-1}XTW'X */
      }
      jacob_solve(&des->xtwx,des->f1); /* (...)(XTWX)^{-1} */
      for (j=0; j<m; j++)
        ulx[j] = innerprod(des->f1,&X[j*p],p); /* (...)XT */
      for (j=0; j<d; j++)
        for (k=0; k<m; k++)
        { ii = des->ind[k];
          dfx[j] = datum(lfd,j,ii)-des->xev[j];
          wdw = des->w[k]*weightd(dfx[j],lfd->sca[j],d,ker(sp),
            kt(sp),h,lfd->sty[j],des->di[ii]);
          lxdd[(i*d+j)*m+k] += ulx[k]*wdw;
          lxdd[(j*d+i)*m+k] += ulx[k]*wdw;
        }
        /* + 2(c'-c(XTWX)^{-1}XTW'X)(XTWX)^{-1}XTW' */
    }
    for (j=0; j<d*d; j++) nnresproj(lfd,sp,des,&lxdd[j*m],m,p);
        /* * (I-X(XTWX)^{-1} XTW */
  }
  if (deg>0)
  { for (j=0; j<d; j++) nnresproj(lfd,sp,des,&lxd[j*m],m,p);
      /* c(0)(XTWX)^{-1}XTW'(I-X(XTWX)^{-1}XTW) */
    for (i=0; i<d; i++)
    { dv->deriv[nd]=i;
      fitfun(lfd,sp,des->xev,des->xev,des->f1,dv);
      jacob_solve(&des->xtwx,des->f1);
      for (k=0; k<m; k++)
        for (l=0; l<p; l++)
          lxd[i*m+k] += des->f1[l]*X[k*p+l]*des->w[k];
            /* add c'(0)(XTWX)^{-1}XTW */
    }
  }

  dv->nd = nd+2;
  if (deg==2)
  { for (i=0; i<d; i++)
    { dv->deriv[nd]=i;
      for (j=0; j<d; j++)
      { dv->deriv[nd+1]=j;
        fitfun(lfd,sp,des->xev,des->xev,des->f1,dv);
        jacob_solve(&des->xtwx,des->f1);
        for (k=0; k<m; k++)
          for (l=0; l<p; l++)
            lxdd[(i*d+j)*m+k] += des->f1[l]*X[k*p+l]*des->w[k];
        /* + c''(x)(XTWX)^{-1}XTW */
      }
    }
  }
  dv->nd = nd;

  k = 1+d*(deg>0)+d*d*(deg==2);

  if (exp) wdexpand(lx,lfd->n,des->ind,m);
 
  if (ty==1) return(m);
  for (i=0; i<m; i++)
  { stdlinks(link,lfd,sp,(int)des->ind[i],des->th[i],robscale);
    link[ZDDLL] = sqrt(fabs(link[ZDDLL]));
    for (j=0; j<k; j++) lx[j*m+i] *= link[ZDDLL];
  }
  return(m);
}
