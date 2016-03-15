/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Evaluate the locfit fitting functions.
 *     calcp(sp,d)
 *       calculates the number of fitting functions.
 *     makecfn(sp,des,dv,d)
 *       makes the coef.number vector.
 *     fitfun(lfd, sp, x,t,f,dv)
 *       lfd is the local fit structure.
 *       sp  smoothing parameter structure.
 *       x is the data point.
 *       t is the fitting point.
 *       f is a vector to return the results.
 *       dv derivative structure.
 *     designmatrix(lfd, sp, des)
 *       is a wrapper for fitfun to build the design matrix.
 *
 */

#include "local.h"

int calcp(sp,d)
smpar *sp;
int d;
{ int i, k;

  if (ubas(sp))
{ printf("calcp-ubas\n");
  return(npar(sp));
}

  switch (kt(sp))
  { case KSPH:
    case KCE:
      k = 1;
      for (i=1; i<=deg(sp); i++) k = k*(d+i)/i;
      return(k);
    case KPROD: return(d*deg(sp)+1);
    case KLM: return(d);
    case KZEON: return(1);
  }
  ERROR(("calcp: invalid kt %d",kt(sp)));
  return(0);
}

int coefnumber(dv,kt,d,deg)
int kt, d, deg;
deriv *dv;
{ int d0, d1, t;

  if (d==1)
  { if (dv->nd<=deg) return(dv->nd);
    return(-1);
  }

  if (dv->nd==0) return(0);
  if (deg==0) return(-1);
  if (dv->nd==1) return(1+dv->deriv[0]);
  if (deg==1) return(-1);
  if (kt==KPROD) return(-1);

  if (dv->nd==2)
  { d0 = dv->deriv[0]; d1 = dv->deriv[1];
    if (d0<d1) { t = d0; d0 = d1; d1 = t; }
    return((d+1)*(d0+1)-d0*(d0+3)/2+d1);
  }
  if (deg==2) return(-1);

  ERROR(("coefnumber not programmed for nd>=3"));
  return(-1);
}

void makecfn(sp,des,dv,d)
smpar *sp;
design *des;
deriv *dv;
int d;
{ int i, nd;
  
  nd = dv->nd;

  des->cfn[0] = coefnumber(dv,kt(sp),d,deg(sp));
  des->ncoef = 1;
  if (nd >= deg(sp)) return;
  if (kt(sp)==KZEON) return;

  if (d>1)
  { if (nd>=2) return;
    if ((nd>=1) && (kt(sp)==KPROD)) return;
  }

  dv->nd = nd+1;
  for (i=0; i<d; i++)
  { dv->deriv[nd] = i;
    des->cfn[i+1] = coefnumber(dv,kt(sp),d,deg(sp));
  }
  dv->nd = nd;

  des->ncoef = 1+d;
}

void fitfunangl(dx,ff,sca,cd,deg)
double dx, *ff, sca;
int deg, cd;
{
  if (deg>=3) WARN(("Can't handle angular model with deg>=3"));

  switch(cd)
  { case 0:
      ff[0] = 1;
      ff[1] = sin(dx/sca)*sca;
      ff[2] = (1-cos(dx/sca))*sca*sca;
      return;
    case 1:
      ff[0] = 0;
      ff[1] = cos(dx/sca);
      ff[2] = sin(dx/sca)*sca;
      return;
    case 2:
      ff[0] = 0;
      ff[1] = -sin(dx/sca)/sca;
      ff[2] = cos(dx/sca);
      return;
    default: WARN(("Can't handle angular model with >2 derivs"));
  }
}

void fitfun(lfd,sp,x,t,f,dv)
lfdata *lfd;
smpar *sp;
double *x, *t, *f;
deriv *dv;
{ int d, deg, nd, m, i, j, k, ct_deriv[MXDIM];
  double ff[MXDIM][1+MXDEG], dx[MXDIM], *xx[MXDIM];

  if (ubas(sp))
  { for (i=0; i<lfd->d; i++) xx[i] = &x[i];
    i = 0;
    sp->vbasis(xx,t,1,lfd->d,&i,1,npar(sp),f);
    return;
  }

  d = lfd->d;
  deg = deg(sp);
  m = 0;
  nd = (dv==NULL) ? 0 : dv->nd;

  if (kt(sp)==KZEON)
  { f[0] = 1.0;
    return;
  }

  if (kt(sp)==KLM)
  { for (i=0; i<d; i++) f[m++] = x[i];
    return;
  }

  f[m++] = (nd==0);
  if (deg==0) return;

  for (i=0; i<d; i++)
  { ct_deriv[i] = 0;
    dx[i] = (t==NULL) ? x[i] : x[i]-t[i];
  }
  for (i=0; i<nd; i++) ct_deriv[dv->deriv[i]]++;

  for (i=0; i<d; i++)
  { switch(lfd->sty[i])
    {
      case STANGL:
        fitfunangl(dx[i],ff[i],lfd->sca[i],ct_deriv[i],deg(sp));
        break;
      default:
        for (j=0; j<ct_deriv[i]; j++) ff[i][j] = 0.0;
        ff[i][ct_deriv[i]] = 1.0;
        for (j=ct_deriv[i]+1; j<=deg; j++)
          ff[i][j] = ff[i][j-1]*dx[i]/(j-ct_deriv[i]);
    }
  }

/*
 *  Product kernels. Note that if ct_deriv[i] != nd, that implies
 *  there is differentiation wrt another variable, and all components
 *  involving x[i] are 0.
 */
  if ((d==1) || (kt(sp)==KPROD))
  { for (j=1; j<=deg; j++)
      for (i=0; i<d; i++)
        f[m++] = (ct_deriv[i]==nd) ? ff[i][j] : 0.0;
    return;
  }

/*
 *  Spherical kernels with the full polynomial basis.
 *  Presently implemented up to deg=3.
 */
  for (i=0; i<d; i++)
    f[m++] = (ct_deriv[i]==nd) ? ff[i][1] : 0.0;
  if (deg==1) return;

  for (i=0; i<d; i++)
  {
    /* xi^2/2 terms. */
    f[m++] = (ct_deriv[i]==nd) ? ff[i][2] : 0.0;

    /* xi xj terms */
    for (j=i+1; j<d; j++)
      f[m++] = (ct_deriv[i]+ct_deriv[j]==nd) ? ff[i][1]*ff[j][1] : 0.0;
  }
  if (deg==2) return;

  for (i=0; i<d; i++)
  { 
    /* xi^3/6 terms */
    f[m++] = (ct_deriv[i]==nd) ? ff[i][3] : 0.0;

    /* xi^2/2 xk terms */
    for (k=i+1; k<d; k++)
      f[m++] = (ct_deriv[i]+ct_deriv[k]==nd) ? ff[i][2]*ff[k][1] : 0.0;

    /* xi xj xk terms */
    for (j=i+1; j<d; j++)
    { f[m++] = (ct_deriv[i]+ct_deriv[j]==nd) ? ff[i][1]*ff[j][2] : 0.0;
      for (k=j+1; k<d; k++)
        f[m++] = (ct_deriv[i]+ct_deriv[j]+ct_deriv[k]==nd) ?
                    ff[i][1]*ff[j][1]*ff[k][1] : 0.0;
    }
  }
  if (deg==3) return;

  ERROR(("fitfun: can't handle deg=%d for spherical kernels",deg));
}

/*
 *  Build the design matrix. Assumes des->ind contains the indices of
 *  the required data points; des->n the number of points; des->xev
 *  the fitting point.
 */
void designmatrix(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int i, ii, j, p;
  double *X, u[MXDIM];

  X = d_x(des);
  p = des->p;

  if (ubas(sp))
  {
    sp->vbasis(lfd->x,des->xev,lfd->n,lfd->d,des->ind,des->n,p,X);
    return;
  }

  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    for (j=0; j<lfd->d; j++) u[j] = datum(lfd,j,ii);
    fitfun(lfd,sp,u,des->xev,&X[i*p],NULL);
  }
}
