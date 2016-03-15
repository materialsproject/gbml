/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  Multivariate integration of a vector-valued function
 *  using Simpson's rule.
 */

#include <math.h>
#include <stdio.h>
#include "mutil.h"
extern void setzero();

static double M[(1+MXIDIM)*MXIDIM*MXIDIM];

/* third order corners */
void simp3(fd,x,d,resd,delta,wt,i0,i1,mg,ct,res2,lfindex)
int (*fd)(), d, wt, i0, i1, *mg, ct, *lfindex;
double *x, *resd, *delta, *res2;
{ int k, l, m, nrd;
  double zb;

  for (k=i1+1; k<d; k++) if ((lfindex[k]==0) | (lfindex[k]==mg[k]))
  {
    setzero(M,d*d);
    m = 0; zb = 1.0;
    for (l=0; l<d; l++)
      if ((l!=i0) & (l!=i1) & (l!=k))
      { M[m*d+l] = 1.0;
        m++;
        zb *= delta[l];
      }
    M[(d-3)*d+i0] = (lfindex[i0]==0) ? -1 : 1;
    M[(d-2)*d+i1] = (lfindex[i1]==0) ? -1 : 1;
    M[(d-1)*d+k] = (lfindex[k]==0) ? -1 : 1;
    nrd = fd(x,d,res2,M);
    if ((ct==0) & (i0==0) & (i1==1) & (k==2)) setzero(resd,nrd);
    for (l=0; l<nrd; l++)
      resd[l] += wt*zb*res2[l];
  }
}

/* second order corners */
void simp2(fc,fd,x,d,resc,resd,delta,wt,i0,mg,ct,res2,lfindex)
int (*fc)(), (*fd)(), d, wt, i0, *mg, ct, *lfindex;
double *x, *resc, *resd, *delta, *res2;
{ int j, k, l, nrc;
  double zb;
  for (j=i0+1; j<d; j++) if ((lfindex[j]==0) | (lfindex[j]==mg[j]))
  { setzero(M,d*d);
    l = 0; zb = 1;
    for (k=0; k<d; k++) if ((k!=i0) & (k!=j))
    { M[l*d+k] = 1.0;
      l++;
      zb *= delta[k];
    }
    M[(d-2)*d+i0] = (lfindex[i0]==0) ? -1 : 1;
    M[(d-1)*d+j] = (lfindex[j]==0) ? -1 : 1;
    nrc = fc(x,d,res2,M);
    if ((ct==0) & (i0==0) & (j==1)) setzero(resc,nrc);
    for (k=0; k<nrc; k++) resc[k] += wt*zb*res2[k];
       
    if (fd!=NULL)
      simp3(fd,x,d,resd,delta,wt,i0,j,mg,ct,res2);
  }
}

/* first order boundary */
void simp1(fb,fc,fd,x,d,resb,resc,resd,delta,wt,mg,ct,res2,lfindex)
int (*fb)(), (*fc)(), (*fd)(), d, wt, *mg, ct, *lfindex;
double *x, *resb, *resc, *resd, *delta, *res2;
{ int i, j, k, nrb;
  double zb;
  for (i=0; i<d; i++) if ((lfindex[i]==0) | (lfindex[i]==mg[i]))
  { setzero(M,(1+d)*d*d);
    k = 0;
    for (j=0; j<d; j++) if (j!=i)
    { M[k*d+j] = 1;
      k++;
    }
    M[(d-1)*d+i] = (lfindex[i]==0) ? -1 : 1;
    nrb = fb(x,d,res2,M);
    zb = 1;
    for (j=0; j<d; j++) if (i!=j) zb *= delta[j];
        if ((ct==0) && (i==0))
          for (j=0; j<nrb; j++) resb[j] = 0.0;
    for (j=0; j<nrb; j++) resb[j] += wt*zb*res2[j];

    if (fc!=NULL)
      simp2(fc,fd,x,d,resc,resd,delta,wt,i,mg,ct,res2,lfindex);
  }
}

void simpson4(f,fb,fc,fd,ll,ur,d,res,resb,resc,resd,mg,res2)
int (*f)(), (*fb)(), (*fc)(), (*fd)(), d, *mg;
double *ll, *ur, *res, *resb, *resc, *resd, *res2;
{ int ct, i, j, nr, wt, lfindex[MXIDIM];
  double x[MXIDIM], delta[MXIDIM], z;

  for (i=0; i<d; i++)
  { lfindex[i] = 0;
    x[i] = ll[i];
    if (mg[i]&1) mg[i]++;
    delta[i] = (ur[i]-ll[i])/(3*mg[i]);
  }
  ct = 0;

  while(1)
  { wt = 1;
    for (i=0; i<d; i++)
      wt *= (4-2*(lfindex[i]%2==0)-(lfindex[i]==0)-(lfindex[i]==mg[i]));
    nr = f(x,d,res2,NULL);
    if (ct==0) setzero(res,nr);
    for (i=0; i<nr; i++) res[i] += wt*res2[i];

    if (fb!=NULL)
      simp1(fb,fc,fd,x,d,resb,resc,resd,delta,wt,mg,ct,res2,lfindex);

    /* compute next grid point */
    for (i=0; i<d; i++)
    { lfindex[i]++;
      if (lfindex[i]>mg[i])
      { lfindex[i] = 0;
        x[i] = ll[i];
        if (i==d-1) /* done */
        { z = 1.0;
          for (j=0; j<d; j++) z *= delta[j];
          for (j=0; j<nr; j++) res[j] *= z;
          return;
        }
      }
      else
      { x[i] = ll[i] + 3*delta[i]*lfindex[i];
        i = d;
      }
    }
    ct++;
  }
}

void simpsonm(f,ll,ur,d,res,mg,res2)
int (*f)(), d, *mg;
double *ll, *ur, *res, *res2;
{ simpson4(f,NULL,NULL,NULL,ll,ur,d,res,NULL,NULL,NULL,mg,res2);
}

double simpson(f,l0,l1,m)
double (*f)(), l0, l1;
int m;
{ double x, sum;
  int i;
  sum = 0;
  for (i=0; i<=m; i++)
  { x = ((m-i)*l0 + i*l1)/m;
    sum += (2+2*(i&1)-(i==0)-(i==m)) * f(x);
  }
  return( (l1-l0) * sum / (3*m) );
}
