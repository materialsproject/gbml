/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *   some miscellaneous entry points.
 */

#include "local.h"

void scritval(k0,d,cov,m,rdf,z,k)
double *k0, *z, *cov, *rdf;
Sint *d, *m, *k;
{ int i;
  lf_error = 0;
  for (i=0; i<*k; i++)
    z[i] = critval(1-cov[i], k0, (int)(*m), (int)(*d), TWO_SIDED,*rdf, (*rdf==0) ? GAUSS : TPROC);
}

void slscv(x,n,h,z)
double *x, *h, *z;
int *n;
{ double res[4];
  kdecri(x,*h,res,0.0,3,WGAUS,*n);
  z[0] = res[0];
  z[1] = res[2];
}

void kdeb(x,mi,band,ind,h0,h1,meth,nmeth,ker)
double *x, *band, *h0, *h1;
Sint *mi, *ind, *meth, *nmeth, *ker;
{ int i, imeth[10];
  for (i=0; i<*nmeth; i++) imeth[i] = meth[i];
  kdeselect(band,x,ind,*h0,*h1,imeth,(int)*nmeth,(int)*ker,(int)mi[MN]);
}
