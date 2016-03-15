/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include <stdio.h>
#include "mutil.h"

/* qr decomposition of X (n*p organized by column).
 * Take w for the ride, if not NULL.
 */
void qr(X,n,p,w)
double *X, *w;
int n, p;
{ int i, j, k, mi;
  double c, s, mx, nx, t;

  for (j=0; j<p; j++)
  { mi = j;
    mx = fabs(X[(n+1)*j]);
    nx = mx*mx;

    /* find the largest remaining element in j'th column, row mi.
     * flip that row with row j.
     */
    for (i=j+1; i<n; i++)
    { nx += X[j*n+i]*X[j*n+i];
      if (fabs(X[j*n+i])>mx)
      { mi = i;
        mx = fabs(X[j*n+i]);
      }
    }
    for (i=j; i<p; i++)
    { t = X[i*n+j];
      X[i*n+j] = X[i*n+mi];
      X[i*n+mi] = t;
    }
    if (w!=NULL) { t = w[j]; w[j] = w[mi]; w[mi] = t; }

    /* want the diag. element -ve, so we do the `good' Householder reflect.
     */
    if (X[(n+1)*j]>0)
    { for (i=j; i<p; i++) X[i*n+j] = -X[i*n+j];
      if (w!=NULL) w[j] = -w[j];
    }

    nx = sqrt(nx);
    c = nx*(nx-X[(n+1)*j]);
    if (c!=0)
    { for (i=j+1; i<p; i++)
      { s = 0;
        for (k=j; k<n; k++)
          s += X[i*n+k]*X[j*n+k];
        s = (s-nx*X[i*n+j])/c;
        for (k=j; k<n; k++)
          X[i*n+k] -= s*X[j*n+k];
        X[i*n+j] += s*nx;
      }
      if (w != NULL)
      { s = 0;
        for (k=j; k<n; k++)
          s += w[k]*X[n*j+k];
        s = (s-nx*w[j])/c;
        for (k=j; k<n; k++)
          w[k] -= s*X[n*j+k];
        w[j] += s*nx;
      }
      X[j*n+j] = nx;
    }
  }
}

void qrinvx(R,x,n,p)
double *R, *x;
int n, p;
{ int i, j;
  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++) x[i] -= R[j*n+i]*x[j];
    x[i] /= R[i*n+i];
  }
}

void qrtinvx(R,x,n,p)
double *R, *x;
int n, p;
{ int i, j;
  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) x[i] -= R[i*n+j]*x[j];
    x[i] /= R[i*n+i];
  }
}

void qrsolv(R,x,n,p)
double *R, *x;
int n, p;
{ qrtinvx(R,x,n,p);
  qrinvx(R,x,n,p);
}
