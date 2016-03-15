/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include "mutil.h"

/* A is a n*p matrix, find the cholesky decomposition
 * of the first p rows. In most applications, will want n=p.
 */
void chol_dec(A,n,p)
double *A;
int n, p;
{ int i, j, k;

  for (j=0; j<p; j++)
  { k = n*j+j;
    for (i=0; i<j; i++) A[k] -= A[n*j+i]*A[n*j+i];
    if (A[k]<=0)
    { for (i=j; i<p; i++) A[n*i+j] = 0.0; }
    else
    { A[k] = sqrt(A[k]);
      for (i=j+1; i<p; i++)
      { for (k=0; k<j; k++)
          A[n*i+j] -= A[n*i+k]*A[n*j+k];
        A[n*i+j] /= A[n*j+j];
      }
    }
  }
  for (j=0; j<p; j++)
    for (i=j+1; i<p; i++) A[n*j+i] = 0.0;
}

int chol_solve(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;

  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[i*n+j]*v[j];
    v[i] /= A[i*n+i];
  }
  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++) v[i] -= A[j*n+i]*v[j];
    v[i] /= A[i*n+i];
  }
  return(p);
}

int chol_hsolve(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;

  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[i*n+j]*v[j];
    v[i] /= A[i*n+i];
  }
  return(p);
}

double chol_qf(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;
  double sum;
 
  sum = 0.0;
  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[i*n+j]*v[j];
    v[i] /= A[i*n+i];
    sum += v[i]*v[i];
  }
  return(sum);
}
