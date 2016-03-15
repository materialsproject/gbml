/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Includes some miscellaneios vector functions:
 *     setzero(v,p)           sets all elements of v to 0.
 *     unitvec(x,k,p)         sets x to k'th unit vector e_k.
 *     innerprod(v1,v2,p)     inner product.
 *     addouter(A,v1,v2,p,c)  A <- A + c * v_1 v2^T
 *     multmatscal(A,z,n)     A <- A*z
 *     transpose(x,m,n)       inline transpose
 *     m_trace(x,n)           trace
 */

#include "mutil.h"

void setzero(v,p)
double *v;
int p;
{ int i;
  for (i=0; i<p; i++) v[i] = 0.0;
}

void unitvec(x,k,p)
double *x;
int k, p;
{ setzero(x,p);
  x[k] = 1.0;
}

double innerprod(v1,v2,p)
double *v1, *v2;
int p;
{ int i;
  double s;
  s = 0;
  for (i=0; i<p; i++) s += v1[i]*v2[i];
  return(s);
}

void addouter(A,v1,v2,p,c)
double *A, *v1, *v2, c;
int p;
{ int i, j;
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      A[i*p+j] += c*v1[i]*v2[j];
}

void multmatscal(A,z,n)
double *A, z;
int n;
{ int i;
  for (i=0; i<n; i++) A[i] *= z;
}

/*
 *  transpose() transposes an m*n matrix in place.
 *  At input, the matrix has n rows, m columns and
 *    x[0..n-1] is the is the first column.
 *  At output, the matrix has m rows, n columns and
 *    x[0..m-1] is the first column.
 */
void transpose(x,m,n)
double *x;
int m, n;
{ int t0, t, ti, tj;
  double z;
  for (t0=1; t0<m*n-2; t0++)
    { ti = t0%m; tj = t0/m;
      do
      { t = ti*n+tj;
        ti= t%m;
        tj= t/m;
      } while (t<t0);
      z = x[t];
      x[t] = x[t0];
      x[t0] = z;
    }
}

/* trace of an n*n square matrix. */
double m_trace(x,n)
double *x;
int n;
{ int i;
  double sum;
  sum = 0;
  for (i=0; i<n; i++)
    sum += x[i*(n+1)];
  return(sum);
}
