/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <stdlib.h>
#include "local.h"
#include "mutil.h"

void svd(x,p,q,d,mxit)  /* svd of square matrix */
double *x, *p, *q;
int d, mxit;
{ int i, j, k, iter, ms, zer;
  double r, u, v, cp, cm, sp, sm, c1, c2, s1, s2, mx;
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) p[i*d+j] = q[i*d+j] = (i==j);
  for (iter=0; iter<mxit; iter++)
  { ms = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<d; j++)
      { s1 = fabs(x[i*d+j]);
        s2 = fabs(x[j*d+i]);
        mx = (s1>s2) ? s1 : s2;
        zer = 1;
        if (mx*mx>1.0e-15*fabs(x[i*d+i]*x[j*d+j]))
        { if (fabs(x[i*(d+1)])<fabs(x[j*(d+1)]))
          { for (k=0; k<d; k++)
            { u = x[i*d+k]; x[i*d+k] = x[j*d+k]; x[j*d+k] = u;
              u = p[k*d+i]; p[k*d+i] = p[k*d+j]; p[k*d+j] = u;
            }
            for (k=0; k<d; k++)
            { u = x[k*d+i]; x[k*d+i] = x[k*d+j]; x[k*d+j] = u;
              u = q[k*d+i]; q[k*d+i] = q[k*d+j]; q[k*d+j] = u;
            }
          }
          cp = x[i*(d+1)]+x[j*(d+1)];
          sp = x[j*d+i]-x[i*d+j];
          r = sqrt(cp*cp+sp*sp);
          if (r>0) { cp /= r; sp /= r; }
              else { cp = 1.0; zer = 0;}
          cm = x[i*(d+1)]-x[j*(d+1)];
          sm = x[i*d+j]+x[j*d+i];
          r = sqrt(cm*cm+sm*sm);
          if (r>0) { cm /= r; sm /= r; }
              else { cm = 1.0; zer = 0;}
          c1 = cm+cp;
          s1 = sm+sp;
          r = sqrt(c1*c1+s1*s1);
          if (r>0) { c1 /= r; s1 /= r; }
              else { c1 = 1.0; zer = 0;}
          if (fabs(s1)>ms) ms = fabs(s1);
          c2 = cm+cp;
          s2 = sp-sm;
          r = sqrt(c2*c2+s2*s2);
          if (r>0) { c2 /= r; s2 /= r; }
              else { c2 = 1.0; zer = 0;}
          for (k=0; k<d; k++)
          { u = x[i*d+k]; v = x[j*d+k];
            x[i*d+k] = c1*u+s1*v;
            x[j*d+k] = c1*v-s1*u;
            u = p[k*d+i]; v = p[k*d+j];
            p[k*d+i] = c1*u+s1*v;
            p[k*d+j] = c1*v-s1*u;
          }
          for (k=0; k<d; k++)
          { u = x[k*d+i]; v = x[k*d+j];
            x[k*d+i] = c2*u-s2*v;
            x[k*d+j] = s2*u+c2*v;
            u = q[k*d+i]; v = q[k*d+j];
            q[k*d+i] = c2*u-s2*v;
            q[k*d+j] = s2*u+c2*v;
          }
          if (zer) x[i*d+j] = x[j*d+i] = 0.0;
          ms = 1;
        }
      }
    if (ms==0) iter=mxit+10;
  }
  if (iter==mxit) WARN(("Warning: svd not converged.\n"));
  for (i=0; i<d; i++)
    if (x[i*d+i]<0)
    { x[i*d+i] = -x[i*d+i];
      for (j=0; j<d; j++) p[j*d+i] = -p[j*d+i];
    }
}

int svdsolve(x,w,P,D,Q,d,tol) /* original X = PDQ^T; comp. QD^{-1}P^T x */
double *x, *w, *P, *D, *Q, tol;
int d;
{ int i, j, rank;
  double mx;
  if (tol>0)
  { mx = D[0];
    for (i=1; i<d; i++) if (D[i*(d+1)]>mx) mx = D[i*(d+1)];
    tol *= mx;
  }
  rank = 0;
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++)
    if (D[i*d+i]>tol)
    { w[i] /= D[i*(d+1)];
      rank++;
    }
  for (i=0; i<d; i++)
  { x[i] = 0.0;
    for (j=0; j<d; j++) x[i] += Q[i*d+j]*w[j];
  }
  return(rank);
}

void hsvdsolve(x,w,P,D,Q,d,tol) /* original X = PDQ^T; comp. D^{-1/2}P^T x */
double *x, *w, *P, *D, *Q, tol;
int d;
{ int i, j;
  double mx;
  if (tol>0)
  { mx = D[0];
    for (i=1; i<d; i++) if (D[i*(d+1)]>mx) mx = D[i*(d+1)];
    tol *= mx;
  }
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++) if (D[i*d+i]>tol) w[i] /= sqrt(D[i*(d+1)]);
  for (i=0; i<d; i++) x[i] = w[i];
}
