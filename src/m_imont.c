/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  Multivariate integration of a vector-valued function
 *  using Monte-Carlo method.
 *
 *  uses drand48() random number generator. Does not seed.
 */

#include <stdlib.h>
//#include "R.h"
#include "mutil.h"
extern void setzero();

/* static int lfindex[MXIDIM];
   static double M[(1+MXIDIM)*MXIDIM*MXIDIM]; */

void monte(f,ll,ur,d,res,n)
int (*f)(), d, n;
double *ll, *ur, *res;
{ int i, j, nr=0;
  double z, x[MXIDIM], tres[MXRESULT];

/* srand48(234L); */
//GetRNGstate(); /* Use R's RNG */  // RJN
  srand48(234L);                    // RJN

  for (i=0; i<n; i++)
//{ for (j=0; j<d; j++) x[j] = ll[j] + (ur[j]-ll[j])*unif_rand(); /* drand48();*/  // RJN
  { for (j=0; j<d; j++) x[j] = ll[j] + (ur[j]-ll[j])*drand48();                    // RJN
    nr = f(x,d,tres,NULL);
    if (i==0) setzero(res,nr);
    for (j=0; j<nr; j++) res[j] += tres[j];
  }

  z = 1;
  for (i=0; i<d; i++) z *= (ur[i]-ll[i]);
  for (i=0; i<nr; i++) res[i] *= z/n;
 
//PutRNGstate();  // RJN
}
