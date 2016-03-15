/*
 *  Integrate a function f over a circle or disc.
 */

#include "mutil.h"
#include <stdio.h>
#ifndef PI
#define PI    3.141592653589793238462643
#endif

void setM(M,r,s,c,b)
double *M, r, s, c;
int b;
{ M[0] =-r*s; M[1] = r*c;
  M[2] = b*c; M[3] = b*s;
  M[4] =-r*c; M[5] = -s;
  M[6] = -s;  M[7] = 0.0;
  M[8] =-r*s; M[9] = c;
  M[10]=   c; M[11]= 0.0;
}

void integ_circ(f,r,orig,res,mint,b)
int (*f)(), mint, b;
double r, *orig, *res;
{ double y, x[2], theta, tres[MXRESULT], M[12], c, s;
  int i, j, nr=0;
  
  y = 0;
  for (i=0; i<mint; i++)
  { theta = 2*PI*(double)i/(double)mint;
    c = cos(theta); s = sin(theta);
    x[0] = orig[0]+r*c;
    x[1] = orig[1]+r*s;

    if (b!=0)
    { M[0] =-r*s; M[1] = r*c;
      M[2] = b*c; M[3] = b*s;
      M[4] =-r*c; M[5] = -s;
      M[6] = -s;  M[7] = 0.0;
      M[8] =-r*s; M[9] = c;
      M[10]=   c; M[11]= 0.0;
    }

    nr = f(x,2,tres,M);
    if (i==0) setzero(res,nr);
    for (j=0; j<nr; j++) res[j] += tres[j];
  }
  y = 2 * PI * ((b==0)?r:1.0) / mint;
  for (j=0; j<nr; j++) res[j] *= y;
}

void integ_disc(f,fb,fl,res,resb,mg)
int (*f)(), (*fb)(), *mg;
double *fl, *res, *resb;
{ double x[2], y, r, tres[MXRESULT], *orig, rmin, rmax, theta, c, s, M[12];
  int ct, ctb, i, j, k, nr, nrb=0, w;

  orig = &fl[2];
  rmax = fl[1];
  rmin = fl[0];
  y = 0.0;
  ct = ctb = 0;

  for (j=0; j<mg[1]; j++)
  { theta = 2*PI*(double)j/(double)mg[1];
    c = cos(theta); s = sin(theta);
    for (i= (rmin>0) ? 0 : 1; i<=mg[0]; i++)
    { r = rmin + (rmax-rmin)*i/mg[0];
      w = (2+2*(i&1)-(i==0)-(i==mg[0]));
      x[0] = orig[0] + r*c;
      x[1] = orig[1] + r*s;
      nr = f(x,2,tres,NULL);
      if (ct==0) setzero(res,nr);
      for (k=0; k<nr; k++) res[k] += w*r*tres[k];
      ct++;
      if (((i==0) | (i==mg[0])) && (fb!=NULL))
      { setM(M,r,s,c,1-2*(i==0));
        nrb = fb(x,2,tres,M);
        if (ctb==0) setzero(resb,nrb);
        ctb++;
        for (k=0; k<nrb; k++) resb[k] += tres[k];
      }
    }
  }


/*  for (i= (rmin>0) ? 0 : 1; i<=mg[0]; i++)
  {
    r = rmin + (rmax-rmin)*i/mg[0];
    w = (2+2*(i&1)-(i==0)-(i==mg[0]));

    for (j=0; j<mg[1]; j++)
    { theta = 2*PI*(double)j/(double)mg[1];
      c = cos(theta); s = sin(theta);
      x[0] = orig[0] + r*c;
      x[1] = orig[1] + r*s;
      nr = f(x,2,tres,NULL);
      if (ct==0) setzero(res,nr);
      ct++;
      for (k=0; k<nr; k++) res[k] += w*r*tres[k];

      if (((i==0) | (i==mg[0])) && (fb!=NULL))
      { setM(M,r,s,c,1-2*(i==0));
        nrb = fb(x,2,tres,M);
        if (ctb==0) setzero(resb,nrb);
        ctb++;
        for (k=0; k<nrb; k++) resb[k] += tres[k];
      }
    }
  } */
  for (j=0; j<nr; j++) res[j] *= 2*PI*(rmax-rmin)/(3*mg[0]*mg[1]);
  for (j=0; j<nrb; j++) resb[j] *= 2*PI/mg[1];
}
