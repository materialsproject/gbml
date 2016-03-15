#include "mutil.h"
#include <stdio.h>

static double *res, *resb, *orig, rmin, rmax;
static int ct0;

void sphM(M,r,u)
double *M, r, *u;
{ double h, u1[3], u2[3];

  /* set the orthogonal unit vectors. */
  h = sqrt(u[0]*u[0]+u[1]*u[1]);
  if (h<=0)
  { u1[0] = u2[1] = 1.0;
    u1[1] = u1[2] = u2[0] = u2[2] = 0.0;
  }
  else
  { u1[0] = u[1]/h; u1[1] = -u[0]/h; u1[2] = 0.0;
    u2[0] = u[2]*u[0]/h; u2[1] = u[2]*u[1]/h; u2[2] = -h;
  }

  /* parameterize the sphere as r(cos(t)cos(v)u + sin(t)u1 + cos(t)sin(v)u2).
   * first layer of M is (dx/dt, dx/dv, dx/dr) at t=v=0.
   */
  M[0] = r*u1[0]; M[1] = r*u1[1]; M[2] = r*u1[2];
  M[3] = r*u2[0]; M[4] = r*u2[1]; M[5] = r*u2[2];
  M[6] = u[0]; M[7] = u[1]; M[8] = u[2];

  /* next layers are second derivative matrix of components of x(r,t,v).
   * d^2x/dt^2 = d^2x/dv^2 = -ru;    d^2x/dtdv = 0;
   * d^2x/drdt = u1; d^2x/drdv = u2; d^2x/dr^2 = 0.
   */

  M[9] = M[13] = -r*u[0];
  M[11]= M[15] = u1[0];
  M[14]= M[16] = u2[0];
  M[10]= M[12] = M[17] = 0.0;

  M[18]= M[22] = -r*u[1];
  M[20]= M[24] = u1[1];
  M[23]= M[25] = u2[1];
  M[19]= M[21] = M[26] = 0.0;

  M[27]= M[31] = -r*u[1];
  M[29]= M[33] = u1[1];
  M[32]= M[34] = u2[1];
  M[28]= M[30] = M[35] = 0.0;

}

double ip3(a,b)
double *a, *b;
{ return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void rn3(a)
double *a;
{ double s;
  s = sqrt(ip3(a,a));
  a[0] /= s; a[1] /= s; a[2] /= s;
}

double sptarea(a,b,c)
double *a, *b, *c;
{ double ea, eb, ec, yab, yac, ybc, sab, sac, sbc;
  double ab[3], ac[3], bc[3], x1[3], x2[3];

  ab[0] = a[0]-b[0]; ab[1] = a[1]-b[1]; ab[2] = a[2]-b[2];
  ac[0] = a[0]-c[0]; ac[1] = a[1]-c[1]; ac[2] = a[2]-c[2];
  bc[0] = b[0]-c[0]; bc[1] = b[1]-c[1]; bc[2] = b[2]-c[2];
 
  yab = ip3(ab,a); yac = ip3(ac,a); ybc = ip3(bc,b);

  x1[0] = ab[0] - yab*a[0]; x2[0] = ac[0] - yac*a[0];
  x1[1] = ab[1] - yab*a[1]; x2[1] = ac[1] - yac*a[1];
  x1[2] = ab[2] - yab*a[2]; x2[2] = ac[2] - yac*a[2];
  sab = ip3(x1,x1); sac = ip3(x2,x2);
  ea = acos(ip3(x1,x2)/sqrt(sab*sac));

  x1[0] = ab[0] + yab*b[0]; x2[0] = bc[0] - ybc*b[0];
  x1[1] = ab[1] + yab*b[1]; x2[1] = bc[1] - ybc*b[1];
  x1[2] = ab[2] + yab*b[2]; x2[2] = bc[2] - ybc*b[2];
  sbc = ip3(x2,x2);
  eb = acos(ip3(x1,x2)/sqrt(sab*sbc));

  x1[0] = ac[0] + yac*c[0]; x2[0] = bc[0] + ybc*c[0];
  x1[1] = ac[1] + yac*c[1]; x2[1] = bc[1] + ybc*c[1];
  x1[2] = ac[2] + yac*c[2]; x2[2] = bc[2] + ybc*c[2];
  ec = acos(ip3(x1,x2)/sqrt(sac*sbc));

/*
 * Euler's formula is a+b+c-PI, except I've cheated...
 * a=ea, c=ec, b=PI-eb, which is more stable.
 */
  return(ea+ec-eb);
}

void li(x,f,fb,mint,ar)
double *x, ar;
int (*f)(), (*fb)(), mint;
{ int i, j, nr=0, nrb, ct1, w;
  double u[3], r, M[36];
  double sres[MXRESULT], tres[MXRESULT];

/* divide mint by 2, and force to even (Simpson's rule...)
 * to make comparable with rectangular interpretation of mint
 */
  mint <<= 1;
  if (mint&1) mint++;

  ct1 = 0;
  for (i= (rmin==0) ? 1 : 0; i<=mint; i++)
  {
    r = rmin + (rmax-rmin)*i/mint;
    w = 2+2*(i&1)-(i==0)-(i==mint);
    u[0] = orig[0]+x[0]*r;
    u[1] = orig[1]+x[1]*r;
    u[2] = orig[2]+x[2]*r;
    nr = f(u,3,tres,NULL);
    if (ct1==0) setzero(sres,nr);
    for (j=0; j<nr; j++)
      sres[j] += w*r*r*tres[j];
    ct1++;

    if ((fb!=NULL) && (i==mint)) /* boundary */
    { sphM(M,rmax,x);
      nrb = fb(u,3,tres,M);
      if (ct0==0) for (j=0; j<nrb; j++) resb[j] = 0.0;
      for (j=0; j<nrb; j++)
        resb[j] += tres[j]*ar;
    }
  }

  if (ct0==0) for (j=0; j<nr; j++) res[j] = 0.0;
  ct0++;

  for (j=0; j<nr; j++)
    res[j] += sres[j] * ar * (rmax-rmin)/(3*mint);
}

void sphint(f,fb,a,b,c,lev,mint,cent)
double *a, *b, *c;
int (*f)(), (*fb)(), lev, mint, cent;
{ double x[3], ab[3], ac[3], bc[3], ar;
  int i;

  if (lev>1)
  { ab[0] = a[0]+b[0]; ab[1] = a[1]+b[1]; ab[2] = a[2]+b[2]; rn3(ab);
    ac[0] = a[0]+c[0]; ac[1] = a[1]+c[1]; ac[2] = a[2]+c[2]; rn3(ac);
    bc[0] = b[0]+c[0]; bc[1] = b[1]+c[1]; bc[2] = b[2]+c[2]; rn3(bc);
    lev >>= 1;
    if (cent==0)
    { sphint(f,fb,a,ab,ac,lev,mint,1);
      sphint(f,fb,ab,bc,ac,lev,mint,0);
    }
    else
    { sphint(f,fb,a,ab,ac,lev,mint,1);
      sphint(f,fb,b,ab,bc,lev,mint,1);
      sphint(f,fb,c,ac,bc,lev,mint,1);
      sphint(f,fb,ab,bc,ac,lev,mint,1);
    }
    return;
  }

  x[0] = a[0]+b[0]+c[0];
  x[1] = a[1]+b[1]+c[1];
  x[2] = a[2]+b[2]+c[2];
  rn3(x);
  ar = sptarea(a,b,c);

  for (i=0; i<8; i++)
  { if (i>0)
    { x[0] = -x[0];
      if (i%2 == 0) x[1] = -x[1];
      if (i==4) x[2] = -x[2];
    }
    switch(cent)
    { case 2: /* the reflection and its 120', 240' rotations */
        ab[0] = x[0]; ab[1] = x[2]; ab[2] = x[1]; li(ab,f,fb,mint,ar);
        ab[0] = x[2]; ab[1] = x[1]; ab[2] = x[0]; li(ab,f,fb,mint,ar);
        ab[0] = x[1]; ab[1] = x[0]; ab[2] = x[2]; li(ab,f,fb,mint,ar);
      case 1: /* and the 120' and 240' rotations */
        ab[0] = x[1]; ab[1] = x[2]; ab[2] = x[0]; li(ab,f,fb,mint,ar);
        ac[0] = x[2]; ac[1] = x[0]; ac[2] = x[1]; li(ac,f,fb,mint,ar);
      case 0: /* and the triangle itself. */
        li( x,f,fb,mint,ar);
    }
  }
}

void integ_sphere(f,fb,fl,Res,Resb,mg)
double *fl, *Res, *Resb;
int (*f)(), (*fb)(), *mg;
{ double a[3], b[3], c[3];

  a[0] = 1; a[1] = a[2] = 0;
  b[1] = 1; b[0] = b[2] = 0;
  c[2] = 1; c[0] = c[1] = 0;
  
  res = Res;
  resb=Resb;
  orig = &fl[2];
  rmin = fl[0];
  rmax = fl[1];

  ct0 = 0;
  sphint(f,fb,a,b,c,mg[1],mg[0],0);
}
