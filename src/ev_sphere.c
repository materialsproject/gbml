/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *   Functions for constructing the fit and
 *   interpolating on the circle/sphere. d=2 only.
 */

#include "local.h"

/*
 * Guess the number of fitting points.
 */
void sphere_guessnv(nvm,ncm,vc,mg)
int *nvm, *ncm, *vc, *mg;
{ *nvm = mg[1]*(mg[0]+1);
  *ncm = 0;
  *vc = 0;
}

void sphere_start(des,lf)
design *des;
lfit *lf;
{ int d, i, j, ct, nv, ncm, vc, *mg;
  double rmin, rmax, *orig, r, th, c, s;

  mg = mg(&lf->evs);
  sphere_guessnv(&nv,&ncm,&vc,mg);
  trchck(lf,nv,0,0);
  d = lf->lfd.d;

  rmin = lf->evs.fl[0];
  rmax = lf->evs.fl[1];
  orig = &lf->evs.fl[2];
rmin = 0; rmax = 1; orig[0] = orig[1] = 0.0;

  ct = 0;
  for (i=0; i<mg[1]; i++)
  { th = 2*PI*i/mg[1];
    c = cos(th);
    s = sin(th);
    for (j=0; j<=mg[0]; j++)
    { r = rmin + (rmax-rmin)*j/mg[0];
      evptx(&lf->fp,ct,0) = orig[0] + r*c;
      evptx(&lf->fp,ct,1) = orig[1] + r*s;
      des->vfun(des,lf,ct);
      ct++;
    }
  }
  lf->fp.nv = ct;
  lf->evs.nce = 0;
}

double sphere_int(lf,x,what)
lfit *lf;
double *x;
int what;
{ double rmin, rmax, *orig, dx, dy, r, th, th0, th1;
  double v[64][64], c0, c1, s0, s1, r0, r1, d0, d1;
  double ll[2], ur[2], xx[2];
  int i0, j0, i1, j1, *mg, nc, ce[4];

  rmin = lf->evs.fl[0];
  rmax = lf->evs.fl[1];
  orig = &lf->evs.fl[2];
rmin = 0; rmax = 1; orig[0] = orig[1] = 0.0;
  mg = mg(&lf->evs);

  dx = x[0] - orig[0];
  dy = x[1] - orig[1];
  r = sqrt(dx*dx+dy*dy);
  th = atan2(dy,dx); /* between -pi and pi */

  i0 = (int)floor(mg[1]*th/(2*PI)) % mg[1];
  j0 = (int)(mg[0]*(r-rmin)/(rmax-rmin));

  i1 = (i0+1) % mg[1];
  j1 = j0+1; if (j1>mg[0]) { j0 = mg[0]-1; j1 = mg[0]; }

  ce[0] = i0*(mg[0]+1)+j0;
  ce[1] = i0*(mg[0]+1)+j1;
  ce[2] = i1*(mg[0]+1)+j0;
  ce[3] = i1*(mg[0]+1)+j1;
  nc = exvval(&lf->fp,v[0],ce[0],2,what,1);
  nc = exvval(&lf->fp,v[1],ce[1],2,what,1);
  nc = exvval(&lf->fp,v[2],ce[2],2,what,1);
  nc = exvval(&lf->fp,v[3],ce[3],2,what,1);

  th0 = 2*PI*i0/mg[1]; c0 = cos(th0); s0 = sin(th0);
  th1 = 2*PI*i1/mg[1]; c1 = cos(th1); s1 = sin(th1);
  r0 = rmin + j0*(rmax-rmin)/mg[0];
  r1 = rmin + j1*(rmax-rmin)/mg[0];
  
  d0 = c0*v[0][1] + s0*v[0][2];
  d1 = r0*(c0*v[0][2]-s0*v[0][1]);
  v[0][1] = d0; v[0][2] = d1;

  d0 = c0*v[1][1] + s0*v[1][2];
  d1 = r1*(c0*v[1][2]-s0*v[1][1]);
  v[1][1] = d0; v[1][2] = d1;

  d0 = c1*v[2][1] + s1*v[2][2];
  d1 = r0*(c1*v[2][2]-s1*v[2][1]);
  v[2][1] = d0; v[2][2] = d1;

  d0 = c1*v[3][1] + s1*v[3][2];
  d1 = r1*(c1*v[3][2]-s1*v[3][1]);
  v[3][1] = d0; v[3][2] = d1;

  xx[0] = r; xx[1] = th;
  ll[0] = r0; ll[1] = th0;
  ur[0] = r1; ur[1] = th1;
  return(rectcell_interp(xx,v,ll,ur,2,nc));
}
