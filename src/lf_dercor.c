/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Derivative corrections. The local slopes are not the derivatives
 *   of the local likelihood estimate; the function dercor() computes
 *   the adjustment to get the correct derivatives under the assumption
 *   that h is constant.
 *
 *   By differentiating the local likelihood equations, one obtains
 *
 *     d ^      ^       T      -1   T  d    .       ^
 *    -- a   =  a  -  (X W V X)    X  -- W  l( Y, X a)
 *    dx  0      1                    dx
 */

#include "local.h"
extern double robscale;

void dercor(lfd,sp,des,coef)
lfdata *lfd;
smpar *sp;
design *des;
double *coef;
{ double s1, dc[MXDIM], wd, link[LLEN];
  int i, ii, j, m, d, p;
  if (fam(sp)<=THAZ) return;
  if (ker(sp)==WPARM) return;

  d = lfd->d;
  p = des->p; m = des->n;

  if (lf_debug>1) printf("  Correcting derivatives\n");
  fitfun(lfd, sp, des->xev,des->xev,des->f1,NULL);
  jacob_solve(&des->xtwx,des->f1);
  setzero(dc,d);

  /* correction term is e1^T (XTWVX)^{-1} XTW' ldot. */
  for (i=0; i<m; i++)
  { s1 = innerprod(des->f1,&des->X[i*p],p);
    ii = des->ind[i];
    stdlinks(link,lfd,sp,ii,des->th[i],robscale);
    for (j=0; j<d; j++)
    { wd = des->w[i]*weightd(datum(lfd,j,ii)-des->xev[j],lfd->sca[j],
        d,ker(sp),kt(sp),des->h,lfd->sty[j],des->di[ii]);
      dc[j] += s1*wd*link[ZDLL];
    }

  }
  for (j=0; j<d; j++) coef[j+1] += dc[j];
}
