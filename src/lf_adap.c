/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

/*
  Functions implementing the adaptive bandwidth selection.
  Will make the final call to nbhd() to set smoothing weights
  for selected bandwidth, But will **not** make the
  final call to locfit().
*/

#include "local.h"

static double hmin;

double adcri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ double y;
/* return(-2*lk/(t0*exp(pen*log(1-t2/t0)))); */
  /* return((-2*lk+pen*t2)/t0); */
  y = (MAX(-2*lk,t0-t2)+pen*t2)/t0;
  return(y);
}

double mmse(lfd,sp,dv,des)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
{ int i, ii, j, p, p1;
  double sv, sb, *l, dp;

  l = des->wd;
  wdiag(lfd, sp, des,l,dv,0,1,0);
  sv = sb = 0;
  p = npar(sp);
  for (i=0; i<des->n; i++)
  { sv += l[i]*l[i];
    ii = des->ind[i];
    dp = des->di[ii];
    for (j=0; j<deg(sp); j++) dp *= des->di[ii];
    sb += fabs(l[i])*dp;
  }
  p1 = factorial(deg(sp)+1);
  return(sv+sb*sb*pen(sp)*pen(sp)/(p1*p1));
}

static double mcp, clo, cup;

/*
  Initial bandwidth will be (by default)
  k-nearest neighbors for k small, just large enough to
  get defined estimate (unless user provided nonzero nn or fix-h components)
*/

int ainitband(lfd,sp,dv,des)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
{ int lf_status=0, p, z, cri, noit, redo;
  double ho, t[6];

  if (lf_debug >= 2) printf("ainitband:\n");
  p = des->p;
  cri = acri(sp);
  noit = (cri!=AOK);
  z = (int)(lfd->n*nn(sp));
  if ((noit) && (z<p+2)) z = p+2;
  redo = 0; ho = -1;
  do
  { 
    nbhd(lfd,des,z,redo,sp);
    if (z<des->n) z = des->n;
    if (des->h>ho) lf_status = locfit(lfd,des,sp,noit,0,0);
    z++;
    redo = 1;
  } while ((z<=lfd->n) && ((des->h==0)||(lf_status!=LF_OK)));
  hmin = des->h;

  switch(cri)
  { case ACP:
      local_df(lfd,sp,des,t);
      mcp = adcri(des->llk,t[0],t[2],pen(sp));
      return(lf_status);
    case AKAT:
      local_df(lfd,sp,des,t);
      clo = des->cf[0]-pen(sp)*t[5];
      cup = des->cf[0]+pen(sp)*t[5];
      return(lf_status);
    case AMDI:
      mcp = mmse(lfd,sp,dv,des);
      return(lf_status);
    case AOK: return(lf_status);
  }
  ERROR(("aband1: unknown criterion"));
  return(LF_ERR);
}

/*
  aband2 increases the initial bandwidth until lack of fit results,
  or the fit is close to a global fit. Increase h by 1+0.3/d at
  each iteration.
*/

double aband2(lfd,sp,dv,des,h0)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
double h0;
{ double t[6], h1, nu1, cp, ncp, tlo, tup;
  int d, inc, n, p, done;

  if (lf_debug >= 2) printf("aband2:\n");
  d = lfd->d; n = lfd->n; p = npar(sp);
  h1 = des->h = h0;
  done = 0; nu1 = 0.0;
  inc = 0; ncp = 0.0;
  while ((!done) & (nu1<(n-p)*0.95))
  { fixh(sp) = (1+0.3/d)*des->h;
    nbhd(lfd,des,0,1,sp);
    if (locfit(lfd,des,sp,1,0,0) > 0) WARN(("aband2: failed fit"));
    local_df(lfd,sp,des,t);
    nu1 = t[0]-t[2]; /* tr(A) */
    switch(acri(sp))
    { case AKAT:
        tlo = des->cf[0]-pen(sp)*t[5];
        tup = des->cf[0]+pen(sp)*t[5];
/* printf("h %8.5f  tlo %8.5f  tup %8.5f\n",des->h,tlo,tup); */
        done = ((tlo>cup) | (tup<clo));
        if (!done)
        { clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
          h1 = des->h;
        }
        break;
      case ACP:
        cp = adcri(des->llk,t[0],t[2],pen(sp));
/* printf("h %8.5f  lk %8.5f  t0 %8.5f  t2 %8.5f  cp %8.5f\n",des->h,des->llk,t[0],t[2],cp); */
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        if (cp>=ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=10) | ((inc>=3) & ((t[0]-t[2])>=10) & (cp>1.5*mcp));
        break;
      case AMDI:
        cp = mmse(lfd,sp,dv,des);
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        if (cp>ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=3);
        break;
    }
  }
  return(h1);
}

/*
  aband3 does a finer search around best h so far. Try
  h*(1-0.2/d), h/(1-0.1/d), h*(1+0.1/d), h*(1+0.2/d)
*/
double aband3(lfd,sp,dv,des,h0)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
double h0;
{ double t[6], h1, cp, tlo, tup;
  int i, i0, d, n;

  if (lf_debug >= 2) printf("aband3:\n");
  d = lfd->d; n = lfd->n;
  h1 = h0;
  i0 = (acri(sp)==AKAT) ? 1 : -2;
  if (h0==hmin) i0 = 1;

  for (i=i0; i<=2; i++)
  { if (i==0) i++;
    fixh(sp) = h0*(1+0.1*i/d);
    nbhd(lfd,des,0,1,sp);
    if (locfit(lfd,des,sp,1,0,0) > 0) WARN(("aband3: failed fit"));
    local_df(lfd,sp,des,t);
    switch (acri(sp))
    { case AKAT:
        tlo = des->cf[0]-pen(sp)*t[5];
        tup = des->cf[0]+pen(sp)*t[5];
        if ((tlo>cup) | (tup<clo)) /* done */
          i = 2;
        else
        { h1 = des->h;
          clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
        }
        break;
      case ACP:
        cp = adcri(des->llk,t[0],t[2],pen(sp));
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        else
        { if (i>0) i = 2; }
        break;
      case AMDI:
        cp = mmse(lfd,sp,dv,des);
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        else
        { if (i>0) i = 2; }
    }
  }
  return(h1);
}

int alocfit(lfd,sp,dv,des)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
{ int lf_status;
  double h0;

  lf_status = ainitband(lfd,sp,dv,des);
  if (lf_error) return(lf_status);
  if (acri(sp) == AOK) return(lf_status);

  h0 = fixh(sp);
  fixh(sp) = aband2(lfd,sp,dv,des,des->h);
  fixh(sp) = aband3(lfd,sp,dv,des,fixh(sp));
  nbhd(lfd,des,0,1,sp);
  lf_status = locfit(lfd,des,sp,0,0,0);
  fixh(sp) = h0;

  return(lf_status);
}
