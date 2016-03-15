/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   This file includes functions to solve for the scale estimate in
 *   local robust regression and likelihood. The main entry point is
 *   lf_robust(lfd,sp,des,mxit),
 *   called from the locfit() function.
 *
 *   The update_rs(x) accepts a residual scale x as the argument (actually,
 *   it works on the log-scale). The function computes the local fit
 *   assuming this residual scale, and re-estimates the scale from this
 *   new fit. The final solution satisfies the fixed point equation
 *   update_rs(x)=x. The function lf_robust() automatically calls
 *   update_rs() through the fixed point iterations.
 *
 *   The estimation of the scale from the fit is based on the sqrt of
 *   the median deviance of observations with non-zero weights (in the
 *   gaussian case, this is the median absolute residual).
 *
 *   TODO:
 *     Should use smoothing weights in the median.
 */

#include "local.h"

extern int lf_status;
double robscale;

static lfdata *rob_lfd;
static smpar *rob_sp;
static design *rob_des;
static int rob_mxit;

double median(x,n)
double *x;
int n;
{ int i, j, lt, eq, gt;
  double lo, hi, s;
  lo = hi = x[0];
  for (i=0; i<n; i++)
  { lo = MIN(lo,x[i]);
    hi = MAX(hi,x[i]);
  }
  if (lo==hi) return(lo);
  lo -= (hi-lo);
  hi += (hi-lo);
  for (i=0; i<n; i++)
  { if ((x[i]>lo) & (x[i]<hi))
    { s = x[i]; lt = eq = gt = 0;
      for (j=0; j<n; j++)
      { lt += (x[j]<s);
        eq += (x[j]==s);
        gt += (x[j]>s);
      }
      if ((2*(lt+eq)>n) && (2*(gt+eq)>n)) return(s);
      if (2*(lt+eq)<=n) lo = s;
      if (2*(gt+eq)<=n) hi = s;
    }
  }
  return((hi+lo)/2);
}

double nrobustscale(lfd,sp,des,rs)
lfdata *lfd;
smpar *sp;
design *des;
double rs;
{ int i, ii, p;
  double link[LLEN], sc, sd, sw, e;
  p = des->p; sc = sd = sw = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = base(lfd,ii)+innerprod(des->cf,d_xi(des,i),p);
    e = resp(lfd,ii)-des->th[i];
    stdlinks(link,lfd,sp,ii,des->th[i],rs);
    sc += des->w[i]*e*link[ZDLL];
    sd += des->w[i]*e*e*link[ZDDLL];
    sw += des->w[i];
  }

  /* newton-raphson iteration for log(s)
     -psi(ei/s) - log(s); s = e^{-th}
  */
  rs *= exp((sc-sw)/(sd+sc));
  return(rs);
}

double robustscale(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int i, ii, p, fam, lin;
  double rs, link[LLEN];
  p = des->p;
  fam = fam(sp);
  lin = link(sp);
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = base(lfd,ii) + innerprod(des->cf,d_xi(des,i),p);
    links(des->th[i],resp(lfd,ii),fam&127,lin,link,cens(lfd,ii),prwt(lfd,ii),1.0);
    des->res[i] = -2*link[ZLIK];
  }
  rs = sqrt(median(des->res,des->n));
  if (rs==0.0) rs = 1.0;
  return(rs);
}

double update_rs(x)
double x;
{ double nx;
  if (lf_status != LF_OK) return(x);
  robscale = exp(x);
  lfiter(rob_des,rob_mxit);
  if (lf_status != LF_OK) return(x);

  nx = log(robustscale(rob_lfd,rob_sp,rob_des));
  if (nx<x-0.2) nx = x-0.2;
  return(nx);
}

void lf_robust(lfd,sp,des,mxit)
lfdata *lfd;
design *des;
smpar *sp;
int mxit;
{ double x;
  rob_lfd = lfd;
  rob_des = des;
  rob_sp = sp;
  rob_mxit = mxit;
  lf_status = LF_OK;

  x = log(robustscale(lfd,sp,des));
  solve_fp(update_rs, x, 1.0e-6, mxit);
}
