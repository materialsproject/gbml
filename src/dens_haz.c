/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Integration for hazard rate estimation. The functions in this
 *   file are used to evaluate
 *      sum int_0^{Ti} W_i(t,x) A()A()' exp( P() ) dt
 *   for hazard rate models.
 *
 *   These routines assume the weight function is supported on [-1,1].
 *   hasint_sph multiplies by exp(base(lf,i)), which allows estimating
 *   the baseline in a proportional hazards model, when the covariate
 *   effect base(lf,i) is known.
 *
 *   TODO:
 *     hazint_sph, should be able to reduce mint in some cases with
 *       small integration range. onedint could be used for beta-family
 *       (RECT,EPAN,BISQ,TRWT) kernels.
 *     hazint_prod, restrict terms from the sum based on x values.
 *       I should count obs >= max, and only do that integration once.
 */

#include "local.h"

static double ilim[2*MXDIM], *ff, tmax;
static lfdata *haz_lfd;
static smpar  *haz_sp;

/*
 *  hrao returns 0 if integration region is empty.
 *               1 otherwise.
 */
int haz_sph_int(dfx,cf,h,r1)
double *dfx, *cf, h, *r1;
{ double s, t0, t1, wt, th;
  int j, dim, p;
  s = 0; p = npar(haz_sp);
  dim = haz_lfd->d;
  for (j=1; j<dim; j++) s += SQR(dfx[j]/(h*haz_lfd->sca[j]));
  if (s>1) return(0);

  setzero(r1,p*p);
  t1 = sqrt(1-s)*h*haz_lfd->sca[0];
  t0 = -t1;
  if (t0<ilim[0])   t0 = ilim[0];
  if (t1>ilim[dim]) t1 = ilim[dim];
  if (t1>dfx[0]) t1 = dfx[0];
  if (t1<t0) return(0);

/*  Numerical integration by Simpson's rule.
 */
  for (j=0; j<=de_mint; j++)
  { dfx[0] = t0+(t1-t0)*j/de_mint;
    wt = weight(haz_lfd, haz_sp, dfx, NULL, h, 0, 0.0);
    fitfun(haz_lfd, haz_sp, dfx,NULL,ff,NULL);
    th = innerprod(cf,ff,p);
    if (link(haz_sp)==LLOG) th = exp(th);
    wt *= 2+2*(j&1)-(j==0)-(j==de_mint);
    addouter(r1,ff,ff,p,wt*th);
  }
  multmatscal(r1,(t1-t0)/(3*de_mint),p*p);

  return(1);
}

int hazint_sph(t,resp,r1,cf,h)
double *t, *resp, *r1, *cf, h;
{ int i, j, n, p, st;
  double dfx[MXDIM], eb, sb;
  p = npar(haz_sp);
  setzero(resp,p*p);
  sb = 0.0;

  n = haz_lfd->n;
  for (i=0; i<=n; i++)
  {
    if (i==n)
    { dfx[0] = tmax-t[0];
      for (j=1; j<haz_lfd->d; j++) dfx[j] = 0.0;
      eb = exp(sb/n);
    }
    else
    { eb = exp(base(haz_lfd,i)); sb += base(haz_lfd,i);
      for (j=0; j<haz_lfd->d; j++) dfx[j] = datum(haz_lfd,j,i)-t[j];
    }

    st = haz_sph_int(dfx,cf,h,r1);
    if (st)
      for (j=0; j<p*p; j++) resp[j] += eb*r1[j];
  }
  return(LF_OK);
}

int hazint_prod(t,resp,x,cf,h)
double *t, *resp, *x, *cf, h;
{ int d, p, i, j, k, st;
  double dfx[MXDIM], t_prev,
         hj, hs, ncf[MXDEG], ef, il1;
  double prod_wk[MXDIM][2*MXDEG+1], eb, sb;

  p = npar(haz_sp);
  d = haz_lfd->d;
  setzero(resp,p*p);
  hj = hs = h*haz_lfd->sca[0];

  ncf[0] = cf[0];
  for (i=1; i<=deg(haz_sp); i++)
  { ncf[i] = hj*cf[(i-1)*d+1]; hj *= hs;
  }

/*   for i=0..n....
 *     First we compute prod_wk[j], j=0..d.
 *     For j=0, this is int_0^T_i (u-t)^k W((u-t)/h) exp(b0*(u-t)) du
 *     For remaining j,   (x(i,j)-x(j))^k Wj exp(bj*(x..-x.))
 *
 *     Second, we add to the integration (exp(a) incl. in integral)
 *     with the right factorial denominators.
 */
  t_prev = ilim[0]; sb = 0.0;
  for (i=0; i<=haz_lfd->n; i++)
  { if (i==haz_lfd->n)
    { dfx[0] = tmax-t[0];
      for (j=1; j<d; j++) dfx[j] = 0.0;
      eb = exp(sb/haz_lfd->n);
    }
    else
    { eb = exp(base(haz_lfd,i)); sb += base(haz_lfd,i);
      for (j=0; j<d; j++) dfx[j] = datum(haz_lfd,j,i)-t[j];
    }

    if (dfx[0]>ilim[0]) /* else it doesn't contribute */
    {
/* time integral */
      il1 = (dfx[0]>ilim[d]) ? ilim[d] : dfx[0];
      if (il1 != t_prev) /* don't repeat! */
      { st = onedint(haz_sp,ncf,ilim[0]/hs,il1/hs,prod_wk[0]);
        if (st>0) return(st);
        hj = eb;
        for (j=0; j<=2*deg(haz_sp); j++)
        { hj *= hs;
          prod_wk[0][j] *= hj;
        }
        t_prev = il1;
      }

/* covariate terms */
      for (j=1; j<d; j++)
      {
        ef = 0.0;
        for (k=deg(haz_sp); k>0; k--) ef = (ef+dfx[j])*cf[1+(k-1)*d+j];
        ef = exp(ef);
        prod_wk[j][0] = ef * W(dfx[j]/(h*haz_lfd->sca[j]),ker(haz_sp));
        for (k=1; k<=2*deg(haz_sp); k++)
          prod_wk[j][k] = prod_wk[j][k-1] * dfx[j];
      }

/*  add to the integration.  */
      prodintresp(resp,prod_wk,d,deg(haz_sp),p);
    } /* if dfx0 > ilim0 */
  } /* n loop */

/* symmetrize */
  for (k=0; k<p; k++)
    for (j=k; j<p; j++)
      resp[j*p+k] = resp[k*p+j];
  return(LF_OK);
}

int hazint(t,resp,resp1,cf,h)
double *t, *resp, *resp1, *cf, h;
{ if (haz_lfd->d==1) return(hazint_prod(t,resp,resp1,cf,h));
  if (kt(haz_sp)==KPROD) return(hazint_prod(t,resp,resp1,cf,h));

  return(hazint_sph(t,resp,resp1,cf,h));
}

void haz_init(lfd,des,sp,il)
lfdata *lfd;
design *des;
smpar *sp;
double *il;
{ int i;
  
  haz_lfd = lfd;
  haz_sp  = sp;

  tmax = datum(lfd,0,0);
  for (i=1; i<lfd->n; i++) tmax = MAX(tmax,datum(lfd,0,i));
  ff = des->xtwx.wk;
  for (i=0; i<2*lfd->d; i++) ilim[i] = il[i];
}
