/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

int lf_maxit = 20;
int lf_debug = 0;

static double s0, s1, tol;
static lfdata *lf_lfd;
static design *lf_des;
static smpar   *lf_sp;
int lf_status;
int ident=0;
int (*like)();
extern double robscale;

void lfdata_init(lfd)
lfdata *lfd;
{ int i;
  for (i=0; i<MXDIM; i++)
  { lfd->sty[i] = 0;
    lfd->sca[i] = 1.0;
    lfd->xl[i] = lfd->xl[i+MXDIM] = 0.0;
  }
  lfd->y = lfd->w = lfd->c = lfd->b = NULL;
  lfd->d = lfd->n = 0;
}

void smpar_init(sp,lfd)
smpar *sp;
lfdata *lfd;
{ nn(sp)  = 0.7;
  fixh(sp)= 0.0;
  pen(sp) = 0.0;
  acri(sp)= ANONE;
  deg(sp) = deg0(sp) = 2;
  ubas(sp) = 0;
  kt(sp) = KSPH;
  ker(sp) = WTCUB;
  fam(sp) = 64+TGAUS;
  link(sp)= LDEFAU;
  npar(sp) = calcp(sp,lfd->d);
}

void deriv_init(dv)
deriv *dv;
{ dv->nd = 0;
}

int des_reqd(n,p)
int n, p;
{
  return(n*(p+5)+2*p*p+4*p + jac_reqd(p));
}
int des_reqi(n,p)
int n, p;
{ return(n+p);
}
 
void des_init(des,n,p)
design *des;
int n, p;
{ double *z;
  int k;

  if (n<=0) WARN(("des_init: n <= 0"));
  if (p<=0) WARN(("des_init: p <= 0"));

  if (des->des_init_id != DES_INIT_ID)
  { des->lwk = des->lind = 0;
    des->des_init_id = DES_INIT_ID;
  }

  k = des_reqd(n,p);
  if (k>des->lwk)
  { des->wk = (double *)calloc(k,sizeof(double));
    des->lwk = k;
  }
  z = des->wk;

  des->X = z; z += n*p;
  des->w = z; z += n;
  des->res=z; z += n;
  des->di =z; z += n;
  des->th =z; z += n;
  des->wd =z; z += n;
  des->V  =z; z += p*p;
  des->P  =z; z += p*p;
  des->f1 =z; z += p;
  des->ss =z; z += p;
  des->oc =z; z += p;
  des->cf =z; z += p;
 
  z = jac_alloc(&des->xtwx,p,z);
 
  k = des_reqi(n,p);
  if (k>des->lind)
  {
    des->ind = (Sint *)calloc(k,sizeof(Sint));
    des->lind = k;
  }
  des->fix = &des->ind[n];
  for (k=0; k<p; k++) des->fix[k] = 0;

  des->n = n; des->p = p;
  des->smwt = n;
  des->xtwx.p = p;                                                              
}

void deschk(des,n,p)
design *des;
int n, p;
{ WARN(("deschk deprecated - use des_init()"));
  des_init(des,n,p);
}

int likereg(coef, lk0, f1, Z)
double *coef, *lk0, *f1, *Z;
{ int i, ii, j, p;
  double lk, ww, link[LLEN], *X;

  if (lf_debug>2) printf("  likereg: %8.5f\n",coef[0]);
  lf_status = LF_OK;
  lk = 0.0; p = lf_des->p;
  setzero(Z,p*p);
  setzero(f1,p);
  for (i=0; i<lf_des->n; i++)
  {
    ii = lf_des->ind[i];
    X = d_xi(lf_des,i);
    lf_des->th[i] = base(lf_lfd,ii)+innerprod(coef,X,p);
    lf_status = stdlinks(link,lf_lfd,lf_sp,ii,lf_des->th[i],robscale);
    if (lf_status == LF_BADP)
    { *lk0 = -1.0e300;
      return(NR_REDUCE);
    }
    if (lf_error) lf_status = LF_ERR;
    if (lf_status != LF_OK) return(NR_BREAK);

    ww = lf_des->w[i];
    lk += ww*link[ZLIK];
    for (j=0; j<p; j++)
      f1[j] += X[j]*ww*link[ZDLL];
    addouter(Z, X, X, p, ww*link[ZDDLL]);
  }
  for (i=0; i<p; i++) if (lf_des->fix[i])
  { for (j=0; j<p; j++) Z[i*p+j] = Z[j*p+i] = 0.0;
    Z[i*p+i] = 1.0;
    f1[i] = 0.0;
  }

  if (lf_debug>4) prresp(coef,Z,p);
  if (lf_debug>3) printf("  likelihood: %8.5f\n",lk);
  *lk0 = lf_des->llk = lk;

  switch (fam(lf_sp)&63) /* parameter checks */
  { case TGAUS: /* prevent iterations! */
      if ((link(lf_sp)==LIDENT)&((fam(lf_sp)&128)==0)) return(NR_BREAK);
        break;
    case TPOIS:
    case TGEOM:
    case TWEIB:
    case TGAMM:
      if ((link(lf_sp)==LLOG) && (fabs(coef[0])>700))
      { lf_status = LF_OOB;
        return(NR_REDUCE);
      }
      if (lk > -1.0e-5*s0)
      { lf_status = LF_PF;
        return(NR_REDUCE);
      }
      break;
    case TRBIN:
    case TLOGT:
      if (lk > -1.0e-5*s0)
      { lf_status = LF_PF;
        return(NR_REDUCE);
      }
      if (fabs(coef[0])>700)
      { lf_status = LF_OOB;
        return(NR_REDUCE);
      }
      break;
  }
  return(NR_OK);
}

int robustinit(lfd,des)
lfdata *lfd;
design *des;
{ int i;
  for (i=0; i<des->n; i++)
    des->res[i] = resp(lfd,(int)des->ind[i]) - base(lfd,(int)des->ind[i]);
  des->cf[0] = median(des->res,des->n);
  for (i=1; i<des->p; i++) des->cf[i] = 0.0;
  tol = 1.0e-6;
  return(LF_OK);
}

int circinit(lfd,des)
lfdata *lfd;
design *des;
{ int i, ii;
  double s0, s1;
  s0 = s1 = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    s0 += des->w[i]*prwt(lfd,ii)*sin(resp(lfd,ii)-base(lfd,ii));
    s1 += des->w[i]*prwt(lfd,ii)*cos(resp(lfd,ii)-base(lfd,ii));
  }
  des->cf[0] = atan2(s0,s1);
  for (i=1; i<des->p; i++) des->cf[i] = 0.0;
  tol = 1.0e-6;
  return(LF_OK);
}

int reginit(lfd,des)
lfdata *lfd;
design *des;
{ int i, ii;
  double sb, link[LLEN];
  s0 = s1 = sb = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    links(base(lfd,ii),resp(lfd,ii),fam(lf_sp),LINIT,link,cens(lfd,ii),prwt(lfd,ii),1.0);
    s1 += des->w[i]*link[ZDLL];
    s0 += des->w[i]*prwt(lfd,ii);
    sb += des->w[i]*prwt(lfd,ii)*base(lfd,ii);
  }
  if (s0==0) return(LF_NOPT); /* no observations with W>0 */
  setzero(des->cf,des->p);
  tol = 1.0e-6*s0;
  switch(link(lf_sp))
  { case LIDENT:
      des->cf[0] = (s1-sb)/s0;
      return(LF_OK);
    case LLOG:
      if (s1<=0.0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      des->cf[0] = log(s1/s0) - sb/s0;
      return(LF_OK);
    case LLOGIT:
      if (s1<=0.0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      if (s1>=s0)
      { des->cf[0] = 1000;
        return(LF_INFA);
      }
      des->cf[0] = logit(s1/s0)-sb/s0;
      return(LF_OK);
    case LINVER:
      if (s1<=0.0)
      { des->cf[0] = 1000;
        return(LF_INFA);
      }
      des->cf[0] = s0/s1-sb/s0;
      return(LF_OK);
    case LSQRT:
      des->cf[0] = sqrt(s1/s0)-sb/s0;
      return(LF_OK);
    case LASIN:
      des->cf[0] = asin(sqrt(s1/s0))-sb/s0;
      return(LF_OK);
    default:
      ERROR(("reginit: invalid link %d",link(lf_sp)));
      return(LF_ERR);
  }
}

int lfinit(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{
  des->xtwx.sm = (deg0(sp)<deg(sp)) ? JAC_CHOL : JAC_EIGD;

  designmatrix(lfd,sp,des);

  like = likereg;
  link(sp) = defaultlink(link(sp),fam(sp));

  switch(fam(sp)&63)
  { case TDEN:
    case TRAT:
    case THAZ:
      like = likeden;
      tol = (link(sp)==LLOG) ? 1.0e-6 : 0.0;
      return(densinit(lfd,des,sp,des->cf));
    case TCAUC:
    case TROBT:
      return(robustinit(lfd,des));
    case TCIRC:
      return(circinit(lfd,des));
    default:
      return(reginit(lfd,des));
  }
}

void lfiter(des,maxit)
design *des;
int maxit;
{ int err;
  if (lf_debug>1) printf(" lfiter: %8.5f\n",des->cf[0]);
  max_nr(like, des->cf, des->oc, des->res, des->f1,
    &des->xtwx, des->p, maxit, tol, &err);
  switch(err)
  { case NR_OK: return;
    case NR_NCON:
      WARN(("max_nr not converged"));
      return;
    case NR_NDIV:
      WARN(("max_nr reduction problem"));
      return;
  }
  WARN(("max_nr return status %d",err));
}

int use_robust_scale(int tg)
{ if ((tg&64)==0) return(0); /* not quasi - no scale */
  if (((tg&128)==0) & (((tg&63)!=TROBT) & ((tg&63)!=TCAUC))) return(0);
  return(1);
}

int locfit(lfd,des,sp,noit,nb,cv)
lfdata *lfd;
design *des;
smpar *sp;
int noit, nb, cv;
{ int i;

  if (des->xev==NULL)
  { ERROR(("locfit: NULL evaluation point?"));
    return(246);
  }

  if (lf_debug>0)
  { printf("locfit: ");
    for (i=0; i<lfd->d; i++) printf(" %10.6f",des->xev[i]);
    printf("\n");
  }

  lf_des = des;
  lf_lfd = lfd;
  lf_sp  = sp;

/* the 1e-12 avoids problems that can occur with roundoff */
  if (nb) nbhd(lfd,des,(int)(lfd->n*nn(sp)+1e-12),0,sp);

  lf_status = lfinit(lfd,sp,des);
  if (lf_status != LF_OK) return(lf_status);

  if (use_robust_scale(fam(sp)))
    lf_robust(lfd,sp,des,lf_maxit);
  else
  { robscale = 1.0;
    lfiter(des,lf_maxit);
  }

  if (lf_status == LF_OOB) setzero(des->cf,des->p);

  if ((fam(sp)&63)==TDEN) /* convert from rate to density */
  { switch(link(sp))
    { case LLOG:
        des->cf[0] -= log(des->smwt);
        break;
      case LIDENT:
        multmatscal(des->cf,1.0/des->smwt,des->p);
        break;
      default: ERROR(("Density adjustment; invalid link"));
    }
  }

  /* variance calculations, if requested */
  if (cv)
    lf_vcov(lfd,sp,des);

  return(lf_status);
}
