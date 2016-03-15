/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 * 
 *
 * functions for computing and subtracting, adding the
 * parametric component
 */

#include "local.h"

int noparcomp(sp,geth)
smpar *sp;
int geth;
{ int tg;
  if (geth==GSMP) return(1);
  if (deg0(sp)<deg(sp)) return(1);
  if (ubas(sp)) return(1);
  tg = fam(sp) & 63;
  if (tg<=THAZ) return(1);
  if (tg==TROBT) return(1);
  if (tg==TCAUC) return(1);
  return(0);
}

int pc_reqd(d,p)
int d, p;
{ return(d + 2*p + jac_reqd(p));
}

void pcchk(pc,d,p,lc)
paramcomp *pc;
int d, p, lc;
{ int rw;
  double *z;

  rw = pc_reqd(d,p);
  if (pc->lwk < rw)
  { pc->wk = (double *)calloc(rw,sizeof(double));
    pc->lwk= rw;
  }
  z = pc->wk;

  pc->xbar = z; z += d;
  pc->coef = z; z += p;
  pc->f    = z; z += p;

  z = jac_alloc(&pc->xtwx,p,z);
  pc->xtwx.p = p;
}

void compparcomp(des,lfd,sp,pc,geth,nopc)
design *des;
lfdata *lfd;
smpar *sp;
paramcomp *pc;
int geth;
int nopc;
{ int i, j, k, p;
  double wt, sw;

  if (lf_debug>1) printf(" compparcomp:\n");
  p = des->p;
  pcchk(pc,lfd->d,p,1);

  for (i=0; i<lfd->d; i++) pc->xbar[i] = 0.0;
  sw = 0.0;
  for (i=0; i<lfd->n; i++)
  { 
    wt = prwt(lfd,i);
    sw += wt;
    for (j=0; j<lfd->d; j++)
      pc->xbar[j] += datum(lfd,j,i)*wt;
    des->ind[i] = i;
    des->w[i] = 1.0;
  }
  for (i=0; i<lfd->d; i++) pc->xbar[i] /= sw;
  if ((nopc) || noparcomp(sp,geth))
  { haspc(pc) = 0;
    return;
  }
  haspc(pc) = 1;
  des->xev = pc->xbar;
  k = locfit(lfd,des,sp,0,0,0);
  if (lf_error) return;
  switch(k)
  { case LF_NOPT:
      ERROR(("compparcomp: no points in dataset?"));
      return;
    case LF_INFA:
      ERROR(("compparcomp: infinite parameters in param. component"));
      return;
    case LF_NCON:
      ERROR(("compparcom: not converged"));
      return;
    case LF_OOB:
      ERROR(("compparcomp: parameters out of bounds"));
      return;
    case LF_PF:
      WARN(("compparcomp: perfect fit"));
    case LF_OK:
      for (i=0; i<p; i++)
      { pc->coef[i] = des->cf[i];
        pc->xtwx.dg[i] = des->xtwx.dg[i];
        pc->xtwx.wk[i] = des->xtwx.wk[i];
      }
      for (i=0; i<p*p; i++)
      { pc->xtwx.Z[i] = des->xtwx.Z[i];
        pc->xtwx.Q[i] = des->xtwx.Q[i];
      }
      pc->xtwx.sm = des->xtwx.sm;
      pc->xtwx.st = des->xtwx.st;
      return;
    default:
      ERROR(("compparcomp: locfit unknown return status %d",k));
      return;
  }
}

void subparcomp(des,lf,coef)
design *des;
lfit *lf;
double *coef;
{ int i, nd;
  deriv *dv;
  paramcomp *pc;

  pc = &lf->pc;
  if (!haspc(pc)) return;

  dv = &lf->dv; nd = dv->nd;
  fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,des->f1,dv);
  coef[0] -= innerprod(pc->coef,des->f1,pc->xtwx.p);
  if (des->ncoef == 1) return;

  dv->nd = nd+1;
  for (i=0; i<lf->lfd.d; i++)
  { dv->deriv[nd] = i;
    fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,des->f1,dv);
    coef[i+1] -= innerprod(pc->coef,des->f1,pc->xtwx.p);
  }
  dv->nd = nd;
}

void subparcomp2(des,lf,vr,il)
design *des;
lfit *lf;
double *vr, *il;
{ double t0, t1;
  int i, nd;
  deriv *dv;
  paramcomp *pc;

  pc = &lf->pc;
  if (!haspc(pc)) return;

  dv = &lf->dv; nd = dv->nd;

  fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,des->f1,dv);
  for (i=0; i<npar(&lf->sp); i++) pc->f[i] = des->f1[i];
  jacob_solve(&pc->xtwx,des->f1);
  t0 = sqrt(innerprod(pc->f,des->f1,pc->xtwx.p));
  vr[0] -= t0;
  il[0] -= t0;
  if ((t0==0) | (des->ncoef==1)) return;

  dv->nd = nd+1;
  for (i=0; i<lf->lfd.d; i++)
  { dv->deriv[nd] = i;
    fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,pc->f,dv);
    t1 = innerprod(pc->f,des->f1,pc->xtwx.p)/t0;
    vr[i+1] -= t1;
    il[i+1] -= t1;
  }
  dv->nd = nd;
}

double addparcomp(lf,x,c)
lfit *lf;
double *x;
int c;
{ double y;
  paramcomp *pc;

  pc = &lf->pc;
  if (!haspc(pc)) return(0.0);
  fitfun(&lf->lfd, &lf->sp, x,pc->xbar,pc->f,&lf->dv);
  if (c==PCOEF) return(innerprod(pc->coef,pc->f,pc->xtwx.p));
  if ((c==PNLX)|(c==PT0)|(c==PVARI))
  { y = sqrt(jacob_qf(&pc->xtwx,pc->f));
    return(y);
  }
  return(0.0);
}
