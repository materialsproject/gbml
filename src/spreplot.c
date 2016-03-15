/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include <stdio.h>
#include "local.h"
#include "spreplot.h"

extern int deitype(char *);  /* in lfstr.c */
static lfit lf;
int lf_error;

void vbasis(){;}  // RJN

static void setevs(evs,mi,cut,mg,flim)
evstruc *evs;
int *mg;
Sint *mi;
double cut, *flim;
{ double *ll, *ur;
  int i, d;

  ev(evs) = mi[MEV];
  mk(evs) = mi[MK];
  d = mi[MDIM];

  if (flim != NULL)
  { ll = flim;
    ur = &flim[d];
    memmove(evs->fl,ll,d*sizeof(double));
    memmove(&evs->fl[d],ur,d*sizeof(double));
  }

  switch(ev(evs))
  { case ETREE:
    case EKDTR:
    case EKDCE:
    case EPHULL:
      cut(evs) = cut;
      return;
    case EGRID:
      for (i=0; i<d; i++)
        evs->mg[i] = mg[i];
      return;
    case ESPHR:
      for (i=0; i<2; i++) evs->mg[i] = mg[i];
      return;
    case EDATA:
    case ECROS:
    case EPRES:
    case EXBAR:
    case ENONE:
      return;
    default:
      printf("setevs: %2d not defined.\n",ev(evs));
  }
}

static void setsmpar(sp,dp,mi)
smpar *sp;
double *dp;
Sint *mi;
{ nn(sp)  = dp[DALP];
  fixh(sp)= dp[DFXH];
  pen(sp) = dp[DADP];
  ker(sp) = mi[MKER];
  kt(sp)  = mi[MKT];
  acri(sp)= mi[MACRI];
  deg(sp) = mi[MDEG];
  deg0(sp) = mi[MDEG0];
  fam(sp)  = mi[MTG];
  link(sp) = mi[MLINK];
  ubas(sp) = mi[MUBAS];
  npar(sp) = mi[MP];
  lf.sp.vbasis = vbasis;
}

static void recoef(xev,coef,cell,nvc,mi,dp)
double *xev, *coef, *dp;
Sint *cell, *nvc, *mi;
{ int d, vc=0;

  d = mi[MDIM];
  lf.fp.nv = lf.fp.nvm = nvc[3];
  lf.fp.xev = xev;
  lf.fp.d   = d;
  lf.fp.coef = coef; coef += lf.fp.nv*(d+1);
  lf.fp.nlx  = coef; coef += lf.fp.nv*(d+1);
  lf.fp.t0   = coef; coef += lf.fp.nv*(d+1);
  lf.fp.lik  = coef; coef += lf.fp.nv*3;
  lf.fp.h    = coef; coef += lf.fp.nv;
  lf.fp.deg  = coef; coef += lf.fp.nv;
  rv(&lf.fp) = dp[DRV];
  rsc(&lf.fp)= dp[DRSC];
  dc(&lf.fp) = mi[MDC];
  lf.fp.hasd = (mi[MDEG]>0) | dc(&lf.fp);

  switch(mi[MEV])
  { case ETREE:
    case EKDTR:
    case EGRID:
    case ESPHR: vc = 1<<d; break;
    case EPHULL: vc = d+1; break;
    case EXBAR:
    case ECROS:
    case EDATA:
    case EPRES:
    case ENONE: vc=0; break;
    default: ERROR(("spreplot: Invalid ev"));
  }

  lf.evs.ce = cell; cell += nvc[4]*vc;
  lf.evs.s  = cell; cell += MAX(nvc[3],nvc[4]);
  lf.evs.lo = cell; cell += MAX(nvc[3],nvc[4]);
  lf.evs.hi = cell; cell += MAX(nvc[3],nvc[4]);
}

void spreplot(xev,coef,sv,cell,x,res,se,wpc,sca,m,nvc,mi,dp,mg,dv,nd,sty,where,what,bs)
double *xev, *coef, *sv, *x, *res, *se, *wpc, *sca, *dp;
Sint *cell, *m, *nvc, *mi, *mg, *dv, *nd, *sty, *where;
char **what;
void **bs;
{ Sint i, p;
  double *xx[MXDIM];

  for (i=0; i<mi[MDIM]; i++)
  { lf.lfd.sty[i] = sty[i];
    lf.lfd.sca[i] = sca[i];
  }
  lf.lfd.d = mi[MDIM];

  setsmpar(&lf.sp,dp,mi);
  setevs(&lf.evs,mi,dp[DCUT],mg,NULL);

//if (mi[MUBAS]) bsfunc = bs[0];

  lf_error = 0; p = mi[MP];
  lf.evs.ncm = nvc[1]; lf.evs.nce = nvc[4];

  recoef(xev,coef,cell,nvc,mi,dp);
  lf.evs.sv = sv;

  lf.pc.wk = wpc;
  lf.pc.lwk = pc_reqd(mi[MDIM],p);
  pcchk(&lf.pc,(int)mi[MDIM],p,0);
  haspc(&lf.pc) = mi[MPC];
  lf.pc.xtwx.st = JAC_EIGD;

  /* set up the data structures right */
  switch (*where)
  { case 2: /* grid */
      for (i=0; i<mi[MDIM]; i++)
      { xx[i] = x;
        x += m[i];
      }
      break;
    case 1: /* vector */
    case 3: /* data */
      for (i=0; i<mi[MDIM]; i++) dvari(&(lf.lfd),i) = xx[i] = &x[i * *m];
      break;
    case 4: /* fit points, need nothing! */
      break;
    default:
      ERROR(("unknown where in spreplot"));
  }

  lf.dv.nd = *nd;
  for (i=0; i<*nd; i++)
    lf.dv.deriv[i] = dv[i]-1;

  if (lf_error) return;
  preplot(&lf,xx,res,se,what[1][0],m,*where,ppwhat(what[0]));
}

static void sfitted(x,y,w,c,ba,fit,cv,st,xev,coef,sv,cell,wpc,sca,nvc,mi,dp,mg,dv,nd,sty,what,bs)
double *x, *y, *w, *c, *ba, *fit, *xev, *coef, *sv, *wpc, *sca, *dp;
Sint *cv, *st, *cell, *nvc, *mi, *mg, *dv, *nd, *sty;
char **what;
void **bs;
{ Sint i, n;

  recoef(xev,coef,cell,nvc,mi,dp);
  setsmpar(&lf.sp,dp,mi);
  setevs(&lf.evs,mi,dp[DCUT],mg,NULL);

//if (mi[MUBAS]) bsfunc = bs[0];
  n = mi[MN];
  lf_error = 0;
  lf.evs.ncm = nvc[1]; lf.evs.nce = nvc[4];

  setdata(&lf.lfd,x,y,c,w,ba,mi[MN],mi[MDIM],sca,sty);

  lf.evs.sv = sv;

  lf.pc.wk = wpc;
  lf.pc.lwk= pc_reqd(mi[MDIM],mi[MP],0);
  pcchk(&lf.pc,mi[MDIM],mi[MP],0);
  haspc(&lf.pc) = mi[MPC];
  lf.pc.xtwx.st = JAC_EIGD;

  lf.dv.nd = *nd;
  for (i=0; i<*nd; i++) lf.dv.deriv[i] = dv[i]-1;

  fitted(&lf,fit,ppwhat(what[0]),*cv,*st,restyp(what[1]));
}

static void triterm(xev,h,ce,lo,hi,sca,nvc,mi,dp,nt,term,box)
double *xev, *h, *sca, *dp, *box;
Sint *ce, *lo, *hi, *nvc, *mi, *nt, *term;
{ int i, d, vc;
  Sint mg;

  lf_error = 0;
  d = mi[MDIM];

  lf.fp.xev = xev;
  lf.fp.h = h;
  lf.fp.d = d;
  lf.fp.nv = lf.fp.nvm = nvc[3];

  memmove(lf.lfd.sca,sca,d*sizeof(double));
  setevs(&lf.evs,mi,dp[DCUT],&mg,NULL);

  lf.evs.ce = ce;
  lf.evs.lo = lo;
  lf.evs.hi = hi;
  *nt = 0;

  if (mi[MEV]==ETREE)
    atree_grow(NULL, &lf, lf.evs.ce, nt, term, box, &box[d]);
  else
  { vc = d+1;
    for (i=0; i<nvc[4]; i++)
      triang_grow(NULL,&lf,&lf.evs.ce[i*vc],nt,term);
  }
}

void guessnv(lw,evt,dp,mi,nvc,mg)
double *dp;
char **evt;
int *lw, *mi, *nvc, *mg;
{ int n, d, i, nvm, ncm, vc;
  smpar sp;
  evstruc evs;

  mi[MEV] = ev(&evs) = lfevstr(evt[0]);
  mi[MKT] = lfketype(evt[1]);
  mk(&evs) = mi[MK];
  d = mi[MDIM];
  n = mi[MN];

  switch(mi[MEV])
  { case ETREE:
      cut(&evs) = dp[DCUT];
      atree_guessnv(&evs,&nvm,&ncm,&vc,d,dp[DALP]);
      break;
    case EPHULL:
      nvm = ncm = mi[MK]*mi[MDIM];
      vc = mi[MDIM]+1;
      break;
    case EDATA:
    case ECROS:
      nvm = mi[MN];
      ncm = vc = 0;
      break;
    case EKDTR:
    case EKDCE:
      cut(&evs) = dp[DCUT];
      kdtre_guessnv(&evs,&nvm,&ncm,&vc,n,d,dp[DALP]);
      break;
    case EGRID:
      nvm = 1;
      for (i=0; i<d; i++) nvm *= mg[i];
      ncm = 0;
      vc = 1<<d;
      break;
    case EXBAR:
    case ENONE:
      nvm = 1;
      ncm = vc = 0;
      break;
    case EPRES:
      nvm = mg[0];
      ncm = vc = 0;
      break;
    default:
      ERROR(("guessnv: I don't know this evaluation structure."));
  }

  ubas(&sp)= mi[MUBAS];
  deg(&sp) = mi[MDEG];
  kt(&sp)  = mi[MKT];
  npar(&sp)= mi[MDEG]; /* for user basis */
  mi[MP] = calcp(&sp,d);
  lw[0] = des_reqd(n,(int)mi[MP]);
  lw[1] = lfit_reqd(d,nvm,ncm,(int)mi[MGETH]);
  lw[2] = lfit_reqi(nvm,ncm,vc);
  lw[6] = des_reqi(n,(int)mi[MP]);
  lw[3] = pc_reqd(d,(int)mi[MP]);
  lw[5] = 1;

  if (mi[MGETH] >= 70)
  { lw[4] = k0_reqd(d,n,0);
    if (lw[4]<2*nvm) lw[4] = 2*nvm;
    lw[5] = d+1;
  }
  else switch(mi[MGETH])
  { case GSTD:  lw[4] = 1; break;                  /* standard fit */
    case GSMP:  lw[4] = 1; break;                  /* simple fit   */
    case GHAT:  lw[4] = nvm*n; break;              /* hat matrix   */
    case GKAP:  lw[4] = k0_reqd(d,n,0);            /* kappa0       */
             lw[5] = 1+d;
             break;
    case GRBD:  lw[5] = 10;                        /* regband      */
    case GAMF:                                     /* gam.lf fit   */
    case GAMP:  lw[4] = 1; break;                  /* gam.lf pred  */
    case GLSC:  lw[4] = 2; break;                  /* lscv         */
    default:
      printf("sguessnv: invalid geth\n");
      lw[4] = 0;
  }

  nvc[0] = nvm;
  nvc[1] = ncm;
  nvc[2] = vc;
  nvc[3] = nvc[4] = 0;

}

/* Registration added Mar 2012 */

/* From smisc.c */
void kdeb(double *x, int *mi, double*band, int *ind, double *h0, double *h1,
	  int *meth, int *nmeth, int *ker);
void scritval(double *k0, int *d, double *cov, int *m, double *rdf, 
	      double *z, int *k);
void slscv(double *x, int *n, double *h, double *z);
