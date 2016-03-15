#include "local.h"

static lfit *lf_scb;
static lfdata *lfd_scb;
static smpar  *scb_sp;
static design *des_scb;

int scbfitter(x,l,reqd)
double *x, *l;
int reqd;
{
  int m;
  des_scb->xev = x;
  if ((ker(scb_sp)!=WPARM) | (!haspc(&lf_scb->pc)))
  { locfit(lfd_scb,des_scb,&lf_scb->sp,1,1);
    m = wdiag(lfd_scb, scb_sp, des_scb,l,&lf_scb->dv,reqd,2,0);
  }
  else
    m = wdiagp(lfd_scb, scb_sp, des_scb,l,&lf_scb->pc,&lf_scb->dv,reqd,2,0);
  return(m);
}

/* function to test tube_constants with covariance.
double ll[5000];
int scbfitter2(x,l,reqd)
double *x, *l;
int reqd;
{ double h;
  int d, m, n, i, j;

  m = scbfitter(x,ll,reqd);

  d = lfd_scb->d;

  n = d*d+d+1;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      l[i*n+j] = innerprod(&ll[i*m],&ll[j*m],m);

  return(n);
}
*/

int constants(des,lf)
design *des;
lfit *lf;
{
  int d, m, nt, rw;
  evstruc *evs;

  lf_scb = lf;
  des_scb = des;
  lfd_scb = &lf->lfd;
  scb_sp  = &lf->sp;

  evs = &lf->evs;
  d = lfd_scb->d;
  m = lfd_scb->n;

  if (lf_error) return(0);
  if ((ker(scb_sp) != WPARM) && (lf->sp.nn>0))
    WARN(("constants are approximate for varying h"));
  npar(scb_sp) = calcp(scb_sp,lf->lfd.d);
  des_init(des,m,npar(scb_sp));
  set_scales(&lf->lfd);
  set_flim(&lf->lfd,&lf->evs);
  compparcomp(des,&lf->lfd,&lf->sp,&lf->pc,geth(&lf->fp),ker(scb_sp)!=WPARM);
  
  rw = k0_reqd(d,m,0);
  if (lf->fp.ll<rw)
  { lf->fp.L = (double *)calloc(rw,sizeof(double));
    lf->fp.ll= rw;
  }

  nt = tube_constants(scbfitter,d,m,ev(evs),mg(evs),evs->fl,
    lf->fp.kap,lf->fp.L,(d>3) ? 4 : d+1,0);
  return(nt);
}
