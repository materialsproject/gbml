/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern double robscale;

/* special version of ressumm to estimate sigma^2, with derivative estimation */
void ressummd(lf)
lfit *lf;
{ int i;
  double s0, s1;
  s0 = s1 = 0.0;
  if ((fam(&lf->sp)&64)==0)
  { rv(&lf->fp) = 1.0;
    return;
  }
  for (i=0; i<lf->fp.nv; i++)
  { s0 += lf->fp.lik[2*lf->fp.nvm+i];
    s1 += lf->fp.lik[i];
  }
  if (s0==0.0)
    rv(&lf->fp) = 0.0;
  else
    rv(&lf->fp) = -2*s1/s0;
}

void ressumm(lf,des)
lfit *lf;
design *des;
{ int i, j, evo, tg, orth;
  double *oy, pw, r1, r2, rdf, t0, t1, u[MXDIM], link[LLEN];
  fitpt *fp;

  fp = &lf->fp;
  llk(fp) = df0(fp) = df1(fp) = 0.0;

  evo = ev(&lf->evs);
  if ((evo==EKDCE) | (evo==EPRES))
  { rv(fp) = 1.0;
    return;
  }
  if (lf->dv.nd>0)
  { ressummd(lf);
    return;
  }
  r1 = r2 = 0.0;
  if ((evo==EDATA) | (evo==ECROS)) evo = EFITP;
  orth = (geth(&lf->fp)==GAMF) | (geth(&lf->fp)==GAMP);
  for (i=0; i<lf->lfd.n; i++)
  { for (j=0; j<lf->lfd.d; j++) u[j] = datum(&lf->lfd,j,i);
    des->th[i] = base(&lf->lfd,i)+dointpoint(lf,u,PCOEF,evo,i);
    des->wd[i] = resp(&lf->lfd,i) - des->th[i];
    des->w[i] = 1.0;
    des->ind[i] = i;
  }

  tg = fam(&lf->sp);
  rsc(&lf->fp) = 1.0;
  if ((tg==TROBT+64) | (tg==TCAUC+64)) /* global robust scale */
  { oy = lf->lfd.y; lf->lfd.y = des->wd;
    des->xev = lf->pc.xbar;
    locfit(&lf->lfd,des,&lf->sp,1,0);
    lf->lfd.y = oy;
    rsc(fp) = robscale;
  }

  if (orth) /* orthog. residuals */
  { int od, op;
    des->n = lf->lfd.n;
    od = deg(&lf->sp); op = npar(&lf->sp);
    deg(&lf->sp) = 1;
    npar(&lf->sp) = des->p = 1+lf->lfd.d;
    oy = lf->lfd.y; lf->lfd.y = des->wd;
    des->xev = lf->pc.xbar;
    locfit(&lf->lfd,des,&lf->sp,1,0);
    for (i=0; i<lf->lfd.n; i++) oy[i] = resp(&lf->lfd,i) - des->th[i];
    lf->lfd.y = oy;
    deg(&lf->sp) = od; npar(&lf->sp) = op;
  }

  for (i=0; i<lf->lfd.n; i++)
  { for (j=0; j<lf->lfd.d; j++) u[j] = datum(&lf->lfd,j,i);
    t0 = dointpoint(lf,u,PT0,evo,i);
    t1 = dointpoint(lf,u,PNLX,evo,i);
    stdlinks(link,&lf->lfd,&lf->sp,i,des->th[i],rsc(fp));
    t1 = t1*t1*link[ZDDLL];
    t0 = t0*t0*link[ZDDLL];
    if (t1>1) t1 = 1;
    if (t0>1) t0 = 1; /* no observation gives >1 deg.free */
    llk(fp) += link[ZLIK];
    df0(fp) += t0;
    df1(fp) += t1;
    pw = prwt(&lf->lfd,i);
    if (pw>0)
    { r1 += link[ZDLL]*link[ZDLL]/pw;
      r2 += link[ZDDLL]/pw;
    }
    if (orth) des->di[i]  = t1;
  }

  if (orth) return;

  rv(fp) = 1.0;
  if ((fam(&lf->sp)&64)==64) /* quasi family */
  { rdf = lf->lfd.n-2*df0(fp)+df1(fp);
    if (rdf<1.0)
    { WARN(("Estimated rdf < 1.0; not estimating variance"));
    }
    else
      rv(fp) = r1/r2 * lf->lfd.n / rdf;
  }

  /* try to ensure consistency for family="circ"! */
  if (((fam(&lf->sp)&63)==TCIRC) & (lf->lfd.d==1))
  { Sint *ind;
    int nv;
    double dlt, th0, th1;
    ind = des->ind;
    nv = fp->nv;
    for (i=0; i<nv; i++) ind[i] = i;
    lforder(ind,evp(fp),0,nv-1);
    for (i=1; i<nv; i++)
    { dlt = evptx(fp,ind[i],0)-evptx(fp,ind[i-1],0);
      th0 = fp->coef[ind[i]]-dlt*fp->coef[ind[i]+nv]-fp->coef[ind[i-1]];
      th1 = fp->coef[ind[i]]-dlt*fp->coef[ind[i-1]+nv]-fp->coef[ind[i-1]];
      if ((th0>PI)&(th1>PI))
      { for (j=0; j<i; j++)
          fp->coef[ind[j]] += 2*PI;
        i--;
      }
      if ((th0<(-PI))&(th1<(-PI)))
      { for (j=0; j<i; j++)
          fp->coef[ind[j]] -= 2*PI;
        i--;
      }
    }
  }
}

double rss(lf,des,df)
lfit *lf;
design *des;
double *df;
{ double ss;
  ss = 0;
  ressumm(lf,des);
  *df = lf->lfd.n - 2*df0(&lf->fp) + df1(&lf->fp);
  return(-2*llk(&lf->fp));
}
