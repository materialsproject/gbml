/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
  startlf(des,lf,vfun,nopc) -- starting point for locfit.
    des and lf are pointers to the design and fit structures.
    vfun is the vertex processing function.
    nopc=1 inhibits computation of parametric component.
  lfit_init(lf) -- initialize the lfit structure.
    lf is pointer to fit.
  preproc() -- fit preprocessing (limits, scales, paramcomp etc.)
  set_scales()
  set_flim()    -- compute bounding box.

  fitoptions()
  clocfit()  -- start point for CLocfit - interpret cmd line etc.
 */

#include "local.h"

void evstruc_init(evs)
evstruc *evs;
{ int i;
  ev(evs) = ETREE;
  mk(evs) = 100;
  cut(evs) = 0.8;
  for (i=0; i<MXDIM; i++)
  { evs->fl[i] = evs->fl[i+MXDIM] = 0.0;
    evs->mg[i] = 10;
  }
  evs->nce = evs->ncm = 0;
}

void fitpt_init(fp)
fitpt *fp;
{ 
  dc(fp) = 0;
  geth(fp) = GSTD;
  fp->nv = fp->nvm = 0;
}

void lfit_init(lf)
lfit *lf;
{
  lfdata_init(&lf->lfd);
  evstruc_init(&lf->evs);
  smpar_init(&lf->sp,&lf->lfd);
  deriv_init(&lf->dv);
  fitpt_init(&lf->fp);
}

void fitdefault(lf)
lfit *lf;
{ WARN(("fitdefault deprecated -- use lfit_init()"));
  lfit_init(lf);
}

void set_flim(lfd,evs)
lfdata *lfd;
evstruc *evs;
{ int i, j, d, n;
  double z, mx, mn, *bx;

  if (ev(evs)==ESPHR) return;
  d = lfd->d; n = lfd->n;
  bx = evs->fl;
  for (i=0; i<d; i++)
    if (bx[i]==bx[i+d])
    { if (lfd->sty[i]==STANGL)
      { bx[i] = 0.0; bx[i+d] = 2*PI*lfd->sca[i];
      }
      else
      { mx = mn = datum(lfd,i,0);
        for (j=1; j<n; j++)
        { mx = MAX(mx,datum(lfd,i,j));
          mn = MIN(mn,datum(lfd,i,j));
        }
        if (lfd->xl[i]<lfd->xl[i+d]) /* user set xlim; maybe use them. */
        { z = mx-mn;
          if (mn-0.2*z < lfd->xl[i]) mn = lfd->xl[i];
          if (mx+0.2*z > lfd->xl[i+d]) mx = lfd->xl[i+d];
        }
        bx[i] = mn;
        bx[i+d] = mx;
      }
    }
}

double vecsum(v,n)
double *v;
int n;
{ int i;
  double sum;
  sum = 0.0;
  for (i=0; i<n; i++) sum += v[i];
  return(sum);
}
 
double vvari(v,n)
double *v;
int n;
{ int i;
  double xb, s2;
  xb = s2 = 0.0;
  for (i=0; i<n; i++) xb += v[i];
  xb /= n;
  for (i=0; i<n; i++) s2 += SQR(v[i]-xb);
  return(s2/(n-1));
}

void set_scales(lfd)
lfdata *lfd;
{ int i;
  for (i=0; i<lfd->d; i++)
    if (lfd->sca[i]<=0) /* set automatic scales */
    { if (lfd->sty[i]==STANGL)
        lfd->sca[i] = 1.0;
      else lfd->sca[i] = sqrt(vvari(lfd->x[i],lfd->n));
    }
}

void startlf(des,lf,vfun,nopc)
design *des;
lfit *lf;
int (*vfun)(), nopc;
{ int i, d, n;

  if (lf_debug>0) printf("startlf\n");
  n = lf->lfd.n;
  d = lf->lfd.d;
  des->vfun = vfun;
  npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);

  des_init(des,n,npar(&lf->sp));
  des->smwt = (lf->lfd.w==NULL) ? n : vecsum(lf->lfd.w,n);
  set_scales(&lf->lfd);
  set_flim(&lf->lfd,&lf->evs);
  compparcomp(des,&lf->lfd,&lf->sp,&lf->pc,geth(&lf->fp),nopc);
  makecfn(&lf->sp,des,&lf->dv,lf->lfd.d);

  lf->lfd.ord = 0;
  if ((d==1) && (lf->lfd.sty[0]!=STANGL))
  { i = 1;
    while ((i<n) && (datum(&lf->lfd,0,i)>=datum(&lf->lfd,0,i-1))) i++;
    lf->lfd.ord = (i==n);
  }
  for (i=0; i<npar(&lf->sp); i++) des->fix[i] = 0;

  lf->fp.d = lf->lfd.d;
  lf->fp.hasd = (des->ncoef==(1+lf->fp.d));

  if (lf_debug>1) printf("call eval structure\n");
  switch(ev(&lf->evs))
  { case EPHULL: triang_start(des,lf); break;
    case EDATA:  dataf(des,lf); break;
    case ECROS:  crossf(des,lf); break;
    case EGRID:  gridf(des,lf); break;
    case ETREE:  atree_start(des,lf); break;
    case EKDCE:  kt(&lf->sp) = KCE;
    case EKDTR:  kdtre_start(des,lf); break;
    case EPRES:  preset(des,lf); break;
    case EXBAR:  xbarf(des,lf); break;
    case ENONE:  lf->fp.nv = lf->evs.nce = 0;
                 return;
    case ESPHR:  sphere_start(des,lf); break;
    case ESPEC:  lf->evs.espec(des,lf); break;
    default: ERROR(("startlf: Invalid evaluation structure %d",ev(&lf->evs)));
  }

  /* renormalize for family=density */
  if ((de_renorm) && (fam(&lf->sp)==TDEN)) dens_renorm(lf,des);
}
