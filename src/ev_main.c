/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

/*
 * trchck checks the working space on the lfit structure 
 * has space for nvm vertices and ncm cells.
 */
void lfit_alloc(lf)
lfit *lf;
{ lf->fp.lwk = lf->fp.lev = lf->fp.ll = lf->evs.liw = lf->pc.lwk = 0;
  lf->lf_init_id = LF_INIT_ID;
}
int lfit_reqd(d,nvm,ncm,geth)
int d, nvm, ncm, geth;
{ int z;
  z = (geth==GSMP) ? d+3 : 3*d+8;
  return(nvm*z+ncm);
}
int lfit_reqi(nvm,ncm,vc)
int nvm, ncm, vc;
{ return(ncm*vc+3*MAX(ncm,nvm));
}

void trchck(lf,nvm,ncm,vc)
lfit *lf;
int nvm, ncm, vc;
{ int rw, d;
  Sint *k;
  double *z;

  if (lf->lf_init_id != LF_INIT_ID) lfit_alloc(lf);

  d = lf->lfd.d;

  if (lf->fp.lev < d*nvm)
  { lf->fp.xev = (double *)calloc(d*nvm,sizeof(double));
    lf->fp.lev = d*nvm;
  }

  rw = lfit_reqd(d,nvm,ncm,geth(&lf->fp));
  if (lf->fp.lwk < rw)
  { lf->fp.coef = (double *)calloc(rw,sizeof(double));
    lf->fp.lwk = rw;
  }
  z = lf->fp.coef;

  lf->fp.coef= z; z += nvm*(d+1);
  if (geth(&lf->fp) != GSMP)
  { lf->fp.nlx = z; z += nvm*(d+1);
    lf->fp.t0  = z; z += nvm*(d+1);
    lf->fp.lik = z; z += 3*nvm;
  }
  lf->fp.h   = z; z += nvm;
  lf->fp.deg = z; z += nvm;
  lf->evs.sv = z; z += ncm;

  rw = lfit_reqi(nvm,ncm,vc);
  if (lf->evs.liw<rw)
  { lf->evs.iwk = (Sint *)calloc(rw,sizeof(Sint));
    lf->evs.liw = rw;
  }
  k = lf->evs.iwk;
  lf->evs.ce = k; k += vc*ncm;
  lf->evs.s  = k; k += MAX(ncm,nvm);
  lf->evs.lo = k; k += MAX(ncm,nvm);
  lf->evs.hi = k; k += MAX(ncm,nvm);

  lf->fp.nvm = nvm; lf->evs.ncm = ncm;
}

void data_guessnv(nvm,ncm,vc,n)
int *nvm, *ncm, *vc, n;
{ *nvm = n;
  *ncm = *vc = 0;
}

void dataf(des,lf)
design *des;
lfit *lf;
{
  int d, i, j, ncm, nv, vc;

  d = lf->lfd.d;
  data_guessnv(&nv,&ncm,&vc,lf->lfd.n);
  trchck(lf,nv,ncm,vc);

  for (i=0; i<nv; i++)
    for (j=0; j<d; j++) evptx(&lf->fp,i,j) = datum(&lf->lfd,j,i);
  for (i=0; i<nv; i++)
  { des->vfun(des,lf,i);
    lf->evs.s[i] = 0;
  }
  lf->fp.nv = lf->fp.nvm = nv; lf->evs.nce = 0;
}

void xbar_guessnv(nvm,ncm,vc)
int *nvm, *ncm, *vc;
{ *nvm = 1;
  *ncm = *vc = 0;
  return;
}

void xbarf(des,lf)
design *des;
lfit *lf;
{ int i, d, nvm, ncm, vc;
  d = lf->lfd.d;
  xbar_guessnv(&nvm,&ncm,&vc);
  trchck(lf,1,0,0);
  for (i=0; i<d; i++) evptx(&lf->fp,0,i) = lf->pc.xbar[i];
  des->vfun(des,lf,0);
  lf->evs.s[0] = 0;
  lf->fp.nv = 1; lf->evs.nce = 0;
}

void preset(des,lf)
design *des;
lfit *lf;
{ int i, nv;

  nv = lf->fp.nvm;
  trchck(lf,nv,0,0);
  for (i=0; i<nv; i++)
  { 
    des->vfun(des,lf,i);
    lf->evs.s[i] = 0;
  }
  lf->fp.nv = nv; lf->evs.nce = 0;
}

void crossf(des,lf)
design *des;
lfit *lf;
{ int d, i, j, n, nv, ncm, vc;
  double w;

  n = lf->lfd.n; d = lf->lfd.d;
  data_guessnv(&nv,&ncm,&vc,n);
  trchck(lf,nv,ncm,vc);

  if (lf->lfd.w==NULL) ERROR(("crossf() needs prior weights"));
  for (i=0; i<n; i++)
    for (j=0; j<d; j++) evptx(&lf->fp,i,j) = datum(&lf->lfd,j,i);
  for (i=0; i<n; i++)
  { lf->evs.s[i] = 0;
    w = prwt(&lf->lfd,i);
    lf->lfd.w[i] = 0;
    des->vfun(des,lf,i);
    lf->lfd.w[i] = w;
  }
  lf->fp.nv = n; lf->evs.nce = 0;
}

void gridf(des,lf)
design *des;
lfit *lf;
{ int d, i, j, nv, u0, u1, z;
  nv = 1; d = lf->lfd.d;
  for (i=0; i<d; i++)
  { if (lf->evs.mg[i]==0)
      lf->evs.mg[i] = 2+(int)((lf->evs.fl[i+d]-lf->evs.fl[i])/(lf->lfd.sca[i]*cut(&lf->evs)));
    nv *= lf->evs.mg[i];
  }
  trchck(lf,nv,0,1<<d);
  for (i=0; i<nv; i++)
  { z = i;
    for (j=0; j<d; j++)
    { u0 = z%lf->evs.mg[j];
      u1 = lf->evs.mg[j]-1-u0;
      evptx(&lf->fp,i,j) = (lf->evs.mg[j]==1) ? lf->evs.fl[j] :
                      (u1*lf->evs.fl[j]+u0*lf->evs.fl[j+d])/(lf->evs.mg[j]-1);
      z = z/lf->evs.mg[j];
    }
    lf->evs.s[i] = 0;
    des->vfun(des,lf,i);
  }
  lf->fp.nv = nv; lf->evs.nce = 0;
}

int findpt(fp,evs,i0,i1)
fitpt *fp;
evstruc *evs;
int i0, i1;
{ int i;
  if (i0>i1) ISWAP(i0,i1);
  for (i=i1+1; i<fp->nv; i++)
    if ((evs->lo[i]==i0) && (evs->hi[i]==i1)) return(i);
  return(-1);
}

/*
  add a new vertex at the midpoint of (x[i0],x[i1]).
  return the vertex number.
*/
int newsplit(des,lf,i0,i1,pv)
design *des;
lfit *lf;
int i0, i1, pv;
{ int i, nv;

  i = findpt(&lf->fp,&lf->evs,i0,i1);
  if (i>=0) return(i);

  if (i0>i1) ISWAP(i0,i1);
  nv = lf->fp.nv;
  
  /* the point is new. Now check we have space for the new point. */
  if (nv==lf->fp.nvm)
  {
    ERROR(("newsplit: out of vertex space"));
    return(-1);
  }

  /* compute the new point, and evaluate the fit */
  lf->evs.lo[nv] = i0;
  lf->evs.hi[nv] = i1;
  for (i=0; i<lf->fp.d; i++)
    evptx(&lf->fp,nv,i) = (evptx(&lf->fp,i0,i)+evptx(&lf->fp,i1,i))/2;
  if (pv) /* pseudo vertex */
  { lf->fp.h[nv] = (lf->fp.h[i0]+lf->fp.h[i1])/2;
    lf->evs.s[nv] = 1; /* pseudo-vertex */
  }
  else /* real vertex */
  {
    des->vfun(des,lf,nv);
    lf->evs.s[nv] = 0;
  }
  lf->fp.nv++;

  return(nv);
}
