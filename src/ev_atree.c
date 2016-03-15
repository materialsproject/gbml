/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *   This file contains functions for constructing and
 *   interpolating the adaptive tree structure. This is
 *   the default evaluation structure used by Locfit.
 */

#include "local.h"

/*
  Guess the number of fitting points.
  Needs improving!
*/
void atree_guessnv(evs,nvm,ncm,vc,d,alp)
evstruc *evs;
double alp;
int *nvm, *ncm, *vc, d;
{ double a0, cu, ifl;
  int i, nv, nc;

  *ncm = 1<<30; *nvm = 1<<30;
  *vc = 1 << d;

  if (alp>0)
  { a0 = (alp > 1) ? 1 : 1/alp;
    if (cut(evs)<0.01)
    { WARN(("guessnv: cut too small."));
      cut(evs) = 0.01;
    }
    cu = 1;
    for (i=0; i<d; i++) cu *= MIN(1.0,cut(evs));
    nv = (int)((5*a0/cu+1)**vc);  /* this allows 10*a0/cu splits */
    nc = (int)(10*a0/cu+1);      /* and 10*a0/cu cells */
    if (nv<*nvm) *nvm = nv;
    if (nc<*ncm) *ncm = nc;
  }

  if (*nvm == 1<<30) /* by default, allow 100 splits */
  { *nvm = 102**vc;
    *ncm = 201;
  }

  /* inflation based on mk */
  ifl = mk(evs)/100.0;
  *nvm = (int)(ifl**nvm);
  *ncm = (int)(ifl**ncm);
  
}

/*
  Determine whether a cell in the tree needs splitting.
  If so, return the split variable (0..d-1).
  Otherwise, return -1.
*/
int atree_split(lf,ce,le,ll,ur)
lfit *lf;
Sint *ce;
double *le, *ll, *ur;
{ int d, vc, i, is;
  double h, hmin, score[MXDIM];
  d = lf->fp.d; vc = 1<<d;

  hmin = 0.0;
  for (i=0; i<vc; i++)
  { h = lf->fp.h[ce[i]];
    if ((h>0) && ((hmin==0)|(h<hmin))) hmin = h;
  }

  is = 0;
  for (i=0; i<d; i++)
  { le[i] = (ur[i]-ll[i])/lf->lfd.sca[i];
    if ((lf->lfd.sty[i]==STCPAR) || (hmin==0))
      score[i] = 2*(ur[i]-ll[i])/(lf->evs.fl[i+d]-lf->evs.fl[i]);
    else
      score[i] = le[i]/hmin;
    if (score[i]>score[is]) is = i;
  }
  if (cut(&lf->evs)<score[is]) return(is);
  return(-1);
}

/*
  recursively grow the tree structure, begining with the parent cell.
*/
void atree_grow(des,lf,ce,ct,term,ll,ur)
design *des;
lfit *lf;
Sint *ce, *ct, *term;
double *ll, *ur;
{ Sint nce[1<<MXDIM];
  int i, i0, i1, d, vc, pv, tk, ns;
  double le[MXDIM], z;
  d = lf->fp.d; vc = 1<<d;

  /* does this cell need splitting?
     If not, wrap up and return.
  */
  ns = atree_split(lf,ce,le,ll,ur);
  if (ns==-1)
  { if (ct != NULL) /* reconstructing terminal cells */
    { for (i=0; i<vc; i++) term[*ct*vc+i] = ce[i];
      (*ct)++;
    }
    return;
  }

  /* split the cell at the midpoint on side ns */
  tk = 1<<ns;
  for (i=0; i<vc; i++)
  { if ((i&tk)==0) nce[i] = ce[i];
    else
    { i0 = ce[i];
      i1 = ce[i-tk];
      pv = (lf->lfd.sty[i]!=STCPAR) &&
           (le[ns] < (cut(&lf->evs)*MIN(lf->fp.h[i0],lf->fp.h[i1])));
      nce[i] = newsplit(des,lf,i0,i1,pv);
      if (lf_error) return;
    }
  }
  z = ur[ns]; ur[ns] = (z+ll[ns])/2;
  atree_grow(des,lf,nce,ct,term,ll,ur);
  if (lf_error) return;
  ur[ns] = z;
  for (i=0; i<vc; i++)
    nce[i] = ((i&tk)== 0) ? nce[i+tk] : ce[i];
  z = ll[ns]; ll[ns] = (z+ur[ns])/2;
  atree_grow(des,lf,nce,ct,term,ll,ur);
  ll[ns] = z;
}

void atree_start(des,lf)
design *des;
lfit *lf;
{ int d, i, j, k, vc, ncm, nvm;
  double ll[MXDIM], ur[MXDIM];

  if (lf_debug>1) printf(" In atree_start\n");
  d = lf->fp.d;
  atree_guessnv(&lf->evs,&nvm,&ncm,&vc,d,nn(&lf->sp));
  if (lf_debug>2) printf(" atree_start: nvm %d ncm %d\n",nvm,ncm);
  trchck(lf,nvm,ncm,vc);

  /* Set the lower left, upper right limits. */
  for (j=0; j<d; j++)
  { ll[j] = lf->evs.fl[j];
    ur[j] = lf->evs.fl[j+d];
  }

  /* Set the initial cell; fit at the vertices. */
  for (i=0; i<vc; i++)
  { j = i;
    for (k=0; k<d; ++k)
    { evptx(&lf->fp,i,k) = (j%2) ? ur[k] : ll[k];
      j >>= 1;
    }
    lf->evs.ce[i] = i;
    des->vfun(des,lf,i);
    if (lf_error) return;
    lf->evs.s[i] = 0;
  }
  lf->fp.nv = vc;

  /* build the tree */
  atree_grow(des,lf,lf->evs.ce,NULL,NULL,ll,ur);
  lf->evs.nce = 1;
}

double atree_int(lf,x,what)
lfit *lf;
double *x;
int what;
{ double vv[64][64], *ll, *ur, h, xx[MXDIM];
  int lo, tk, ns, nv, nc=0, d, i, vc;
  Sint ce[64];

fitpt *fp;
evstruc *evs;
fp = &lf->fp;
evs= &lf->evs;

  d = fp->d;
  vc = 1<<d;
  for (i=0; i<vc; i++)
  { setzero(vv[i],vc);
    nc = exvval(fp,vv[i],i,d,what,1);
    ce[i] = evs->ce[i];
  }
  ns = 0;
  while(ns!=-1)
  { ll = evpt(fp,ce[0]); ur = evpt(fp,ce[vc-1]);
    ns = atree_split(lf,ce,xx,ll,ur);
    if (ns!=-1)
    { tk = 1<<ns;
      h = ur[ns]-ll[ns];
      lo = (2*(x[ns]-ll[ns])) < h;
      for (i=0; i<vc; i++) if ((tk&i)==0)
      { nv = findpt(fp,evs,(int)ce[i],(int)ce[i+tk]);
        if (nv==-1) ERROR(("Descend tree problem"));
        if (lf_error) return(0.0);
        if (lo)
        { ce[i+tk] = nv;
          if (evs->s[nv]) exvvalpv(vv[i+tk],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(fp,vv[i+tk],nv,d,what,1);
        }
        else
        { ce[i] = nv;
          if (evs->s[nv]) exvvalpv(vv[i],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(fp,vv[i],nv,d,what,1);
      } }
  } }
  ll = evpt(fp,ce[0]); ur = evpt(fp,ce[vc-1]);
  return(rectcell_interp(x,vv,ll,ur,d,nc));
}
