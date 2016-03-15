/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

/*
  Functions for computing residuals and fitted values from
  the locfit object.

  fitted(lf,fit,what,cv,ty) computes fitted values from the
    fit structure in lf. 
  resid(y,c,w,th,mi,ty) converts fitted values to residuals
*/

#include "local.h"

double resid(y,w,th,fam,ty,res)
int fam, ty;
double y, w, th, *res;
{ double raw;

  fam = fam & 63;
  if ((fam==TGAUS) | (fam==TROBT) | (fam==TCAUC))
    raw = y-res[ZMEAN];
  else
    raw = y-w*res[ZMEAN];
  switch(ty)
  { case RDEV:
      if (res[ZDLL]>0) return(sqrt(-2*res[ZLIK]));
            else return(-sqrt(-2*res[ZLIK]));
    case RPEAR:
      if (res[ZDDLL]<=0)
      { if (res[ZDLL]==0) return(0);
        return(NOSLN);
      }
      return(res[ZDLL]/sqrt(res[ZDDLL]));
    case RRAW:  return(raw);
    case RLDOT: return(res[ZDLL]);
    case RDEV2: return(-2*res[ZLIK]);
    case RLDDT: return(res[ZDDLL]);
    case RFIT:  return(th);
    case RMEAN: return(res[ZMEAN]);
    default: ERROR(("resid: unknown residual type %d",ty));
  }
  return(0.0);
}

double studentize(res,inl,var,ty,link)
double res, inl, var, *link;
int ty;
{ double den;
  inl *= link[ZDDLL];
  var = var*var*link[ZDDLL];
  if (inl>1) inl = 1;
  if (var>inl) var = inl;
  den = 1-2*inl+var;
  if (den<0) return(0.0);
  switch(ty)
  { case RDEV:
    case RPEAR:
    case RRAW:
    case RLDOT:
      return(res/sqrt(den));
    case RDEV2:
      return(res/den);
    default: return(res);
  }
}

void fitted(lf,fit,what,cv,st,ty)
lfit *lf;
double *fit;
int what, cv, st, ty;
{ int i, j, d, n, evo;
  double xx[MXDIM], th, inl=0.0, var, link[LLEN];
  n = lf->lfd.n;
  d = lf->lfd.d;
  evo = ev(&lf->evs);
  cv &= (evo!=ECROS);
  if ((evo==EDATA)|(evo==ECROS)) evo = EFITP;
  for (i=0; i<n; i++)
  { for (j=0; j<d; j++) xx[j] = datum(&lf->lfd,j,i);
    th = dointpoint(lf,xx,what,evo,i);
    if ((what==PT0)|(what==PVARI)) th = th*th;
    if (what==PCOEF)
    { th += base(&lf->lfd,i);
      stdlinks(link,&lf->lfd,&lf->sp,i,th,rsc(&lf->fp));
      if ((cv)|(st))
      { inl = dointpoint(lf,xx,PT0,evo,i);
        inl = inl*inl;
        if (cv)
        { th -= inl*link[ZDLL];
          stdlinks(link,&lf->lfd,&lf->sp,i,th,rsc(&lf->fp));
        }
        if (st) var = dointpoint(lf,xx,PNLX,evo,i);
      }
      fit[i] = resid(resp(&lf->lfd,i),prwt(&lf->lfd,i),th,fam(&lf->sp),ty,link);
      if (st) fit[i] = studentize(fit[i],inl,var,ty,link);
    } else fit[i] = th;
    if (lf_error) return;
  }
}
