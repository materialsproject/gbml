/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern double robscale;

double vocri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ if (pen==0) return(-2*t0*lk/((t0-t2)*(t0-t2)));
  return((-2*lk+pen*t2)/t0);
}

int procvraw(des,lf,v)
design *des;
lfit *lf;
int v;
{ int i, lf_status;
  double coef[1+MXDIM];

  if (lf_debug>1) printf(" procvraw: %d\n",v);
  des->xev = evpt(&lf->fp,v);

  if (acri(&lf->sp)==ANONE)
    lf_status = locfit(&lf->lfd,des,&lf->sp,0,1,0);
  else
    lf_status = alocfit(&lf->lfd,&lf->sp,&lf->dv,des);

  lf->fp.h[v] = des->h;
  for (i=0; i<des->ncoef; i++) coef[i] = des->cf[cfn(des,i)];

  if (!lf_error)
  { if (dc(&lf->fp)) dercor(&lf->lfd,&lf->sp,des,coef);
    subparcomp(des,lf,coef);
    for (i=0; i<des->ncoef; i++) lf->fp.coef[i*lf->fp.nvm+v] =  coef[i];
  }

  lf->fp.deg[v] = deg(&lf->sp);

  return(lf_status);
}

/*
 * Set default values for the likelihood e.t.c. This
 * is called in cases where the optimization for the fit
 * has failed.
 */

void set_default_like(fp,v)
fitpt *fp;
int v;
{ int i, nvm, d;
  nvm = fp->nvm;
  d = fp->d;
  fp->lik[v] = fp->lik[nvm+v] = 0;
  fp->lik[2*nvm+v] = 0; /* should use sum of weights here? */
  for (i=0; i<=d; i++)
    fp->t0[i*nvm+v] = fp->nlx[i*nvm+v] = 0.0;
}

int procv(des,lf,v)
design *des;
lfit *lf;
int v;
{ int d, p, nvm, i, k;
  double trc[6], t0[1+MXDIM], vari[1+MXDIM];
  k = procvraw(des,lf,v);
  if (lf_error) return(k);
   
  d = lf->lfd.d;
  p = npar(&lf->sp);
  nvm = lf->fp.nvm;

  switch(k)
  { case LF_OK: break;
    case LF_NCON:
      WARN(("procv: locfit did not converge"));
      break;
    case LF_OOB:
      WARN(("procv: parameters out of bounds"));
      break;
    case LF_PF:
      if (lf_debug>1) WARN(("procv: perfect fit"));
      set_default_like(&lf->fp,v);
      return(k);
    case LF_NOPT:
      WARN(("procv: no points with non-zero weight"));
      set_default_like(&lf->fp,v);
      return(k);
    case LF_INFA:
      if (lf_debug>1) WARN(("procv: initial value problem"));
      set_default_like(&lf->fp,v);
      return(k);
    case LF_DEMP:
      WARN(("procv: density estimate, empty integration region"));
      set_default_like(&lf->fp,v);
      return(k);
    case LF_XOOR:
      WARN(("procv: fit point outside xlim region"));
      set_default_like(&lf->fp,v);
      return(k);
    case LF_DNOP:
      if (lf_debug>1)
        WARN(("density estimation -- insufficient points in smoothing window"));
      set_default_like(&lf->fp,v);
      return(k);
    case LF_FPROB:
      WARN(("procv: f problem; likelihood failure"));
      set_default_like(&lf->fp,v);
      return(k);
    default:
      WARN(("procv: unknown return code %d",k));
      set_default_like(&lf->fp,v);
      return(k);
  }

  comp_vari(&lf->lfd,&lf->sp,des,trc,t0);
  lf->fp.lik[v] = des->llk;
  lf->fp.lik[nvm+v] = trc[2];
  lf->fp.lik[2*nvm+v] = trc[0]-trc[2];

  for (i=0; i<des->ncoef; i++)
    vari[i] = des->V[p*cfn(des,0) + cfn(des,i)];
  vari[0] = sqrt(vari[0]);
  if (vari[0]>0) for (i=1; i<des->ncoef; i++) vari[i] /= vari[0];
  t0[0] = sqrt(t0[0]);
  if (t0[0]>0) for (i=1; i<des->ncoef; i++) t0[i] /= t0[0];

  subparcomp2(des,lf,vari,t0);
  for (i=0; i<des->ncoef; i++)
  { lf->fp.nlx[i*nvm+v] = vari[i];
    lf->fp.t0[i*nvm+v]  = t0[i];
  }

  return(k);
}

double intvo(des,lf,c0,c1,a,p,t0,t20,t21)
design *des;
lfit *lf;
double *c0, *c1, a, t0, t20, t21;
int p;
{ double th, lk, link[LLEN];
  int i;
  lk = 0;
  for (i=0; i<des->n; i++)
  { th = (1-a)*innerprod(c0,&des->X[i*p],p) + a*innerprod(c1,&des->X[i*p],p);
    stdlinks(link,&lf->lfd,&lf->sp,(int)des->ind[i],th,robscale);
    lk += des->w[i]*link[ZLIK];
  }
  des->llk = lk;
  return(vocri(des->llk,t0,(1-a)*t20+a*t21,pen(&lf->sp)));
}

int procvvord(des,lf,v)
design *des;
lfit *lf;
int v;
{ double tr[6], gcv, g0, ap, coef[4][10], t2[4], th, md=0.0;
  int i, j, k=0, d1, i0, p1, ip;
  des->xev = evpt(&lf->fp,v);

  ap = pen(&lf->sp);
  if ((ap==0) & ((fam(&lf->sp)&63)!=TGAUS)) ap = 2.0;
  d1 = deg(&lf->sp); p1 = npar(&lf->sp);
  for (i=0; i<p1; i++) coef[0][i] = coef[1][i] = coef[2][i] = coef[3][i] = 0.0;
  i0 = 0; g0 = 0;
  ip = 1;

  for (i=deg0(&lf->sp); i<=d1; i++)
  { deg(&lf->sp) = i;
    des->p = npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);
    k = locfit(&lf->lfd,des,&lf->sp,0, i==deg0(&lf->sp),0);

    local_df(&lf->lfd,&lf->sp,des,tr);
    gcv = vocri(des->llk,tr[0],tr[2],ap);
    if ((i==deg0(&lf->sp)) || (gcv<g0)) { i0 = i; g0 = gcv; md = i; }

    for (j=0; j<des->p; j++) coef[i][j] = des->cf[j];
    t2[i] = tr[2];

#ifdef RESEARCH
printf("variable order\n");
    if ((ip) && (i>deg0(&lf->sp)))
    { for (j=1; j<10; j++)
      { gcv = intvo(des,lf,coef[i-1],coef[i],j/10.0,des->p,tr[0],t2[i-1],t2[i]);
        if (gcv<g0) { g0 = gcv; md = i-1+j/10.0; }
      }
    }
#endif
  }
  lf->fp.h[v] = des->h;
  if (lf->fp.h[v]<=0) WARN(("zero bandwidth in procvvord"));

  if (i0<d1) /* recompute the best fit */
  { deg(&lf->sp) = i0;
    des->p = npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);
    k = locfit(&lf->lfd,des,&lf->sp,0,0,0);
    for (i=npar(&lf->sp); i<p1; i++) des->cf[i] = 0.0;
    i0 = md; if (i0==d1) i0--;
    th = md-i0;
    for (i=0; i<p1; i++) des->cf[i] = (1-th)*coef[i0][i]+th*coef[i0+1][i];
    deg(&lf->sp) = d1; npar(&lf->sp) = p1;
  }

  for (i=0; i<p1; i++) lf->fp.coef[i*lf->fp.nvm+v] = des->cf[i];
  lf->fp.deg[v] = md;
  return(k);
}

int procvhatm(des,lf,v)
design *des;
lfit *lf;
int v;
{ int k=0;
  double *l;
  l = &lf->fp.L[v*lf->lfd.n];
  if ((ker(&lf->sp)!=WPARM) | (!haspc(&lf->pc)))
  { k = procvraw(des,lf,v);
    wdiag(&lf->lfd,&lf->sp,des,l,&lf->dv,0,1,1);
  }
  else
    wdiagp(&lf->lfd,&lf->sp,des,l,&lf->pc,&lf->dv,0,1,1);
  return(k);
}
