/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *   Routines for building and interpolating the kd tree.
 *   Initially, this started from the loess code.
 *
 *   Todo: EKDCE isn't working.
 */

#include "local.h"

void newcell();
static int nterm;

void kdtre_guessnv(evs,nvm,ncm,vc,n,d,alp)
evstruc *evs;
double alp;
int *nvm, *ncm, *vc, n, d;
{ int k;
  if (ev(evs) == EKDTR)
  { nterm = (int)(cut(evs)/4 * n * MIN(alp,1.0) );
    k = 2*n/nterm;
    *vc = 1<<d;
    *ncm = 2*k+1;
    *nvm = (k+2)**vc/2;
    return;
  }
  if (ev(evs) == EKDCE)
  { nterm = (int)(n * alp);
    *vc = 1;
    *nvm = 1+(int)(2*n/nterm);
    *ncm = 2**nvm+1;
    return;
  }
  *nvm = *ncm = *vc = 0;
  return;
}

/*
  Split x[pi[l..r]] into two equal sized sets.

  Let m=(l+r)/2.
  At return,
    x[pi[l..m]] < x[pi[m+1..r]].
    Assuming no ties:
      If l+r is odd, the sets have the same size.
      If l+r is even, the low set is larger by 1.
    If there are ties, all ties go in the low set.
*/      
int ksmall(l, r, m, x, pi)
Sint *pi;
int l, r, m;
double *x;
{
  int il, ir, jl, jr;
  double t;


  while(l<r)
  { t = x[pi[m]];

    /*
      permute the observations so that
        x[pi[l..il]] < t <= x[pi[ir..r]].
    */
    ir = l; il = r;
    while (ir<il)
    { while ((ir<=r) && (x[pi[ir]] < t)) ir++;
      while ((il>=l) && (x[pi[il]]>= t)) il--;
      if (ir<il) ISWAP(pi[ir],pi[il]);
    }

    /*
      move  = t to the middle:
        x[pi[l..il]]  < t
        x[pi[jl..jr]] = t
        x[pi[ir..r]] > t
    */
    jl = ir; jr = r;
    while (ir<jr)
    { while ((ir<=r)  && (x[pi[ir]]== t)) ir++;
      while ((jr>=jl) && (x[pi[jr]] > t)) jr--;
      if (ir<jr) ISWAP(pi[ir],pi[jr]);
    }

    /*
      we're done if m is in the middle, jl <= m <= jr.
    */
    if ((jl<=m) & (jr>=m)) return(jr);

    /*
      update l or r.
    */
    if (m>=ir) l = ir;
    if (m<=il) r = il;
  }
  if (l==r) return(l);
  ERROR(("ksmall failure"));
  return(0);
}

int terminal(lf,p,pi,fc,d,m,split_val)
lfit *lf;
Sint *pi;
int p, d, fc, *m;
double *split_val;
{ int i, k, lo, hi, split_var;
  double max, min, score, max_score, t;

  /*
    if there are fewer than fc points in the cell, this cell
    is terminal.
  */
  lo = lf->evs.lo[p]; hi = lf->evs.hi[p];
  if (hi-lo < fc) return(-1);

  /* determine the split variable */
  max_score = 0.0; split_var = 0;
  for (k=0; k<d; k++)
  { max = min = datum(&lf->lfd, k, pi[lo]);
    for (i=lo+1; i<=hi; i++)
    { t = datum(&lf->lfd,k,pi[i]);
      if (t<min) min = t;
      if (t>max) max = t;
    }
    score = (max-min) / lf->lfd.sca[k];
    if (score > max_score)
    { max_score = score;
      split_var = k;
    }
  }
  if (max_score==0) /* all points in the cell are equal */
    return(-1);

  *m = ksmall(lo,hi,(lo+hi)/2, dvari(&lf->lfd,split_var), pi);
  *split_val = datum(&lf->lfd, split_var, pi[*m]);

  if (*m==hi) /* all observations go lo */
    return(-1);
  return(split_var);
}

void kdtre_start(des,lf)
design *des;
lfit *lf;
{ Sint *pi;
  int i, j, vc, d, nc, nv, ncm, nvm, k, m, n, p;
  double sv;
  d = lf->lfd.d; n = lf->lfd.n; pi = des->ind;
  kdtre_guessnv(&lf->evs,&nvm,&ncm,&vc,n,d,nn(&lf->sp));
  trchck(lf,nvm,ncm,vc);

  nv = 0;
  if (ev(&lf->evs) != EKDCE)
  { for (i=0; i<vc; i++)
    { j = i;
      for (k=0; k<d; ++k)
      { evptx(&lf->fp,i,k) = lf->evs.fl[d*(j%2)+k];
        j >>= 1;
      }
    }
    nv = vc;
    for (j=0; j<vc; j++) lf->evs.ce[j] = j;
  }

  for (i=0; i<n; i++) pi[i] = i;
  p = 0; nc = 1;
  lf->evs.lo[p] = 0; lf->evs.hi[p] = n-1;
  lf->evs.s[p] = -1;
  while (p<nc)
  { k = terminal(lf,p,pi,nterm,d,&m,&sv);
    if (k>=0)
    {
      if ((ncm<nc+2) | (2*nvm<2*nv+vc))
      { WARN(("Insufficient space for full tree"));
        lf->evs.nce = nc; lf->fp.nv = nv;
        return;
      }

      /* new lo cell has obsn's lo[p]..m */
      lf->evs.lo[nc] = lf->evs.lo[p];
      lf->evs.hi[nc] = m;
      lf->evs.s[nc] = -1;

      /* new hi cell has obsn's m+1..hi[p] */
      lf->evs.lo[nc+1] = m+1;
      lf->evs.hi[nc+1] = lf->evs.hi[p];
      lf->evs.s[nc+1] = -1;

      /* cell p is split on variable k, value sv */
      lf->evs.s[p] = k;
      lf->evs.sv[p] = sv;
      lf->evs.lo[p] = nc; lf->evs.hi[p] = nc+1;

      nc=nc+2; i = nv;

      /* now compute the new vertices. */
      if (ev(&lf->evs) != EKDCE)
        newcell(&nv,vc,evp(&lf->fp), d, k, sv,
             &lf->evs.ce[p*vc], &lf->evs.ce[(nc-2)*vc], &lf->evs.ce[(nc-1)*vc]);

    }
    else if (ev(&lf->evs)==EKDCE) /* new vertex at cell center */
    { sv = 0;
      for (i=0; i<d; i++) evptx(&lf->fp,nv,i) = 0;
      for (j=lf->evs.lo[p]; j<=lf->evs.hi[p]; j++)
      { sv += prwt(&lf->lfd,(int)pi[j]);
        for (i=0; i<d; i++)
          evptx(&lf->fp,nv,i) += datum(&lf->lfd,i,pi[j])*prwt(&lf->lfd,(int)pi[j]);
      }
      for (i=0; i<d; i++) evptx(&lf->fp,nv,i) /= sv;
      lf->lfd.n = lf->evs.hi[p] - lf->evs.lo[p] + 1;
      des->ind = &pi[lf->evs.lo[p]]; /* why? */
      des->vfun(des,lf,nv);
      lf->lfd.n = n; des->ind = pi;
      nv++;
    }
    p++;
  }

  /* We've built the tree. Now do the fitting. */
  if (ev(&lf->evs)==EKDTR)
    for (i=0; i<nv; i++) des->vfun(des,lf,i);

  lf->evs.nce = nc; lf->fp.nv = nv;
  return;
}

void newcell(nv,vc,xev, d, k, split_val, cpar, clef, crig)
double *xev, split_val;
Sint *cpar, *clef, *crig;
int *nv, vc, d, k;
{ int i, ii, j, j2, tk, match;
  tk = 1<<k;
  for (i=0; i<vc; i++)
  { if ((i&tk) == 0)
    { for (j=0; j<d; j++) xev[*nv*d+j] = xev[d*cpar[i]+j];
      xev[*nv*d+k] = split_val;
      match = 0; j = vc; /* no matches in first vc points */
      while ((j<*nv) && (!match))
      { j2 = 0;
        while ((j2<d) && (xev[*nv*d+j2] == xev[j*d+j2])) j2++;
        match = (j2==d);
        if (!match) j++;
      }
      ii = i+tk;
      clef[i] = cpar[i];
      clef[ii]= crig[i] = j;
      crig[ii]= cpar[ii];
      if (!match) (*nv)++;
    }
  }
  return;
}

extern void hermite2();

double blend(fp,evs,s,x,ll,ur,j,nt,t,what)
fitpt *fp;
evstruc *evs;
double s, *x, *ll, *ur;
int j, nt, *t, what;
{ Sint *ce;
  int k, k1, m, nc, j0, j1;
  double v0, v1, xibar, g0[3], g1[3], gg[4], gp[4], phi[4];
  ce = evs->ce;
  for (k=0; k<4; k++)  /* North South East West */
  { k1 = (k>1);
    v0 = ll[k1]; v1 = ur[k1];
    j0 = ce[j+2*(k==0)+(k==2)];
    j1 = ce[j+3-2*(k==1)-(k==3)];
    xibar = (k%2==0) ? ur[k<2] : ll[k<2];
    m = nt;
    while ((m>=0) && ((evs->s[t[m]] != (k<=1)) | (evs->sv[t[m]] != xibar))) m--;
    if (m >= 0)
    { m = (k%2==1) ? evs->lo[t[m]] : evs->hi[t[m]];
      while (evs->s[m] != -1)
        m = (x[evs->s[m]] < evs->sv[m]) ? evs->lo[m] : evs->hi[m];
      if (v0 < evptx(fp,ce[4*m+2*(k==1)+(k==3)],k1))
      { j0 = ce[4*m+2*(k==1)+(k==3)];
        v0 = evptx(fp,j0,k1);
      }
      if (evptx(fp,ce[4*m+3-2*(k==0)-(k==2)],k1) < v1)
      { j1 = ce[4*m+3-2*(k==0)-(k==2)];
        v1 = evptx(fp,j1,k1);
      }
    }
    nc = exvval(fp,g0,j0,2,what,0);
    nc = exvval(fp,g1,j1,2,what,0);
    if (nc==1)
      gg[k] = linear_interp((x[(k>1)]-v0),v1-v0,g0[0],g1[0]);
    else
    { hermite2(x[(k>1)]-v0,v1-v0,phi);
      gg[k] = phi[0]*g0[0]+phi[1]*g1[0]+(phi[2]*g0[1+k1]+phi[3]*g1[1+k1])*(v1-v0);
      gp[k] = phi[0]*g0[2-k1] + phi[1]*g1[2-k1];
    }
  }
  s = -s;
  if (nc==1)
    for (k=0; k<2; k++)
      s += linear_interp(x[k]-ll[k],ur[k]-ll[k],gg[3-2*k],gg[2-2*k]);
    else
    for (k=0; k<2; k++) /* EW NS */
    { hermite2(x[k]-ll[k],ur[k]-ll[k],phi);
      s += phi[0]*gg[3-2*k] + phi[1]*gg[2-2*k]
          +(phi[2]*gp[3-2*k] + phi[3]*gp[2-2*k]) * (ur[k]-ll[k]);
    }
  return(s);
}

double kdtre_int(fp,evs,x,what)
fitpt *fp;
evstruc *evs;
double *x;
int what;
{ Sint *ce;
  int k, vc, t[20], nt, nc, j, d;
  double *ll, *ur, ff, vv[64][64];
  d = fp->d;
  vc = 1<<d;
  if (d > 6) ERROR(("d too large in kdint"));

  /* descend the tree to find the terminal cell */
  nt = 0; t[nt] = 0; k = 0;
  while (evs->s[k] != -1)
  { nt++;
    if (nt>=20) { ERROR(("Too many levels in kdint")); return(NOSLN); }
    k = t[nt] = (x[evs->s[k]] < evs->sv[k]) ? evs->lo[k] : evs->hi[k];
  }

  ce = &evs->ce[k*vc];
  ll = evpt(fp,ce[0]);
  ur = evpt(fp,ce[vc-1]);
  nc = 0;
  for (j=0; j<vc; j++) nc = exvval(fp,vv[j],(int)ce[j],d,what,0);
  ff = rectcell_interp(x,vv,ll,ur,d,nc);

  if (d==2) ff = blend(fp,evs,ff,x,ll,ur,k*vc,nt,t,what);
  return(ff);
}
