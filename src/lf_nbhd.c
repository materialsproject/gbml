/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *  Functions for determining bandwidth; smoothing neighborhood
 *  and smoothing weights.
 */

#include "local.h"

double rho(x,sc,d,kt,sty) /* ||x|| for appropriate distance metric */
double *x, *sc;
int d, kt, *sty;
{ double rhoi[MXDIM], s;
  int i;
  for (i=0; i<d; i++)
  { if (sty!=NULL)
    { switch(sty[i])
      { case STANGL:  rhoi[i] = 2*sin(x[i]/(2*sc[i])); break;
        case STCPAR: rhoi[i] = 0; break;
        default: rhoi[i] = x[i]/sc[i];
    } }
    else rhoi[i] = x[i]/sc[i];
  }

  if (d==1) return(fabs(rhoi[0]));

  s = 0;
  if (kt==KPROD)
  { for (i=0; i<d; i++)
    { rhoi[i] = fabs(rhoi[i]);
      if (rhoi[i]>s) s = rhoi[i];
    }
    return(s);
  }

  if (kt==KSPH)
  { for (i=0; i<d; i++)
      s += rhoi[i]*rhoi[i];
    return(sqrt(s));
  }

  ERROR(("rho: invalid kt"));
  return(0.0);
}

double kordstat(x,k,n,ind)
double *x;
int k, n;
Sint *ind;
{ int i, i0, i1, l, r;
  double piv;
  if (k<1) return(0.0);
  i0 = 0; i1 = n-1;
  while (1)
  { piv = x[ind[(i0+i1)/2]];
    l = i0; r = i1;
    while (l<=r)
    { while ((l<=i1) && (x[ind[l]]<=piv)) l++;
      while ((r>=i0) && (x[ind[r]]>piv)) r--;
      if (l<=r) ISWAP(ind[l],ind[r]);
    } /* now, x[ind[i0..r]] <= piv < x[ind[l..i1]] */
    if (r<k-1) i0 = l;  /* go right */
    else /* put pivots in middle */
    { for (i=i0; i<=r; )
        if (x[ind[i]]==piv) { ISWAP(ind[i],ind[r]); r--; }
        else i++;
      if (r<k-1) return(piv);
      i1 = r;
    }
  }
}

/* check if i'th data point is in limits */
int inlim(lfd,i)
lfdata *lfd;
int i;
{ int d, j, k;
  double *xlim;

  xlim = lfd->xl;
  d = lfd->d;
  k = 1;
  for (j=0; j<d; j++)
  { if (xlim[j]<xlim[j+d])
      k &= ((datum(lfd,j,i)>=xlim[j]) & (datum(lfd,j,i)<=xlim[j+d]));
  }
  return(k);
}

double compbandwid(di,ind,x,n,d,nn,fxh)
double *di, *x, fxh;
Sint *ind;
int n, d, nn;
{ int i;
  double nnh;

  if (nn==0) return(fxh);

  if (nn<n)
    nnh = kordstat(di,nn,n,ind);
  else
  { nnh = 0;
    for (i=0; i<n; i++) nnh = MAX(nnh,di[i]);
    nnh = nnh*exp(log(1.0*nn/n)/d);
  }
  return(MAX(fxh,nnh));
}

/*
  fast version of nbhd for ordered 1-d data
*/
void nbhd1(lfd,sp,des,k)
lfdata *lfd;
smpar *sp;
design *des;
int k;
{ double x, h, *xd, sc;
  int i, l, r, m, n, z;

  n = lfd->n;
  x = des->xev[0];
  xd = dvari(lfd,0);
  sc = lfd->sca[0];

  /* find closest data point to x */
  if (x<=xd[0]) z = 0;
  else
  if (x>=xd[n-1]) z = n-1;
  else
  { l = 0; r = n-1;
    while (r-l>1)
    { z = (r+l)/2;
      if (xd[z]>x) r = z;
              else l = z;
    }
    /* now, xd[0..l] <= x < x[r..n-1] */
    if ((x-xd[l])>(xd[r]-x)) z = r; else z = l;
  }
  /* closest point to x is xd[z] */

  if (nn(sp)<0)  /* user bandwidth */
    h = sp->vb(des->xev);
  else
  { if (k>0) /* set h to nearest neighbor bandwidth */
    { l = r = z;
      if (l==0) r = k-1;
      if (r==n-1) l = n-k;
      while (r-l<k-1)
      { if ((x-xd[l-1])<(xd[r+1]-x)) l--; else r++;
        if (l==0) r = k-1;
        if (r==n-1) l = n-k;
      }
      h = x-xd[l];
      if (h<xd[r]-x) h = xd[r]-x;
    }
    else h = 0;
    h /= sc;
    if (h<fixh(sp)) h = fixh(sp);
  }

  m = 0;
  if (xd[z]>x) z--; /* so xd[z]<=x */
  /* look left */
  for (i=z; i>=0; i--) if (inlim(lfd,i))
  { des->di[i] = (x-xd[i])/sc;
    des->w[m] = weight(lfd, sp, &xd[i], &x, h, 1, des->di[i]);
    if (des->w[m]>0)
    { des->ind[m] = i;
      m++; 
    } else i = 0;
  }
  /* look right */
  for (i=z+1; i<n; i++) if (inlim(lfd,i))
  { des->di[i] = (xd[i]-x)/sc;
    des->w[m] = weight(lfd, sp, &xd[i], &x, h, 1, des->di[i]);
    if (des->w[m]>0)
    { des->ind[m] = i;
      m++; 
    } else i = n;
  }

  des->n = m;
  des->h = h;
}

void nbhd_zeon(lfd,des)
lfdata *lfd;
design *des;
{ int i, j, m, eq;

  m = 0;
  for (i=0; i<lfd->n; i++)
  { eq = 1;
    for (j=0; j<lfd->d; j++) eq = eq && (des->xev[j] == datum(lfd,j,i));
    if (eq)
    { des->w[m] = 1;
      des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  des->h = 1.0;
}

void nbhd(lfd,des,nn,redo,sp)
lfdata *lfd;
design *des;
int redo, nn;
smpar *sp;
{ int d, i, j, m, n;
  double h, u[MXDIM];

  if (lf_debug>1) printf("nbhd: nn %d  fixh %8.5f\n",nn,fixh(sp));
  
  d = lfd->d; n = lfd->n;

  if (ker(sp)==WPARM)
  { for (i=0; i<n; i++)
    { des->w[i] = 1.0;
      des->ind[i] = i;
    }
    des->n = n;
    return;
  }

  if (kt(sp)==KZEON)
  { nbhd_zeon(lfd,des);
    return;
  }

  if (kt(sp)==KCE)
  { des->h = 0.0;
    return;
  }

  /* ordered 1-dim; use fast searches */
  if ((nn<=n) & (lfd->ord) & (ker(sp)!=WMINM) & (lfd->sty[0]!=STANGL))
  { nbhd1(lfd,sp,des,nn);
    return;
  }

  if (!redo)
  { for (i=0; i<n; i++)
    { for (j=0; j<d; j++) u[j] = datum(lfd,j,i)-des->xev[j];
      des->di[i] = rho(u,lfd->sca,d,kt(sp),lfd->sty);
      des->ind[i] = i;
    }
  }
  else
    for (i=0; i<n; i++) des->ind[i] = i;

  if (ker(sp)==WMINM)
  { des->h = minmax(lfd,des,sp);
    return;
  }

  if (nn<0)
    h = sp->vb(des->xev);
  else
    h = compbandwid(des->di,des->ind,des->xev,n,lfd->d,nn,fixh(sp));
  m = 0;
  for (i=0; i<n; i++) if (inlim(lfd,i))
  { for (j=0; j<d; j++) u[j] = datum(lfd,j,i);
    des->w[m] = weight(lfd, sp, u, des->xev, h, 1, des->di[i]);
    if (des->w[m]>0)
    { des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  des->h = h;
}
