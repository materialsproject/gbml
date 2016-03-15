/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

double linear_interp(h,d,f0,f1)
double h, d, f0, f1;
{ if (d==0) return(f0);
  return( ( (d-h)*f0 + h*f1 ) / d );
}

void hermite2(x,z,phi)
double x, z, *phi;
{ double h;
  if (z==0)
  { phi[0] = 1.0; phi[1] = phi[2] = phi[3] = 0.0;
    return;
  }
  h = x/z;
  if (h<0)
  { phi[0] = 1; phi[1] = 0;
    phi[2] = h; phi[3] = 0;
    return;
  }
  if (h>1)
  { phi[0] = 0; phi[1] = 1;
    phi[2] = 0; phi[3] = h-1;
    return;
  }
  phi[1] = h*h*(3-2*h);
  phi[0] = 1-phi[1];
  phi[2] = h*(1-h)*(1-h);
  phi[3] = h*h*(h - 1);
}

double cubic_interp(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  hermite2(h,1.0,phi);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

double cubintd(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  phi[1] = 6*h*(1-h);
  phi[0] = -phi[1];
  phi[2] = (1-h)*(1-3*h);
  phi[3] = h*(3*h-2);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

/*
  interpolate over a rectangular cell.
    x = interpolation point. 
    vv = array of vertex values.
    ll = lower left corner.
    ur = upper right corner.
    d = dimension.
    nc = no of coefficients.
*/
double rectcell_interp(x,vv,ll,ur,d,nc)
double *x, vv[64][64], *ll, *ur;
int d, nc;
{ double phi[4];
  int i, j, k, tk;

  tk = 1<<d;
  for (i=0; i<tk; i++) if (vv[i][0]==NOSLN) return(NOSLN);

  /* no derivatives - use multilinear interpolation */
  if (nc==1)
  { for (i=d-1; i>=0; i--)
    { tk = 1<<i;
      for (j=0; j<tk; j++)
        vv[j][0] = linear_interp(x[i]-ll[i],ur[i]-ll[i],vv[j][0],vv[j+tk][0]);
    }
    return(vv[0][0]);
  }

  /* with derivatives -- use cubic */
  if (nc==d+1)
  { for (i=d-1; i>=0; i--)
    { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
      tk = 1<<i;
      phi[2] *= ur[i]-ll[i];
      phi[3] *= ur[i]-ll[i];
      for (j=0; j<tk; j++)
      { vv[j][0] = phi[0]*vv[j][0] + phi[1]*vv[j+tk][0]
                 + phi[2]*vv[j][i+1] + phi[3]*vv[j+tk][i+1];
        for (k=1; k<=i; k++)
          vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k];
      }
    }
    return(vv[0][0]); 
  }

  /* with all coefs -- use multicubic */
  for (i=d-1; i>=0; i--)
  { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
    tk = 1<<i;
    phi[2] *= ur[i]-ll[i];
    phi[3] *= ur[i]-ll[i];
    for (j=0; j<tk; j++)
      for (k=0; k<tk; k++)
        vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k]
                 + phi[2]*vv[j][k+tk] + phi[3]*vv[j+tk][k+tk];
  }
  return(vv[0][0]);
}

int exvval(fp,vv,nv,d,what,z)
fitpt *fp;
double *vv;
int nv, d, z, what;
{ int i, k;
  double *values;

  k = (z) ? 1<<d : d+1;
  for (i=1; i<k; i++) vv[i] = 0.0;
  switch(what)
  { case PCOEF:
      values = fp->coef;
      break;
    case PVARI:
    case PNLX:
      values = fp->nlx;
      break;
    case PT0:
      values = fp->t0;
      break;
    case PBAND:
      vv[0] = fp->h[nv];
      return(1);
    case PDEGR:
      vv[0] = fp->deg[nv];
      return(1);
    case PLIK:
      vv[0] = fp->lik[nv];
      return(1);
    case PRDF:
      vv[0] = fp->lik[2*fp->nvm+nv];
      return(1);
    default:
      ERROR(("Invalid what in exvval"));
      return(0);
  }
  vv[0] = values[nv];
  if (!fp->hasd) return(1);
  if (z)
  { for (i=0; i<d; i++) vv[1<<i] = values[(i+1)*fp->nvm+nv];
    return(1<<d);
  }
  else
  { for (i=1; i<=d; i++) vv[i] = values[i*fp->nvm+nv];
    return(d+1);
  }
}

void exvvalpv(vv,vl,vr,d,k,dl,nc)
double *vv, *vl, *vr, dl;
int d, k, nc;
{ int i, tk, td;
  double f0, f1;
  if (nc==1)
  { vv[0] = (vl[0]+vr[0])/2;
    return;
  }
  tk = 1<<k;
  td = 1<<d;
  for (i=0; i<td; i++) if ((i&tk)==0)
  { f0 = (vl[i]+vr[i])/2 + dl*(vl[i+tk]-vr[i+tk])/8;
    f1 = 1.5*(vr[i]-vl[i])/dl - (vl[i+tk]+vr[i+tk])/4;
    vv[i] = f0;
    vv[i+tk] = f1;
  }
} 

double grid_int(fp,evs,x,what)
fitpt *fp;
evstruc *evs;
double *x;
int what;
{ int d, i, j, jj, nc=0, sk, v[MXDIM], vc, z0, nce[1<<MXDIM], *mg;
  double *ll, *ur, vv[64][64], z;

  d = fp->d;
  ll = evpt(fp,0); ur = evpt(fp,fp->nv-1);
  mg = mg(evs);

  z0 = 0; vc = 1<<d;
  for (j=d-1; j>=0; j--)
  { v[j] = (int)((mg[j]-1)*(x[j]-ll[j])/(ur[j]-ll[j]));
    if (v[j]<0) v[j]=0;
    if (v[j]>=mg[j]-1) v[j] = mg[j]-2;
    z0 = z0*mg[j]+v[j];
  }
  nce[0] = z0; nce[1] = z0+1; sk = jj = 1; 
  for (i=1; i<d; i++)
  { sk *= mg[i-1];
    jj<<=1;
    for (j=0; j<jj; j++)
      nce[j+jj] = nce[j]+sk;
  }
  for (i=0; i<vc; i++)
    nc = exvval(fp,vv[i],nce[i],d,what,1);
  ll = evpt(fp,nce[0]);
  ur = evpt(fp,nce[vc-1]);
  z = rectcell_interp(x,vv,ll,ur,d,nc);
  return(z);
}

double fitp_int(fp,x,what,i)
fitpt *fp;
double *x;
int what, i;
{ double vv[1+MXDIM];
  exvval(fp,vv,i,fp->d,what,0);
  return(vv[0]);
}

double xbar_int(fp,x,what)
fitpt *fp;
double *x;
int what;
{ int i, nc;
  double vv[1+MXDIM], f;
  nc = exvval(fp,vv,0,fp->d,what,0);
  f = vv[0];
  if (nc>1)
    for (i=0; i<fp->d; i++)
      f += vv[i+1]*(x[i]-evptx(fp,0,i));
  return(f);
}

double dointpoint(lf,x,what,ev,j)
lfit *lf;
double *x;
int what, ev, j;
{ double xf, f=0.0;
  int i;
  fitpt *fp;
  evstruc *evs;

  fp = &lf->fp; evs = &lf->evs;
  for (i=0; i<fp->d; i++) if (lf->lfd.sty[i]==STANGL)
  { xf = floor(x[i]/(2*PI*lf->lfd.sca[i]));
    x[i] -= xf*2*PI*lf->lfd.sca[i];
  }

  switch(ev)
  { case EGRID: f = grid_int(fp,evs,x,what); break;
    case EKDTR: f = kdtre_int(fp,evs,x,what); break;
    case ETREE: f = atree_int(lf,x,what); break;
    case EPHULL: f = triang_int(lf,x,what); break;
    case EFITP: f = fitp_int(fp,x,what,j); break;
    case EXBAR: f = xbar_int(fp,x,what); break;
    case ENONE: f = 0; break;
    case ESPHR: f = sphere_int(lf,x,what); break;
    default: ERROR(("dointpoint: cannot interpolate structure %d",ev));
  }
  if (((what==PT0)|(what==PNLX)) && (f<0)) f = 0.0;
  f += addparcomp(lf,x,what);
  return(f);
}
