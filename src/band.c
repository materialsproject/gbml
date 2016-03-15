/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern void fitoptions();

static double hmin, gmin, sig2, pen, vr, tb;
static lfit *blf;
static design *bdes;

int procvbind(des,lf,v)
design *des;
lfit *lf;
int v;
{ double s0, s1, bi;
  int i, ii, k;
  k = procvraw(des,lf,v);
  wdiag(&lf->lfd, &lf->sp, des,des->wd,&lf->dv,0,1,0);
  s0 = s1 = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    s0+= prwt(&lf->lfd,ii)*des->wd[i]*des->wd[i];
    bi = prwt(&lf->lfd,ii)*fabs(des->wd[i]*ipower(des->di[ii],deg(&lf->sp)+1));
    s1+= bi*bi;
  }
  vr += s0;
  tb += s1;
  return(k);
}

double bcri(h,c,cri)
double h;
int c, cri;
{ double num, den;
  int (*pv)();
  if (c==DALP)
    blf->sp.nn = h;
  else
    blf->sp.fixh = h;
  if ((cri&63)==BIND)
  { pv = procvbind;
    vr = tb = 0.0;
  }
  else pv = procv;
  if (cri<64) startlf(bdes,blf,pv,0);
  switch(cri&63)
  { case BGCV:
      ressumm(blf,bdes);
      num = -2*blf->lfd.n*llk(&blf->fp);
      den = blf->lfd.n-df0(&blf->fp);
      return(num/(den*den));
    case BCP:
      ressumm(blf,bdes);
      return(-2*llk(&blf->fp)/sig2-blf->lfd.n+pen*df0(&blf->fp));
    case BIND:
      return(vr+pen*pen*tb);
  } 
  ERROR(("bcri: unknown criterion"));
  return(0.0);
}

void bsel2(h0,g0,ifact,c,cri)
double h0, g0, ifact;
int c, cri;
{ int done, inc;
  double h1, g1;
  h1 = h0; g1 = g0;
  done = inc = 0;
  while (!done)
  { h1 *= 1+ifact;
    g0 = g1;
    g1 = bcri(h1,c,cri);
    if (g1<gmin) { hmin = h1; gmin = g1; }
    if (g1>g0) inc++; else inc = 0;
    switch(cri)
    { case BIND:
        done = (inc>=4) & (vr<blf->fp.nv);
        break;
      default:
        done = (inc>=4);
    }
  }
}

void bsel3(h0,g0,ifact,c,cri)
double h0, g0, ifact;
int c, cri;
{ double h1, g1;
  int i;
  hmin = h0; gmin = g0;
  for (i=-1; i<=1; i++) if (i!=0)
  { h1 = h0*(1+i*ifact);
    g1 = bcri(h1,c,cri);
    if (g1<gmin) { hmin = h1; gmin = g1; }
  }
  return;
}

void bselect(lf,des,c,cri,pn)
lfit *lf;
design *des;
int c, cri;
double pn;
{ double h0, g0, ifact;
  int i;
  pen = pn;
  blf = lf;
  bdes = des;
  if (cri==BIND) pen /= factorial(deg(&lf->sp)+1);
  hmin = h0 = (c==DFXH) ? fixh(&lf->sp) : nn(&lf->sp);
  if (h0==0) ERROR(("bselect: initial bandwidth is 0"));
  if (lf_error) return;
  sig2 = 1.0;

  gmin = g0 = bcri(h0,c,cri);
  if (cri==BCP)
  { sig2 = rv(&lf->fp);
    g0 = gmin = bcri(h0,c,cri+64);
  }
  
  ifact = 0.3;
  bsel2(h0,g0,ifact,c,cri);

  for (i=0; i<5; i++)
  { ifact = ifact/2;
    bsel3(hmin,gmin,ifact,c,cri);
  }
  if (c==DFXH)
    fixh(&lf->sp) = hmin;
  else
    nn(&lf->sp) = hmin;
  startlf(des,lf,procv,0);
  ressumm(lf,des);
}

double compsda(x,h,n)
double *x, h;
int n;
/* n/(n-1) * int( fhat''(x)^2 dx ); bandwidth h */
{ int i, j;
  double ik, sd, z;
  ik = wint(1,NULL,0,WGAUS);
  sd = 0;

  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
    { z = (x[i]-x[j])/h;
      sd += (2-(i==j))*Wconv4(z,WGAUS)/(ik*ik);
    }
  sd = sd/(n*(n-1)*h*h*h*h*h);
  return(sd);
}

double widthsj(x,lambda,n)
double *x, lambda;
int n;
{ double ik, a, b, td, sa, z, c, c1, c2, c3;
  int i, j;
  a = GFACT*0.920*lambda*exp(-log((double)n)/7)/SQRT2;
  b = GFACT*0.912*lambda*exp(-log((double)n)/9)/SQRT2;
  ik = wint(1,NULL,0,WGAUS);

  td = 0;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
    { z = (x[i]-x[j])/b;
      td += (2-(i==j))*Wconv6(z,WGAUS)/(ik*ik);
    }

  td = -td/(n*(n-1));
  j = 2.0;
  c1 = Wconv4(0.0,WGAUS);
  c2 = wint(1,&j,1,WGAUS);
  c3 = Wconv(0.0,WGAUS);  /* (2*c1/(c2*c3))^(1/7)=1.357 */
  sa = compsda(x,a,n);
  c = b*exp(log(c1*ik/(c2*c3*GFACT*GFACT*GFACT*GFACT)*sa/td)/7)*SQRT2;
  return(c);
}

void kdecri(x,h,res,c,k,ker,n)
double *x, h, *res, c;
int k, ker, n;
{ int i, j;
  double degfree, dfd, pen, s, r0, r1, d0, d1, ik, wij;

  if (h<=0) WARN(("kdecri, h = %6.4f",h));

  res[0] = res[1] = 0.0;
  ik = wint(1,NULL,0,ker);
  switch(k)
  { case 1: /* aic */
      pen = 2.0;
      for (i=0; i<n; i++)
      { r0 = d0 = 0.0;
        for (j=0; j<n; j++)
        { s = (x[i]-x[j])/h;
          r0 += W(s,ker);
          d0 += s*s*Wd(s,ker);
        }
        d0 = -(d0+r0)/(n*h*h*ik);  /* d0 = d/dh fhat(xi) */
        r0 /= n*h*ik;              /* r0 = fhat(xi) */
        res[0] += -2*log(r0)+pen*W(0.0,ker)/(n*h*ik*r0);
        res[1] += -2*d0/r0-pen*W(0.0,ker)/(n*h*ik*r0)*(d0/r0+1.0/h);
      }
      return;
    case 2: /* ocv */
      for (i=0; i<n; i++)
      { r0 = 0.0; d0 = 0.0;
        for (j=0; j<n; j++) if (i!=j)
        { s = (x[i]-x[j])/h;
          r0 += W(s,ker);
          d0 += s*s*Wd(s,ker);
        }
        d0 = -(d0+r0)/((n-1)*h*h);
        r0 = r0/((n-1)*h);
        res[0] -= log(r0);
        res[1] -= d0/r0;
      }
      return;
    case 3: /* lscv */
      r0 = r1 = d0 = d1 = degfree = 0.0;
      for (i=0; i<n; i++)
      { dfd = 0.0;
        for (j=0; j<n; j++)
        { s = (x[i]-x[j])/h;
          wij = W(s,ker);
          dfd += wij;
/* 
 *  r0 = \int fhat * fhat = sum_{i,j} W*W( (Xi-Xj)/h ) / n^2 h.
 *  d0 is it's derivative wrt h.
 *
 *  r1 = 1/n sum( f_{-i}(X_i) ).
 *  d1 is  it's derivative wrt h.
 *
 *  degfree = sum_i ( W_0 / sum_j W( (Xi-Xj)/h ) ) is fitted d.f.
 */
          r0 += Wconv(s,ker);
          d0 += -s*s*Wconv1(s,ker);
          if (i != j)
          { r1 += wij;
            d1 += -s*s*Wd(s,ker);
          }
        }
        degfree += 1.0/dfd;
      }
      d1 -= r1;
      d0 -= r0;
      res[0] = r0/(n*n*h*ik*ik)   - 2*r1/(n*(n-1)*h*ik);
      res[1] = d0/(n*n*h*h*ik*ik) - 2*d1/(n*(n-1)*h*h*ik);
      res[2] = degfree;
      return;
    case 4: /* bcv */
      r0 = d0 = 0.0;
      for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
        { s = (x[i]-x[j])/h;
          r0 += 2*Wconv4(s,ker);
          d0 += 2*s*Wconv5(s,ker);
        }
      d0 = (-d0-r0)/(n*n*h*h*ik*ik);
      r0 = r0/(n*n*h*ik*ik);
      j = 2.0;
      d1 = wint(1,&j,1,ker);
      r1 = Wconv(0.0,ker);
      res[0] = (d1*d1*r0/4+r1/(n*h))/(ik*ik);
      res[1] = (d1*d1*d0/4-r1/(n*h*h))/(ik*ik);
      return;
    case 5: /* sjpi */
      s = c*exp(5*log(h)/7)/SQRT2;
      d0 = compsda(x,s,n);
      res[0] = d0; /* this is S(alpha) in SJ */
      res[1] = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
      return;
    case 6: /* gas-k-k */
      s = exp(log(1.0*n)/10)*h;
      d0 = compsda(x,s,n);
      res[0] = d0;
      res[1] = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
      return;
  }
  ERROR(("kdecri: what???"));
  return;
}

double esolve(x,j,h0,h1,k,c,ker,n)
double *x, h0, h1, c;
int j, k, ker, n;
{ double h[7], d[7], r[7], res[4], min, minh, fact;
  int i, nc;
  for ( i = 0; i < 7; ++i) {
      h[i] = 0.0;
      d[i] = 0.0;
      r[i] = 0.0;
  }
  for ( i = 0; i < 4; ++i) res[i] = 0.0;
  min = 1.0e30; minh = 0.0;
  fact = 1.00001;
  h[6] = h0; kdecri(x,h[6],res,c,j,ker,n);
  r[6] = res[0]; d[6] = res[1];
  if (lf_error) return(0.0);
  nc = 0;
  for (i=0; i<k; i++)
  { h[5] = h[6]; r[5] = r[6]; d[5] = d[6];
    h[6] = h0*exp((i+1)*log(h1/h0)/k);
    kdecri(x,h[6],res,c,j,ker,n);
    r[6] = res[0]; d[6] = res[1];
    if (lf_error) return(0.0);
    if (d[5]*d[6]<0)
    { h[2] = h[0] = h[5]; d[2] = d[0] = d[5]; r[2] = r[0] = r[5];
      h[3] = h[1] = h[6]; d[3] = d[1] = d[6]; r[3] = r[1] = r[6];
      while ((h[3]>fact*h[2])|(h[2]>fact*h[3]))
      { h[4] = h[3]-d[3]*(h[3]-h[2])/(d[3]-d[2]);
        if ((h[4]<h[0]) | (h[4]>h[1])) h[4] = (h[0]+h[1])/2;
        kdecri(x,h[4],res,c,j,ker,n);
        r[4] = res[0]; d[4] = res[1];
        if (lf_error) return(0.0);
        h[2] = h[3]; h[3] = h[4];
        d[2] = d[3]; d[3] = d[4];
        r[2] = r[3]; r[3] = r[4];
        if (d[4]*d[0]>0) { h[0] = h[4]; d[0] = d[4]; r[0] = r[4]; }
                    else { h[1] = h[4]; d[1] = d[4]; r[1] = r[4]; }
      }
      if (j>=4) return(h[4]); /* first min for BCV etc */
      if (r[4]<=min) { min = r[4]; minh = h[4]; }
      nc++;
    }
  }
  if (nc==0) minh = (r[5]<r[6]) ? h0 : h1;
  return(minh);
}

void kdeselect(band,x,ind,h0,h1,meth,nm,ker,n)
double h0, h1, *band, *x;
Sint *ind;
int nm, ker, n, *meth;
{ double scale, c;
  int i, k;
  k = n/4;
  for (i=0; i<n; i++) ind[i] = i;
  scale = kordstat(x,n+1-k,n,ind) - kordstat(x,k,n,ind);
  c = widthsj(x,scale,n);
  for (i=0; i<nm; i++)
    band[i] = esolve(x,meth[i],h0,h1,10,c,ker,n);
}
