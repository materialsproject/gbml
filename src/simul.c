/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

static double pen, sig2;

void goldensec(f,des,tr,eps,xm,ym,meth)
double (*f)(), eps, *xm, *ym;
int meth;
design *des;
lfit *tr;
{ double x[4], y[4], xx[11], yy[11];
  int i, im=0;
  xx[0] = tr->sp.fixh;
  if (xx[0]<=0)
  { ERROR(("regband: initialize h>0"));
    return;
  }
  for (i=0; i<=10; i++)
  { if (i>0) xx[i] = (1+GOLDEN)*xx[i-1];
    yy[i] = f(xx[i],des,tr,meth);
    if ((i==0) || (yy[i]<yy[im])) im = i;
  }
  if (im==0) im = 1;
  if (im==10)im = 9;
  x[0] = xx[im-1]; y[0] = yy[im-1];
  x[1] = xx[im];   y[1] = yy[im];
  x[3] = xx[im+1]; y[3] = yy[im+1];
  x[2] = GOLDEN*x[3]+(1-GOLDEN)*x[0];
  y[2] = f(x[2],des,tr,meth);
  while (x[3]-x[0]>eps)
  { if (y[1]<y[2])
    { x[3] = x[2]; y[3] = y[2];
      x[2] = x[1]; y[2] = y[1];
      x[1] = GOLDEN*x[0]+(1-GOLDEN)*x[3];
      y[1] = f(x[1],des,tr,meth);
    }
    else
    { x[0] = x[1]; y[0] = y[1];
      x[1] = x[2]; y[1] = y[2];
      x[2] = GOLDEN*x[3]+(1-GOLDEN)*x[0];
      y[2] = f(x[2],des,tr,meth);
    }
  }
  im = 0;
  for (i=1; i<4; i++) if (y[i]<y[im]) im = i;
  *xm = x[im]; *ym = y[im];
}

double dnk(x,k)
double x;
int k;
{ double f;
  switch(k)
  { case 0: f = 1; break;
    case 1: f = -x; break;
    case 2: f = x*x-1; break;
    case 3: f = x*(x*x-3); break;
    case 4: f = 3-x*x*(6-x*x); break;
    case 5: f = -x*(15-x*x*(10-x*x)); break;
    case 6: f = -15+x*x*(45-x*x*(15-x*x)); break;
    default: ERROR(("dnk: k=%d too large",k)); return(0.0);
  }
  return(f*exp(-x*x/2)/S2PI);
}

double locai(h,des,lf)
double h;
design *des;
lfit *lf;
{ double cp;
  nn(&lf->sp) = h;
  startlf(des,lf,procv,0);
  ressumm(lf,des);
  cp = -2*llk(&lf->fp) + pen*df0(&lf->fp);
  return(cp);
}

double loccp(h,des,lf,m) /* m=1: cp    m=2: gcv */
double h;
design *des;
lfit *lf;
int m;
{ double cp;
  int dg, n;

  n = lf->lfd.n;
  nn(&lf->sp) = 0;
  fixh(&lf->sp) = h;
  dg = deg(&lf->sp);
  deg(&lf->sp) = deg0(&lf->sp);
  startlf(des,lf,procv,0);
  ressumm(lf,des);
  if (m==1)
    cp = -2*llk(&lf->fp)/sig2 - n + 2*df0(&lf->fp);
  else
    cp = -2*n*llk(&lf->fp)/((n-df0(&lf->fp))*(n-df0(&lf->fp)));
  printf("h %8.5f  deg %2d  rss %8.5f  trl %8.5f  cp: %8.5f\n",h,deg(&lf->sp),-2*llk(&lf->fp),df0(&lf->fp),cp);
  deg0(&lf->sp) = deg(&lf->sp);
  deg(&lf->sp) = dg;
  return(cp);
}

double cp(des,lf,meth)
design *des;
lfit *lf;
int meth;
{ double hm, ym;
  goldensec(loccp,des,lf,0.001,&hm,&ym,meth);
  return(hm);
}

double gkk(des,lf)
design *des;
lfit *lf;
{ double h, h5, nf, th;
  int i, j, n, dg0, dg1;
  ev(&lf->evs)= EDATA;
  nn(&lf->sp) = 0;
  n = lf->lfd.n;
  dg0 = deg0(&lf->sp);     /* target degree */
  dg1 = dg0+1+(dg0%2==0);  /* pilot degree */
  nf = exp(log(1.0*n)/10); /* bandwidth inflation factor */
  h = lf->sp.fixh;        /* start bandwidth */
  for (i=0; i<=10; i++)
  { deg(&lf->sp) = dg1;
    lf->sp.fixh = h*nf;
    startlf(des,lf,procv,0);
    th = 0;
    for (j=10; j<n-10; j++)
      th += lf->fp.coef[dg1*n+j]*lf->fp.coef[dg1*n+j];
th *= n/(n-20.0);
    h5 = sig2 * Wikk(ker(&lf->sp),dg0) / th;
    h = exp(log(h5)/(2*dg1+1));
/* printf("pilot %8.5f  sel %8.5f\n",lf->sp.fixh,h); */
  }
  return(h);
}

double rsw(des,lf)
design *des;
lfit *lf;
{ int i, j, k, nmax, nvm, n, mk, evo, dg0, dg1;
  double rss[6], cp[6], th22, dx, d2, hh;
  nmax = 5;
  evo = ev(&lf->evs); ev(&lf->evs) = EGRID;
  mk = ker(&lf->sp);  ker(&lf->sp) = WRECT;
  dg0 = deg0(&lf->sp);
  dg1 = 1 + dg0 + (dg0%2==0);
  deg(&lf->sp) = 4;
  for (k=nmax; k>0; k--)
  { lf->evs.mg[0] = k;
    lf->evs.fl[0] = 1.0/(2*k);
    lf->evs.fl[1] = 1-1.0/(2*k);
    nn(&lf->sp) = 0;
    fixh(&lf->sp) = 1.0/(2*k);
    startlf(des,lf,procv,0);
    nvm = lf->fp.nvm;
    rss[k] = 0;
    for (i=0; i<k; i++) rss[k] += -2*lf->fp.lik[i];
  }
  n = lf->lfd.n; k = 1;
  for (i=1; i<=nmax; i++)
  { /* cp[i] = (n-5*nmax)*rss[i]/rss[nmax]-(n-10*i); */
    cp[i] = rss[i]/sig2-(n-10*i);
    if (cp[i]<cp[k]) k = i;
  }
  lf->evs.mg[0] = k;
  lf->evs.fl[0] = 1.0/(2*k);
  lf->evs.fl[1] = 1-1.0/(2*k);
  nn(&lf->sp) = 0;
  fixh(&lf->sp) = 1.0/(2*k);
  startlf(des,lf,procv,0);
  ker(&lf->sp) = mk; ev(&lf->evs) = evo;
  nvm = lf->fp.nvm;
  th22 = 0;
  for (i=10; i<n-10; i++)
  { j = floor(k*datum(&lf->lfd,0,i));
    if (j>=k) j = k-1;
    dx = datum(&lf->lfd,0,i)-evptx(&lf->fp,0,j);
    if (dg1==2)
      d2 = lf->fp.coef[2*nvm+j]+dx*lf->fp.coef[3*nvm+j]+dx*dx*lf->fp.coef[4*nvm+j]/2;
    else d2 = lf->fp.coef[4*nvm+j];
    th22 += d2*d2;
  }
  hh = Wikk(mk,dg0)*sig2/th22*(n-20.0)/n;
  return(exp(log(hh)/(2*dg1+1)));
}

void rband(des,lf,hhat,meth,nmeth)
design *des;
lfit *lf;
double *hhat;
int *meth, nmeth;
{ int i, dg;
  double h0;

  /* first, estimate sigma^2 */
  dg = deg(&lf->sp); deg(&lf->sp) = 2;
  h0 = lf->sp.fixh;  lf->sp.fixh = 0.05;
printf("alp: %8.5f  h: %8.5f  deg %2d  ev %2d\n",nn(&lf->sp),fixh(&lf->sp),deg(&lf->sp),ev(&lf->evs));
  startlf(des,lf,procv,0);
  ressumm(lf,des);
  deg(&lf->sp) = dg; lf->sp.fixh = h0;
  sig2 = rv(&lf->fp); 
  printf("sd est: %8.5f\n",sqrt(sig2));

  for (i=0; i<nmeth; i++)
  { switch(meth[i])
    { case 1: hhat[i] = cp(des,lf,1);
              break;
      case 2: hhat[i] = cp(des,lf,2);
              break;
      case 3: hhat[i] = gkk(des,lf);
              break;
      case 4: hhat[i] = rsw(des,lf);
              break;
      default: hhat[i] = 0;
    }
    lf->sp.fixh = h0;
    deg(&lf->sp) = dg;
  }
}
