#include "local.h"
static double scb_crit, *x, c[10], kap[5], kaq[5], max_p2;
/* static int side, type; */
static int type;
design *scb_des;

double covar_par(lf,des,x1,x2)
lfit *lf;
design *des;
double x1, x2;
{ double *v1, *v2, *wk;
  paramcomp *pc;
  int i, j, p, ispar;

  v1 = des->f1; v2 = des->ss; wk = des->oc;
  ispar = (ker(&lf->sp)==WPARM) && (haspc(&lf->pc));
  p = npar(&lf->sp);

/*  for parametric models, the covariance is
 *  A(x1)^T (X^T W V X)^{-1} A(x2)
 *  which we can find easily from the parametric component.
 */
  if (ispar)
  { pc = &lf->pc;
    fitfun(&lf->lfd, &lf->sp, &x1,pc->xbar,v1,NULL);
    fitfun(&lf->lfd, &lf->sp, &x2,pc->xbar,v2,NULL);
    jacob_hsolve(&lf->pc.xtwx,v1);
    jacob_hsolve(&lf->pc.xtwx,v2);
  }

/*  for non-parametric models, we must use the cholseky decomposition
 *  of M2 = X^T W^2 V X. Courtesy of comp_vari, we already have
 *  des->P = M2^{1/2} M1^{-1}.
 */
  if (!ispar)
  { fitfun(&lf->lfd, &lf->sp, &x1,des->xev,wk,NULL);
    for (i=0; i<p; i++)
    { v1[i] = 0;
      for (j=0; j<p; j++) v1[i] += des->P[i*p+j]*wk[j];
    }
    fitfun(&lf->lfd, &lf->sp, &x2,des->xev,wk,NULL);
    for (i=0; i<p; i++)
    { v2[i] = 0;
      for (j=0; j<p; j++) v2[i] += des->P[i*p+j]*wk[j];
    }
  }

  return(innerprod(v1,v2,p));
}

void cumulant(lf,des,sd)
lfit *lf;
design *des;
double sd;
{ double b2i, b3i, b3j, b4i;
  double ss, si, sj, uii, uij, ujj, k1;
  int ii, i, j, jj;
  for (i=1; i<10; i++) c[i] = 0.0;
  k1 = 0;

  /* ss = sd*sd; */
  ss = covar_par(lf,des,des->xev[0],des->xev[0]);

/*
 * this isn't valid for nonparametric models. At a minimum,
 * the sums would have to include weights. Still have to work
 * out the right way.
 */
  for (i=0; i<lf->lfd.n; i++)
  { ii = des->ind[i];
    b2i = b2(des->th[i],fam(&lf->sp),prwt(&lf->lfd,ii));
    b3i = b3(des->th[i],fam(&lf->sp),prwt(&lf->lfd,ii));
    b4i = b4(des->th[i],fam(&lf->sp),prwt(&lf->lfd,ii));
    si = covar_par(lf,des,des->xev[0],datum(&lf->lfd,0,ii));
    uii= covar_par(lf,des,datum(&lf->lfd,0,ii),datum(&lf->lfd,0,ii));
    if (lf_error) return;

    c[2] += b4i*si*si*uii;
    c[6] += b4i*si*si*si*si;
    c[7] += b3i*si*uii;
    c[8] += b3i*si*si*si;
    /* c[9] += b2i*si*si*si*si;
       c[9] += b2i*b2i*si*si*si*si; */
    k1 += b3i*si*(si*si/ss-uii);

    /* i=j components */
    c[1] += b3i*b3i*si*si*uii*uii;
    c[3] += b3i*b3i*si*si*si*si*uii;
    c[4] += b3i*b3i*si*si*uii*uii;

    for (j=i+1; j<lf->lfd.n; j++)
    { jj = des->ind[j];
      b3j = b3(des->th[j],fam(&lf->sp),prwt(&lf->lfd,jj));
      sj = covar_par(lf,des,des->xev[0],datum(&lf->lfd,0,jj));
      uij= covar_par(lf,des,datum(&lf->lfd,0,ii),datum(&lf->lfd,0,jj));
      ujj= covar_par(lf,des,datum(&lf->lfd,0,jj),datum(&lf->lfd,0,jj));

      c[1] += 2*b3i*b3j*si*sj*uij*uij;
      c[3] += 2*b3i*b3j*si*si*sj*sj*uij;
      c[4] += b3i*b3j*uij*(si*si*ujj+sj*sj*uii);
      if (lf_error) return;
    }
  }
  c[5] = c[1];
  c[7] = c[7]*c[8];
  c[8] = c[8]*c[8];

  c[1] /= ss; c[2] /= ss; c[3] /= ss*ss; c[4] /= ss;
  c[5] /= ss; c[6] /= ss*ss; c[7] /= ss*ss;
  c[8] /= ss*ss*ss; c[9] /= ss*ss;

/* constants used in p(x,z) computation */
  kap[1] = k1/(2*sqrt(ss));
  kap[2] = 1 + 0.5*(c[1]-c[2]+c[4]-c[7]) - 3*c[3] + c[6] + 1.75*c[8];
  kap[4] = -9*c[3] + 3*c[6] + 6*c[8] + 3*c[9];

/* constants used in q(x,u) computation */
  kaq[2] = c[3] - 1.5*c[8] - c[5] - c[4] + 0.5*c[7] + c[6] - c[2];
  kaq[4] = -3*c[3] - 6*c[4] - 6*c[5] + 3*c[6] + 3*c[7] - 3*c[8] + 3*c[9];
}

/* q2(u) := u+q2(x,u) in paper */
double q2(u)
double u;
{ return(u-u*(36.0*kaq[2] + 3*kaq[4]*(u*u-3) + c[8]*((u*u-10)*u*u+15))/72.0);
}

/*  p2(u) := p2(x,u) in paper */
double p2(u)
double u;
{ return( -u*( 36*(kap[2]-1+kap[1]*kap[1])
     + 3*(kap[4]+4*kap[1]*sqrt(kap[3]))*(u*u-3)
     + c[8]*((u*u-10)*u*u+15) ) / 72 );
}

extern int likereg();
double gldn_like(a)
double a;
{ int err;

  scb_des->fix[0] = 1;
  scb_des->cf[0] = a;
  max_nr(likereg, scb_des->cf, scb_des->oc, scb_des->res, scb_des->f1,
    &scb_des->xtwx, scb_des->p, lf_maxit, 1.0e-6, &err); 
  scb_des->fix[0] = 0;

  return(scb_des->llk);
}

/* v1/v2 is correct for deg=0 only */
void get_gldn(fp,des,lo,hi,v)
fitpt *fp;
design *des;
double *lo, *hi;
int v;
{ double v1, v2, c, tlk;
  int err;

  v1 = fp->nlx[v];
  v2 = fp->t0[v];
  c = scb_crit * v1 / v2;
  tlk = des->llk - c*c/2;
printf("v %8.5f %8.5f  c %8.5f  tlk %8.5f   llk %8.5f\n",v1,v2,c,tlk,des->llk);

  /* want: { a : l(a) >= l(a-hat) - c*c/2 } */
  lo[v] = fp->coef[v] - scb_crit*v1;
  hi[v] = fp->coef[v] + scb_crit*v1;

  err = 0;

printf("lo %2d\n",v);
  lo[v] = solve_secant(gldn_like,tlk,lo[v],fp->coef[v],1e-8,BDF_EXPLEFT,&err);
  if (err>0) printf("solve_secant error\n");
printf("hi %2d\n",v);
  hi[v] = solve_secant(gldn_like,tlk,fp->coef[v],hi[v],1e-8,BDF_EXPRIGHT,&err);
  if (err>0) printf("solve_secant error\n");
}

int procvscb2(des,lf,v)
design *des;
lfit *lf;
int v;
{ double thhat, sd, *lo, *hi, u;
  int err, st, tmp;
  x = des->xev = evpt(&lf->fp,v);
  tmp = haspc(&lf->pc);
  /* if ((ker(&lf->sp)==WPARM) && (haspc(&lf->pc)))
  { lf->coef[v] = thhat = addparcomp(lf,des->xev,PCOEF);
    lf->nlx[v] = lf->t0[v] = sd = addparcomp(lf,des->xev,PNLX);
  }
  else */
  { haspc(&lf->pc) = 0;
    st = procv(des,lf,v);
    thhat = lf->fp.coef[v];
    sd = lf->fp.nlx[v];
  }
  if ((type==GLM2) | (type==GLM3) | (type==GLM4))
  { if (ker(&lf->sp) != WPARM)
      WARN(("nonparametric fit; correction is invalid"));
    cumulant(lf,des,sd);
  }
  haspc(&lf->pc) = tmp;
  lo = lf->fp.L;
  hi = &lo[lf->fp.nvm];
  switch(type)
  {
    case GLM1:
      return(st);
    case GLM2: /* centered scr */
      lo[v] = kap[1];
      hi[v] = sqrt(kap[2]);
      return(st);
    case GLM3: /* corrected 2 */
      lo[v] = solve_secant(q2,scb_crit,0.0,2*scb_crit,0.000001,BDF_NONE,&err);
      return(st);
    case GLM4: /* corrected 2' */
      u = fabs(p2(scb_crit));
      max_p2 = MAX(max_p2,u);
      return(st);
    case GLDN:
      get_gldn(&lf->fp,des,lo,hi,v);
      return(st);
  }
  ERROR(("procvscb2: invalid type"));
  return(st);
}

void scb(des,lf)
design *des;
lfit *lf;
{ double k1, k2; /* kap[10], */
  double *lo, *hi, sig, thhat, nlx;
  int i, nterms;

  scb_des= des;

  npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);
  des_init(des,lf->lfd.n,npar(&lf->sp));
  link(&lf->sp) = defaultlink(link(&lf->sp),fam(&lf->sp));

  type = geth(&lf->fp);

  if (type >= 80) /* simultaneous */
  {
    nterms = constants(des,lf);
    scb_crit = critval(0.05,lf->fp.kap,nterms,lf->lfd.d,TWO_SIDED,0.0,GAUSS);
    type -= 10;
  }
  else /* pointwise */
  { lf->fp.kap[0] = 1;
    scb_crit = critval(0.05,lf->fp.kap,1,lf->lfd.d,TWO_SIDED,0.0,GAUSS);
  }

  max_p2 = 0.0;
  startlf(des,lf,procvscb2,0);

  if ((fam(&lf->sp)&64)==64)
  { i = haspc(&lf->pc); haspc(&lf->pc) = 0;
    ressumm(lf,des);
    haspc(&lf->pc) = i;
    sig = sqrt(rv(&lf->fp));
  }
  else sig = 1.0;

  lo = lf->fp.L;
  hi = &lo[lf->fp.nvm];
  for (i=0; i<lf->fp.nv; i++)
  { thhat = lf->fp.coef[i];
    nlx = lf->fp.nlx[i];
    switch(type)
    {
      case GLM1:  /* basic scb */
        lo[i] = thhat - scb_crit * sig * nlx;
        hi[i] = thhat + scb_crit * sig * nlx;
        break;
      case GLM2:
        k1 = lo[i];
        k2 = hi[i];
        lo[i] = thhat - k1*nlx - scb_crit*nlx*k2;
        hi[i] = thhat - k1*nlx + scb_crit*nlx*k2;
        break;
      case GLM3:
        k1 = lo[i];
        lo[i] = thhat - k1*nlx;
        hi[i] = thhat + k1*nlx;
      case GLM4:  /* corrected 2' */
        lo[i] = thhat - (scb_crit-max_p2)*lf->fp.nlx[i];
        hi[i] = thhat + (scb_crit-max_p2)*lf->fp.nlx[i];
        break;
      case GLDN:
        break;
    }
  }
}
