/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern int lf_status;
static double u[MXDIM], ilim[2*MXDIM], *ff, hh, *cff;
static lfdata *den_lfd;
static design *den_des;
static smpar *den_sp;
int fact[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
int de_mint  = 20;
int de_itype = IDEFA;
int de_renorm= 0;

int multint(), prodint(), gausint(), mlinint();

void prresp(coef,resp,p)
double *coef, *resp;
int p;
{ int i, j;
  printf("Coefficients:\n");
  for (i=0; i<p; i++) printf("%8.5f ",coef[i]);
  printf("\n");
  printf("Response matrix:\n");
  for (i=0; i<p; i++)
  { for (j=0; j<p; j++) printf("%9.6f, ",resp[i+j*p]);
    printf("\n");
  }
}

int mif(u,d,resp,M)
double *u, *resp, *M;
int d;
{ double wt;
  int i, j, p;

  p = den_des->p;
  wt = weight(den_lfd, den_sp, u, NULL, hh, 0, 0.0);
  if (wt==0)
  { setzero(resp,p*p);
    return(p*p);
  }

  fitfun(den_lfd, den_sp, u,NULL,ff,NULL);
  if (link(den_sp)==LLOG)
    wt *= lf_exp(innerprod(ff,cff,p));
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      resp[i*p+j] = wt*ff[i]*ff[j];
  return(p*p);
}

int multint(t,resp1,resp2,cf,h)
double *t, *resp1, *resp2, *cf, h;
{ int d, i, mg[MXDIM];

  if (ker(den_sp)==WGAUS) return(gausint(t,resp1,resp2,cf,h,den_lfd->sca));

  d = den_lfd->d;
  for (i=0; i<d; i++) mg[i] = de_mint;

  hh = h;
  cff= cf;
  simpsonm(mif,ilim,&ilim[d],d,resp1,mg,resp2);
  return(LF_OK);
}

int mlinint(t,resp1,resp2,cf,h)
double *t, *resp1, *resp2, *cf, h;
{
  double hd, nb, wt, wu, g[4], w0, w1, v, *sca;
  int d, p, i, j, jmax, k, l, z, jj[2];

  d = den_lfd->d; p = den_des->p; sca = den_lfd->sca;
  hd = 1;
  for (i=0; i<d; i++) hd *= h*sca[i];

  if (link(den_sp)==LIDENT)
  { setzero(resp1,p*p);
    resp1[0] = wint(d,NULL,0,ker(den_sp))*hd;
    if (deg(den_sp)==0) return(LF_OK);
    jj[0] = 2; w0 = wint(d,jj,1,ker(den_sp))*hd*h*h;
    for (i=0; i<d; i++) resp1[(i+1)*p+i+1] = w0*sca[i]*sca[i];
    if (deg(den_sp)==1) return(LF_OK);
    for (i=0; i<d; i++)
    { j = p-(d-i)*(d-i+1)/2;
      resp1[j] = resp1[p*j] = w0*sca[i]*sca[i]/2;
    }
    if (d>1)
    { jj[1] = 2;
      w0 = wint(d,jj,2,ker(den_sp)) * hd*h*h*h*h;
    }
    jj[0] = 4;
    w1 = wint(d,jj,1,ker(den_sp)) * hd*h*h*h*h/4;
    z = d+1;
    for (i=0; i<d; i++)
    { k = p-(d-i)*(d-i+1)/2;
      for (j=i; j<d; j++)
      { l = p-(d-j)*(d-j+1)/2;
        if (i==j) resp1[z*p+z] = w1*SQR(sca[i])*SQR(sca[i]);
        else
        { resp1[z*p+z] = w0*SQR(sca[i])*SQR(sca[j]);
          resp1[k*p+l] = resp1[k+p*l] = w0/4*SQR(sca[i])*SQR(sca[j]);
        }
        z++;
    } }
    return(LF_OK);
  }
  switch(deg(den_sp))
  { case 0:
      resp1[0] = lf_exp(cf[0])*wint(d,NULL,0,ker(den_sp))*hd;
      return(LF_OK);
    case 1:
      nb = 0.0;
      for (i=1; i<=d; i++)
      { v = h*cf[i]*sca[i-1];
        nb += v*v;
      }
      if (ker(den_sp)==WGAUS)
      { w0 = 1/(GFACT*GFACT);
        g[0] = lf_exp(cf[0]+w0*nb/2+d*log(S2PI/2.5));
        g[1] = g[3] = g[0]*w0;
        g[2] = g[0]*w0*w0;
      }
      else
      { wt = wu = lf_exp(cf[0]);
        w0 = wint(d,NULL,0,ker(den_sp)); g[0] = wt*w0;
        g[1] = g[2] = g[3] = 0.0;
        j = 0; jmax = (d+2)*de_mint;
        while ((j<jmax) && (wt*w0/g[0]>1.0e-8))
        { j++;
          jj[0] = 2*j; w0 = wint(d,jj,1,ker(den_sp));
          if (d==1) g[3] += wt * w0;
          else
          { jj[0] = 2; jj[1] = 2*j-2; w1 = wint(d,jj,2,ker(den_sp));
            g[3] += wt*w1;
            g[2] += wu*(w0-w1);
          }
          wt /= (2*j-1.0); g[1] += wt*w0;
          wt *= nb/(2*j); g[0] += wt*w0;
          wu /= (2*j-1.0)*(2*j);
          if (j>1) wu *= nb;
        }
        if (j==jmax) WARN(("mlinint: series not converged"));
      }
      g[0] *= hd; g[1] *= hd;
      g[2] *= hd; g[3] *= hd;
      resp1[0] = g[0];
      for (i=1; i<=d; i++)
      { resp1[i] = resp1[(d+1)*i] = cf[i]*SQR(h*sca[i-1])*g[1];
        for (j=1; j<=d; j++)
        { resp1[(d+1)*i+j] = (i==j) ? g[3]*SQR(h*sca[i-1]) : 0;
          resp1[(d+1)*i+j] += g[2]*SQR(h*h*sca[i-1]*sca[j-1])*cf[i]*cf[j];
        }
      }
      return(LF_OK);
  }
  ERROR(("mlinint: deg=0,1 only"));
  return(LF_ERR);
}

void prodintresp(resp,prod_wk,dim,deg,p)
double *resp, prod_wk[MXDIM][2*MXDEG+1];
int dim, deg, p;
{ double prod;
  int i, j, k, j1, k1;

  prod = 1.0;
  for (i=0; i<dim; i++) prod *= prod_wk[i][0];
  resp[0] += prod;
  if (deg==0) return;

  for (j1=1; j1<=deg; j1++)
  { for (j=0; j<dim; j++)
    { prod = 1.0;
      for (i=0; i<dim; i++) prod *= prod_wk[i][j1*(j==i)];
      prod /= fact[j1];
      resp[1 + (j1-1)*dim +j] += prod;
    }
  }

  for (k1=1; k1<=deg; k1++)
    for (j1=k1; j1<=deg; j1++)
    { for (k=0; k<dim; k++)
        for (j=0; j<dim; j++)
        { prod = 1.0;
          for (i=0; i<dim; i++) prod *= prod_wk[i][k1*(k==i) + j1*(j==i)];
          prod /= fact[k1]*fact[j1];
          resp[ (1+(k1-1)*dim+k)*p + 1+(j1-1)*dim+j] += prod;
        }
    }
}

int prodint(t,resp,resp2,coef,h)
double *t, *resp, *resp2, *coef, h;
{ int dim, p, i, j, k, st=0;
  double cf[MXDEG+1], hj, hs, prod_wk[MXDIM][2*MXDEG+1];

  dim = den_lfd->d;
  p = den_des->p;
  for (i=0; i<p*p; i++) resp[i] = 0.0;
  cf[0] = coef[0];

/*  compute the one dimensional terms
 */
  for (i=0; i<dim; i++)
  { hj = 1; hs = h*den_lfd->sca[i];
    for (j=0; j<deg(den_sp); j++)
    { hj *= hs;
      cf[j+1] = hj*coef[ j*dim+i+1 ];
    }
    st = onedint(den_sp,cf,ilim[i]/hs,ilim[i+dim]/hs,prod_wk[i]);
    if (st==LF_BADP) return(st);
    hj = 1;
    for (j=0; j<=2*deg(den_sp); j++)
    { hj *= hs;
      prod_wk[i][j] *= hj;
    }
    cf[0] = 0.0; /* so we only include it once, when d>=2 */
  }

/*  transfer to the resp array
 */
  prodintresp(resp,prod_wk,dim,deg(den_sp),p);

/* Symmetrize.
*/
  for (k=0; k<p; k++)
    for (j=k; j<p; j++)
      resp[j*p+k] = resp[k*p+j];

  return(st);
}

int gausint(t,resp,C,cf,h,sca)
double *t, *resp, *C, *cf, h, *sca;
{ double nb, det, z, *P;
  int d, p, i, j, k, l, m1, m2, f;
  d = den_lfd->d; p = den_des->p;
  m1 = d+1; nb = 0;
  P = &C[d*d];
  resp[0] = 1;
  for (i=0; i<d; i++)
  { C[i*d+i] = SQR(GFACT/(h*sca[i]))-cf[m1++];
    for (j=i+1; j<d; j++) C[i*d+j] = C[j*d+i] = -cf[m1++];
  }
  eig_dec(C,P,d);
  det = 1;
  for (i=1; i<=d; i++)
  { det *= C[(i-1)*(d+1)];
    if (det <= 0) return(LF_BADP);
    resp[i] = cf[i];
    for (j=1; j<=d; j++) resp[j+i*p] = 0;
    resp[i+i*p] = 1;
    svdsolve(&resp[i*p+1],u,P,C,P,d,0.0);
  }
  svdsolve(&resp[1],u,P,C,P,d,0.0);
  det = sqrt(det);
  for (i=1; i<=d; i++)
  { nb += cf[i]*resp[i];
    resp[i*p] = resp[i];
    for (j=1; j<=d; j++)
      resp[i+p*j] += resp[i]*resp[j];
  }
  m1 = d;
  for (i=1; i<=d; i++)
    for (j=i; j<=d; j++)
    { m1++; f = 1+(i==j);
      resp[m1] = resp[m1*p] = resp[i*p+j]/f;
      m2 = d;
      for (k=1; k<=d; k++)
      { resp[m1+k*p] = resp[k+m1*p] =
        ( resp[i]*resp[j*p+k] + resp[j]*resp[i*p+k]
        + resp[k]*resp[i*p+j] - 2*resp[i]*resp[j]*resp[k] )/f;
        for (l=k; l<=d; l++)
        { m2++; f = (1+(i==j))*(1+(k==l));
          resp[m1+m2*p] = resp[m2+m1*p] = ( resp[i+j*p]*resp[k+l*p]
            + resp[i+k*p]*resp[j+l*p] + resp[i+l*p]*resp[j+k*p]
            - 2*resp[i]*resp[j]*resp[k]*resp[l] )/f;
    } } }
  z = lf_exp(d*0.918938533+cf[0]+nb/2)/det;
  multmatscal(resp,z,p*p);
  return(LF_OK);
}

int likeden(coef, lk0, f1, A)
double *coef, *lk0, *f1, *A;
{ double lk=0.0, r;
  int i, j, p, rstat;

  lf_status = LF_OK;
  p = den_des->p;
  if ((link(den_sp)==LIDENT) && (coef[0] != 0.0)) return(NR_BREAK);
  lf_status = (den_des->itype)(den_des->xev,A,den_des->xtwx.Q,coef,den_des->h);
  if (lf_error) lf_status = LF_ERR;
  if (lf_status==LF_BADP)
  { *lk0 = -1.0e300;
    return(NR_REDUCE);
  }
  if (lf_status!=LF_OK) return(NR_BREAK);
  if (lf_debug>2) prresp(coef,A,p);

  den_des->xtwx.p = p;
  rstat = NR_OK;
  switch(link(den_sp))
  { case LLOG:
      r = den_des->ss[0]/A[0];
      coef[0] += log(r);
      multmatscal(A,r,p*p);
      A[0] = den_des->ss[0];
      lk = -A[0];
      if (fabs(coef[0]) > 700)
      { lf_status = LF_OOB;
        rstat = NR_REDUCE;
      }
      for (i=0; i<p; i++)
      { lk += coef[i]*den_des->ss[i];
        f1[i] = den_des->ss[i]-A[i];
      }
      break;
    case LIDENT:
      lk = 0.0;
      for (i=0; i<p; i++)
      { f1[i] = den_des->ss[i];
        for (j=0; j<p; j++)
          den_des->res[i] -= A[i*p+j]*coef[j];
      }
      break;
  }
  *lk0 = den_des->llk = lk;

  return(rstat);
}

int inre(x,bound,d)
double *x, *bound;
int d;
{ int i, z;
  z = 1;
  for (i=0; i<d; i++)
    if (bound[i]<bound[i+d])
      z &= (x[i]>=bound[i]) & (x[i]<=bound[i+d]);
  return(z);
}

int setintlimits(lfd, x, h, ang, lset)
lfdata *lfd;
int *ang, *lset;
double *x, h;
{ int d, i;
  d = lfd->d;
  *ang = *lset = 0;
  for (i=0; i<d; i++)
  { if (lfd->sty[i]==STANGL)
    { ilim[i+d] = ((h<2) ? 2*asin(h/2) : PI)*lfd->sca[i];
      ilim[i] = -ilim[i+d];
      *ang = 1;
    }
    else
    { ilim[i+d] = h*lfd->sca[i];
      ilim[i] = -ilim[i+d];

      if (lfd->sty[i]==STLEFT) { ilim[i+d] = 0; *lset = 1; }
      if (lfd->sty[i]==STRIGH) { ilim[i] = 0;   *lset = 1; }

      if (lfd->xl[i]<lfd->xl[i+d]) /* user limits for this variable */
      { if (lfd->xl[i]-x[i]> ilim[i])
        { ilim[i] = lfd->xl[i]-x[i]; *lset=1; }
        if (lfd->xl[i+d]-x[i]< ilim[i+d])
        { ilim[i+d] = lfd->xl[i+d]-x[i]; *lset=1; }
      }
    }
    if (ilim[i]==ilim[i+d]) return(LF_DEMP); /* empty integration */
  }
  return(LF_OK);
}

int selectintmeth(itype,lset,ang)
int itype, lset, ang;
{
  if (itype==IDEFA) /* select the default method */
  { if (fam(den_sp)==THAZ)
    { if (ang) return(IDEFA);
      return( IHAZD );
    }

    if (ubas(den_sp)) return(IMULT);

    if (ang) return(IMULT);

    if (iscompact(ker(den_sp)))
    { if (kt(den_sp)==KPROD) return(IPROD);
      if (lset)
        return( (den_lfd->d==1) ? IPROD : IMULT );
      if (deg(den_sp)<=1) return(IMLIN);
      if (den_lfd->d==1) return(IPROD);
      return(IMULT);
    }

    if (ker(den_sp)==WGAUS)
    { if (lset) WARN(("Integration for Gaussian weights ignores limits"));
      if ((den_lfd->d==1)|(kt(den_sp)==KPROD)) return(IPROD);
      if (deg(den_sp)<=1) return(IMLIN);
      if (deg(den_sp)==2) return(IMULT);
    }

    return(IDEFA);
  }

  /* user provided an integration method, check it is valid */

  if (fam(den_sp)==THAZ)
  { if (ang) return(INVLD);
    if (!iscompact(ker(den_sp))) return(INVLD);
    return( ((kt(den_sp)==KPROD) | (kt(den_sp)==KSPH)) ? IHAZD : INVLD );
  }

  if ((ang) && (itype != IMULT)) return(INVLD);

  switch(itype)
  { case IMULT:
      if (ker(den_sp)==WGAUS) return(deg(den_sp)==2);
      return( iscompact(ker(den_sp)) ? IMULT : INVLD );
    case IPROD: return( ((den_lfd->d==1) | (kt(den_sp)==KPROD)) ? IPROD : INVLD );
    case IMLIN: return( ((kt(den_sp)==KSPH) && (!lset) &&
      (deg(den_sp)<=1)) ? IMLIN : INVLD );
  }

  return(INVLD);
}

int densinit(lfd,des,sp,cf)
lfdata *lfd;
design *des;
smpar *sp;
double *cf;
{ int p, i, ii, j, nnz, rnz, ang, lset, status;
  double w;

  den_lfd = lfd;
  den_des = des;
  den_sp  = sp;

  p = des->p;
  ff = des->xtwx.wk;
  cf[0] = NOSLN;
  for (i=1; i<p; i++) cf[i] = 0.0;

  if (!inre(des->xev,lfd->xl,lfd->d)) return(LF_XOOR);

  status = setintlimits(lfd,des->xev,des->h,&ang,&lset);
  if (status != LF_OK) return(status);

  switch(selectintmeth(de_itype,lset,ang))
  { case IMULT: des->itype = multint; break;
    case IPROD: des->itype = prodint; break;
    case IMLIN: des->itype = mlinint; break;
    case IHAZD: des->itype = hazint; break;
    case INVLD: ERROR(("Invalid integration method %d",de_itype));
                break;
    case IDEFA: ERROR(("No integration type available for this model"));
                break;
    default: ERROR(("densinit: unknown integral type"));
  }

  switch(deg(den_sp))
  { case 0: rnz = 1; break;
    case 1: rnz = 1; break;
    case 2: rnz = lfd->d+1; break;
    case 3: rnz = lfd->d+2; break;
    default: ERROR(("densinit: invalid degree %d",deg(den_sp)));
  }
  if (lf_error) return(LF_ERR);

  setzero(des->ss,p);
  nnz = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    if (!cens(lfd,ii))
    { w = des->w[i]*prwt(lfd,ii);
      for (j=0; j<p; j++) des->ss[j] += d_xij(des,i,j)*w;
      if (des->w[i]>0.00001) nnz++;
  } }

  if (fam(den_sp)==THAZ) haz_init(lfd,des,sp,ilim);

  if (lf_debug>2)
  { printf("    LHS: ");
    for (i=0; i<p; i++) printf(" %8.5f",des->ss[i]);
    printf("\n");
  }

  switch(link(den_sp))
  { case LIDENT:
      cf[0] = 0.0;
      return(LF_OK);
    case LLOG:
      if (nnz<rnz) { cf[0] = -1000; return(LF_DNOP); }
      cf[0] = 0.0;
      return(LF_OK);
    default:
      ERROR(("unknown link in densinit"));
      return(LF_ERR);
  }
}
