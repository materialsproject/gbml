/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Post-fitting functions to compute the local variance and
 *   influence functions. Also the local degrees of freedom
 *   calculations for adaptive smoothing.
 */

#include "local.h"

extern double robscale;
static double tr0, tr1, tr2;

/*
  vmat() computes (after the local fit..) the matrix 
  M2  = X^T W^2 V X.
  M12 = (X^T W V X)^{-1} M2
  Also, for convenience, tr[0] = sum(wi) tr[1] = sum(wi^2).
*/
void vmat(lfd, sp, des, M12, M2)
lfdata *lfd;
smpar *sp;
design *des;
double *M12, *M2;
{ int i, p, nk, ok;
  double link[LLEN], h, ww;
  p = des->p;
  setzero(M2,p*p);

  nk = -1;

  /* for density estimation, use integral rather than
     sum form, if W^2 is programmed...
  */
  if ((fam(sp)<=THAZ) && (link(sp)==LLOG))
  { switch(ker(sp))
    { case WGAUS: nk = WGAUS; h = des->h/SQRT2; break;
      case WRECT: nk = WRECT; h = des->h; break;
      case WEPAN: nk = WBISQ; h = des->h; break;
      case WBISQ: nk = WQUQU; h = des->h; break;
      case WTCUB: nk = W6CUB; h = des->h; break;
      case WEXPL: nk = WEXPL; h = des->h/2; break;
    }
  }

  tr0 = tr1 = 0.0;
  if (nk != -1)
  { ok = ker(sp); ker(sp) = nk;
/* compute M2 using integration. Use M12 as work matrix. */
    (des->itype)(des->xev, M2, M12, des->cf, h);
    ker(sp) = ok;
    if (fam(sp)==TDEN) multmatscal(M2,des->smwt,p*p);
    tr0 = des->ss[0];
    tr1 = M2[0]; /* n int W e^<a,A> */
  }
  else
  { for (i=0; i<des->n; i++)
    { stdlinks(link,lfd,sp,(int)des->ind[i],des->th[i],robscale);
      ww = SQR(des->w[i])*link[ZDDLL];
      tr0 += des->w[i];
      tr1 += SQR(des->w[i]);
      addouter(M2,d_xi(des,i),d_xi(des,i),p,ww);
    }
  }

  memmove(M12,M2,p*p*sizeof(double));
  for (i=0; i<p; i++)
    jacob_solve(&des->xtwx,&M12[i*p]);
}

void lf_vcov(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int i, j, k, p;
  double *M12, *M2;
  M12 = des->V; M2 = des->P; p = des->p;
  vmat(lfd,sp,des,M12,M2); /* M2 = X^T W^2 V X  tr0=sum(W) tr1=sum(W*W) */
  tr2 = m_trace(M12,p);   /* tr (XTWVX)^{-1}(XTW^2VX) */

/*
 * Covariance matrix is M1^{-1} * M2 * M1^{-1}
 * We compute this using the cholesky decomposition of
 * M2; premultiplying by M1^{-1} and squaring. This
 * is more stable than direct computation in near-singular cases.
 */
  chol_dec(M2,p,p);
  for (i=0; i<p; i++)
    for (j=0; j<i; j++)
    { M2[j*p+i] = M2[i*p+j];
      M2[i*p+j] = 0.0;
    }
  for (i=0; i<p; i++) jacob_solve(&des->xtwx,&M2[i*p]);
  for (i=0; i<p; i++)
  { for (j=0; j<p; j++)
    { M12[i*p+j] = 0;
      for (k=0; k<p; k++)
        M12[i*p+j] += M2[k*p+i]*M2[k*p+j]; /* ith column of covariance */
    }
  }
  if ((fam(sp)==TDEN) && (link(sp)==LIDENT))
    multmatscal(M12,1/SQR(des->smwt),p*p);
}

void comp_vari(lfd,sp,des,tr,t0)
lfdata *lfd;
smpar *sp;
design *des;
double *tr, *t0;
{ int i;
  lf_vcov(lfd,sp,des);
  tr[0] = tr0;
  tr[1] = tr1;
  tr[2] = tr2;
  /* influence components */
  unitvec(des->f1,0,des->p);
  jacob_solve(&des->xtwx,des->f1);
  for (i=0; i<=lfd->d; i++) t0[i] = des->f1[i];
}

/* local_df computes:
 *   tr[0] = trace(W)
 *   tr[1] = trace(W*W)
 *   tr[2] = trace( M1^{-1} M2 )
 *   tr[3] = trace( M1^{-1} M3 )
 *   tr[4] = trace( (M1^{-1} M2)^2 )
 *   tr[5] = var(theta-hat).
 */
void local_df(lfd,sp,des,tr)
lfdata *lfd;
smpar *sp;
design *des;
double *tr;
{ int i, j, p;
  double *m2, *V, ww, link[LLEN];

  tr[0] = tr[1] = tr[2] = tr[3] = tr[4] = tr[5] = 0.0;
  m2 = des->V; V = des->P; p = des->p;

  vmat(lfd,sp,des,m2,V);  /* M = X^T W^2 V X  tr0=sum(W) tr1=sum(W*W) */
  tr[0] = tr1;
  tr[1] = tr1;
  tr[2] = m_trace(m2,p);   /* tr (XTWVX)^{-1}(XTW^2VX) */

  unitvec(des->f1,0,p);
  jacob_solve(&des->xtwx,des->f1);
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
    { tr[4] += m2[i*p+j]*m2[j*p+i];  /* tr(M^2) */
      tr[5] += des->f1[i]*V[i*p+j]*des->f1[j]; /* var(thetahat) */
  }
  tr[5] = sqrt(tr[5]);

  setzero(m2,p*p);
  for (i=0; i<des->n; i++)
  { stdlinks(link,lfd,sp,(int)des->ind[i],des->th[i],robscale);
    ww = SQR(des->w[i])*des->w[i]*link[ZDDLL];
    addouter(m2,d_xi(des,i),d_xi(des,i),p,ww);
  }
  for (i=0; i<p; i++)
  { jacob_solve(&des->xtwx,&m2[i*p]);
    tr[3] += m2[i*(p+1)];
  }

  return;
}
