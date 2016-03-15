/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

/*
  preplot():  interpolates the fit to a new set of points.
  lf  -- the fit structure.
  x   -- the points to predict at.
  f   -- vector to return the predictions.
  se  -- vector to return std errors (NULL if not req'd)
  band-- char for conf band type. ('n'=none, 'g'=global etc.)
  n   -- no of predictions (or vector of margin lengths for grid)
  where -- where to predict:
           1 = points in the array x.
           2 = grid defined by margins in x.
           3 = data points from lf (ignore x).
           4 = fit points from lf (ignore x).
  what -- what to predict.
           (PCOEF etc; see lfcons.h file)

*/

static char cb;
double *sef, *fit, sigmahat;

void predptall(lf,x,what,ev,i)
lfit *lf;
double *x;
int what, ev, i;
{ double lik, rdf;
  fit[i] = dointpoint(lf,x,what,ev,i);
  if (cb=='n') return;
  sef[i] = dointpoint(lf,x,PNLX,ev,i);
  if (cb=='g')
  { sef[i] *= sigmahat;
    return;
  }
  if (cb=='l')
  { lik = dointpoint(lf,x,PLIK,ev,i);
    rdf = dointpoint(lf,x,PRDF,ev,i);
    sef[i] *= sqrt(-2*lik/rdf);
    return;
  }
  if (cb=='p')
  { sef[i] = sigmahat*sqrt(1+sef[i]*sef[i]);
    return;
  }
}

void prepvector(lf,x,n,what) /* interpolate a vector */
lfit *lf;
double **x;
int n, what;
{ int i, j;
  double xx[MXDIM];
  for (i=0; i<n; i++)
  { for (j=0; j<lf->fp.d; j++) xx[j] = x[j][i];
    predptall(lf,xx,what,ev(&lf->evs),i);
    if (lf_error) return;
  }
}

void prepfitp(lf,what)
lfit *lf;
int what;
{ int  i;
  for (i=0; i<lf->fp.nv; i++)
  { predptall(lf,evpt(&lf->fp,i),what,EFITP,i);
    if (lf_error) return;
  }
}

void prepgrid(lf,x,mg,n,what) /* interpolate a grid given margins */
lfit *lf;
double **x;
Sint *mg;
int n, what;
{ int i, ii, j, d;
  double xv[MXDIM];
  d = lf->fp.d;
  for (i=0; i<n; i++)
  { ii = i;
    for (j=0; j<d; j++)
    { xv[j] = x[j][ii%mg[j]];
      ii /= mg[j];
    }
    predptall(lf,xv,what,ev(&lf->evs),i);
    if (lf_error) return;
  }
}

void preplot(lf,x,f,se,band,mg,where,what)
lfit *lf;
double **x, *f, *se;
Sint *mg;
int where, what;
char band;
{ int d, i, n=0;
  double *xx[MXDIM];
  d = lf->fp.d;
  fit = f;
  sef = se;
  cb = band;
  if (cb!='n') sigmahat = sqrt(rv(&lf->fp));

//printf("      where = %d\n", where);
//printf("      what  = %d\n", what );

  switch(where)
  { case 1: /* vector */
      n = mg[0];
      prepvector(lf,x,n,what);
      break;
    case 2: /* grid */
      n = 1;
      for (i=0; i<d; i++) n *= mg[i];
      prepgrid(lf,x,mg,n,what);
      break;
    case 3: /* data */
      n = lf->lfd.n;
      if ((ev(&lf->evs)==EDATA) | (ev(&lf->evs)==ECROS))
        prepfitp(lf,what);
      else
      { for (i=0; i<d; i++) xx[i] = dvari(&lf->lfd,i);
        prepvector(lf,xx,n,what);
      }
      break;
    case 4: /* fit points */
      n = lf->fp.nv;
      prepfitp(lf,what);
      break;
    default:
      ERROR(("unknown where in preplot"));
  }

  if ((what==PT0)|(what==PVARI))
    for (i=0; i<n; i++) f[i] = f[i]*f[i];
}
