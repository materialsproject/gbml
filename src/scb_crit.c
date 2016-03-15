/*
 *   Copyright (c) 1996-2004 Catherine Loader.
 *
 *   Computes the critical values from constants kappa0 etc
 *   and significance level.
 */

#include <math.h>
#include "local.h"
#include "tube.h"

/*
 * some old math libraries choke on lgamma()...
 */
/* #define LGAMMA(arg) lgamma(arg) */
#define LOGPI 1.144729885849400174143427

/* area(d) = 2 pi^(d/2) / Gamma(d/2)
 *         = surface area of unit sphere in R^d
 */
static double A[10] =
  { 1, /* d=0, whatever */
    2,
    6.2831853071795864770, /* 2*pi */
    12.566370614359172954, /* 4*pi */
    19.739208802178717238, /* 2*pi^2 */
    26.318945069571622985, /* 8/3*pi^2 */
    31.006276680299820177, /* pi^3 */
    33.073361792319808190, /* 16/15*pi^3 */
    32.469697011334145747, /* 1/3*pi^4 */
    29.686580124648361825  /* 32/105*pi^4 */
  };

double area(d)
int d;
{ if (d<10) return(A[d]);
  return(2*exp(d*LOGPI/2.0-LGAMMA(d/2.0)));
}

double tailp_uniform(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * ibeta(1-c*c,(n-d+i-1)/2.0,(d+1-i)/2.0) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double tailp_gaussian(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * (1-pchisq(c*c,(double) d+1-i)) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double tailp_tprocess(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * (1-pf(c*c/(d+1-i),(double) d+1-i, n)) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double taild_uniform(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * 2*c*dbeta(1-c*c,(n-d+i-1)/2.0,(d+1-i)/2.0,0) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double taild_gaussian(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * 2*c*dchisq(c*c,(double) d+1-i,0) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double taild_tprocess(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * 2*c*df(c*c/(d+1-i),(double) d+1-i, n,0) / ((d+1-i)*area(d+1-i));
  return( (s==TWO_SIDED) ? 2*p : p );
}

double tailp(c,k0,m,d,s,nu, process)
double c, *k0, nu;
int m, d, s, process;
{ switch(process)
  { case UNIF:  return(tailp_uniform(c,k0,m,d,s,nu));
    case GAUSS: return(tailp_gaussian(c,k0,m,d,s,nu));
    case TPROC: return(tailp_tprocess(c,k0,m,d,s,nu));
  }
  printf("taild: unknown process.\n");
  return(0.0);
}

double taild(c,k0,m,d,s,nu, process)
double c, *k0, nu;
int m, d, s, process;
{ switch(process)
  { case UNIF:  return(taild_uniform(c,k0,m,d,s,nu));
    case GAUSS: return(taild_gaussian(c,k0,m,d,s,nu));
    case TPROC: return(taild_tprocess(c,k0,m,d,s,nu));
  }
  printf("taild: unknown process.\n");
  return(0.0);
}

double critval(alpha,k0,m,d,s,nu,process)
double alpha, *k0, nu;
int m, d, s, process;
{ double c, cn, c0, c1, tp, td;
  int j, maxit;
  double (*tpf)(), (*tdf)();

  maxit = 20;
  if (m<0)
  { printf("critval: no terms?\n");
    return(2.0);
  }
  if (m>d+1) m = d+1;
  if ((alpha<=0) | (alpha>=1))
  { printf("critval: invalid alpha %8.5f\n",alpha);
    return(2.0);
  }
  if (alpha>0.5)
    printf("critval: A mighty large tail probability alpha=%8.5f\n",alpha);
  if (m==0) { d = 0; k0[0] = 1; m = 1; }

  switch(process)
  { case UNIF:
      c = 0.5; c0 = 0.0; c1 = 1.0;
      tpf = tailp_uniform;
      tdf = taild_uniform;
      break;
    case GAUSS:
      c = 2.0; c0 = 0.0; c1 = 0.0;
      tpf = tailp_gaussian;
      tdf = taild_gaussian;
      break;
    case TPROC:
      c = 2.0; c0 = 0.0; c1 = 0.0;
      tpf = tailp_tprocess;
      tdf = taild_tprocess;
      break;
    default:
      printf("critval: unknown process.\n");
      return(0.0);
  }

  for (j=0; j<maxit; j++)
  { tp = tpf(c,k0,m,d,s,nu)-alpha;
    td = tdf(c,k0,m,d,s,nu);
    if (tp>0) c0 = c;
    if (tp<0) c1 = c;
    cn = c + tp/td;
    if (cn<c0) cn = (c+c0)/2;
    if ((c1>0.0) && (cn>c1)) cn = (c+c1)/2;
    c = cn;
    if (fabs(tp/alpha)<1.0e-10) return(c);
  }
  return(c);
}
