/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include "mutil.h"

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) */

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508 /* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
static double sferr_halves[31] = {
0.0, /* n=0 - wrong, place holder only */
0.1534264097200273452913848,  /* 0.5 */
0.0810614667953272582196702,  /* 1.0 */
0.0548141210519176538961390,  /* 1.5 */
0.0413406959554092940938221,  /* 2.0 */
0.03316287351993628748511048, /* 2.5 */
0.02767792568499833914878929, /* 3.0 */
0.02374616365629749597132920, /* 3.5 */
0.02079067210376509311152277, /* 4.0 */
0.01848845053267318523077934, /* 4.5 */
0.01664469118982119216319487, /* 5.0 */
0.01513497322191737887351255, /* 5.5 */
0.01387612882307074799874573, /* 6.0 */
0.01281046524292022692424986, /* 6.5 */
0.01189670994589177009505572, /* 7.0 */
0.01110455975820691732662991, /* 7.5 */
0.010411265261972096497478567, /* 8.0 */
0.009799416126158803298389475, /* 8.5 */
0.009255462182712732917728637, /* 9.0 */
0.008768700134139385462952823, /* 9.5 */
0.008330563433362871256469318, /* 10.0 */
0.007934114564314020547248100, /* 10.5 */
0.007573675487951840794972024, /* 11.0 */
0.007244554301320383179543912, /* 11.5 */
0.006942840107209529865664152, /* 12.0 */
0.006665247032707682442354394, /* 12.5 */
0.006408994188004207068439631, /* 13.0 */
0.006171712263039457647532867, /* 13.5 */
0.005951370112758847735624416, /* 14.0 */
0.005746216513010115682023589, /* 14.5 */
0.005554733551962801371038690  /* 15.0 */
};

double stirlerr(n)
double n;
{ double nn;

  if (n<15.0)
  { nn = 2.0*n;
    if (nn==(int)nn) return(sferr_halves[(int)nn]);
    return(lgamma(n+1.0) - (n+0.5)*log((double)n)+n - HF_LG_PIx2);
  }

  nn = (double)n;
  nn = nn*nn;
  if (n>500) return((S0-S1/nn)/n);
  if (n>80) return((S0-(S1-S2/nn)/nn)/n);
  if (n>35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
  return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

double bd0(x,np)
double x, np;
{ double ej, s, s1, v;
  int j;
  if (fabs(x-np)<0.1*(x+np))
  {
    s = (x-np)*(x-np)/(x+np);
    v = (x-np)/(x+np);
    ej = 2*x*v; v = v*v;
    for (j=1; ;++j)
    { ej *= v;
      s1 = s+ej/((j<<1)+1);
      if (s1==s) return(s1);
      s = s1;
    }
  }
  return(x*log(x/np)+np-x);
}

/*
  Raw binomial probability calculation.
  (1) This has both p and q arguments, when one may be represented
      more accurately than the other (in particular, in df()).
  (2) This should NOT check that inputs x and n are integers. This
      should be done in the calling function, where necessary.
  (3) Does not check for 0<=p<=1 and 0<=q<=1 or NaN's. Do this in
      the calling function.
*/
double dbinom_raw(x,n,p,q,give_log)
double x, n, p, q;
int give_log;
{ double f, lc;

  if (p==0.0) return((x==0) ? D_1 : D_0);
  if (q==0.0) return((x==n) ? D_1 : D_0);

  if (x==0)
  { lc = (p<0.1) ? -bd0(n,n*q) - n*p : n*log(q);
    return( DEXP(lc) );
  }

  if (x==n)
  { lc = (q<0.1) ? -bd0(n,n*p) - n*q : n*log(p);
    return( DEXP(lc) );
  }

  if ((x<0) | (x>n)) return( D_0 );

  lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x)
         - bd0(x,n*p) - bd0(n-x,n*q);
  f = (PIx2*x*(n-x))/n;

  return( FEXP(f,lc) );
}

double dbinom(x,n,p,give_log)
int x, n;
double p;
int give_log;
{ 
  if ((p<0) | (p>1) | (n<0)) return(INVALID_PARAMS);
  if (x<0) return( D_0 );
  
  return( dbinom_raw((double)x,(double)n,p,1-p,give_log) );
}

/*
  Poisson probability  lb^x exp(-lb) / x!.
  I don't check that x is an integer, since other functions
  that call dpois_raw() (i.e. dgamma) may use a fractional
  x argument.
*/
double dpois_raw(x,lambda,give_log)
int give_log;
double x, lambda;
{
  if (lambda==0) return( (x==0) ? D_1 : D_0 );
  if (x==0) return( DEXP(-lambda) );
  if (x<0) return( D_0 );

  return(FEXP( PIx2*x, -stirlerr(x)-bd0(x,lambda) ));
}

double dpois(x,lambda,give_log)
int x, give_log;
double lambda;
{
  if (lambda<0) return(INVALID_PARAMS);
  if (x<0) return( D_0 );

  return( dpois_raw((double)x,lambda,give_log) );
}

double dbeta(x,a,b,give_log)
double x, a, b;
int give_log;
{ double f, p;

  if ((a<=0) | (b<=0)) return(INVALID_PARAMS);
  if ((x<=0) | (x>=1)) return(D_0);

  if (a<1)
  { if (b<1)                                    /* a<1, b<1 */
    { f = a*b/((a+b)*x*(1-x));
      p = dbinom_raw(a,a+b,x,1-x,give_log);
    }
    else                                        /* a<1, b>=1 */
    { f = a/x;
      p = dbinom_raw(a,a+b-1,x,1-x,give_log);
    }
  }
  else
  { if (b<1)                                    /* a>=1, b<1 */
    { f = b/(1-x);
      p = dbinom_raw(a-1,a+b-1,x,1-x,give_log);
    }
    else                                        /* a>=1, b>=1 */
    { f = a+b-1;
      p = dbinom_raw(a-1,(a-1)+(b-1),x,1-x,give_log);
    }
  }

  return( (give_log) ? p + log(f) : p*f );
}

/*
 *   To evaluate the F density, write it as a Binomial probability
 *   with p = x*m/(n+x*m). For m>=2, use the simplest conversion.
 *   For m<2, (m-2)/2<0 so the conversion will not work, and we must use
 *   a second conversion. Note the division by p; this seems unavoidable
 *   for m < 2, since the F density has a singularity as x (or p) -> 0.
 */
double df(x,m,n,give_log)
double x, m, n;
int give_log;
{ double p, q, f, dens;

  if ((m<=0) | (n<=0)) return(INVALID_PARAMS);
  if (x <= 0.0) return(D_0);

  f = 1.0/(n+x*m);
  q = n*f;
  p = x*m*f;

  if (m>=2)
  { f = m*q/2;
    dens = dbinom_raw((m-2)/2.0, (m+n-2)/2.0, p, q, give_log);
  }
  else
  { f = m*m*q / (2*p*(m+n));
    dens = dbinom_raw(m/2.0, (m+n)/2.0, p, q, give_log);
  }

  return((give_log) ? log(f)+dens : f*dens);
}

/*
 * Gamma density,
 *                  lb^r x^{r-1} exp(-lb*x)
 *      p(x;r,lb) = -----------------------
 *                          (r-1)!
 *
 * If USE_SCALE is defined below, the lb argument will be interpreted
 * as a scale parameter (i.e. replace lb by 1/lb above). Otherwise,
 * it is interpreted as a rate parameter, as above.
 */

/* #define USE_SCALE */

double dgamma(x,r,lambda,give_log)
int give_log;
double x, r, lambda;
{ double pr;

  if ((r<=0) | (lambda<0)) return(INVALID_PARAMS);
  if (x<=0.0) return( D_0 );

#ifdef USE_SCALE
  lambda = 1.0/lambda;
#endif

  if (r<1)
  { pr = dpois_raw(r,lambda*x,give_log);
    return( (give_log) ?  pr + log(r/x) : pr*r/x );
  }

  pr = dpois_raw(r-1.0,lambda*x,give_log);
  return( (give_log) ? pr + log(lambda) : lambda*pr);
}

double dchisq(x, df, give_log)
double x, df;
int give_log;
{
 return(dgamma(x, df/2.0,
  0.5
  ,give_log));
/*
#ifdef USE_SCALE
  2.0
#else
  0.5
#endif
  ,give_log));
*/
}

/*
 * Given a sequence of r successes and b failures, we sample n (\le b+r)
 * items without replacement. The hypergeometric probability is the
 * probability of x successes:
 *
 *                dbinom(x,r,p) * dbinom(n-x,b,p)
 *   p(x;r,b,n) = ---------------------------------
 *                          dbinom(n,r+b,p)
 *
 * for any p. For numerical stability, we take p=n/(r+b); with this choice,
 * the denominator is not exponentially small.
 */
double dhyper(x,r,b,n,give_log)
int x, r, b, n, give_log;
{ double p, q, p1, p2, p3;

  if ((r<0) | (b<0) | (n<0) | (n>r+b))
    return( INVALID_PARAMS );

  if (x<0) return(D_0);

  if (n==0) return((x==0) ? D_1 : D_0);

  p = ((double)n)/((double)(r+b));
  q = ((double)(r+b-n))/((double)(r+b));

  p1 = dbinom_raw((double)x,(double)r,p,q,give_log);
  p2 = dbinom_raw((double)(n-x),(double)b,p,q,give_log);
  p3 = dbinom_raw((double)n,(double)(r+b),p,q,give_log);

  return( (give_log) ? p1 + p2 - p3 : p1*p2/p3 );
}

/*
  probability of x failures before the nth success.
*/
double dnbinom(x,n,p,give_log)
double n, p;
int x, give_log;
{ double prob, f;

  if ((p<0) | (p>1) | (n<=0)) return(INVALID_PARAMS);

  if (x<0) return( D_0 );

  prob = dbinom_raw(n,x+n,p,1-p,give_log);
  f = n/(n+x);

  return((give_log) ? log(f) + prob : f*prob);
}

double dt(x, df, give_log)
double x, df;
int give_log;
{ double t, u, f;

  if (df<=0.0) return(INVALID_PARAMS);

  /*
     exp(t) = Gamma((df+1)/2) /{ sqrt(df/2) * Gamma(df/2) }
            = sqrt(df/2) / ((df+1)/2) * Gamma((df+3)/2) / Gamma((df+2)/2).
     This form leads to a computation that should be stable for all
     values of df, including df -> 0 and df -> infinity.
  */
  t = -bd0(df/2.0,(df+1)/2.0) + stirlerr((df+1)/2.0) - stirlerr(df/2.0);

  if (x*x>df)
    u = log( 1+ x*x/df ) * df/2;
  else
    u = -bd0(df/2.0,(df+x*x)/2.0) + x*x/2.0;

  f = PIx2*(1+x*x/df);

  return( FEXP(f,t-u) );
}
