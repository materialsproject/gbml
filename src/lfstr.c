/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *  Functions for converting string arguments to Locfit's numeric values.
 *  Typically, these will be assigned to appopriate place on one of locfit's structures:
 *    fam(sp) = lffamily(z)
 *    ker(sp) = lfkernel(z)
 *    kt(sp)  = lfketype(z)
 *    link(sp)= lflink(z)
 *    de_itype= deitype(z)
 *    ev(evs) = lfevstr(z)
 *    acri(sp)= lfacri(z)
 *  sp is a pointer to the smpar structure, &lf->sp.
 *  evs is a pointer to the evaluation structure, &lf->evs.
 *  int ppwhat(str) interprets the preplot what argument.
 *  int restyp(str) interprets the residual type argument.
 *
 *  return values of -1 indicate failure/unknown string.
 */

#include "local.h"

int ct_match(z1, z2)
char *z1, *z2;
{ int ct = 0;
  while (z1[ct]==z2[ct])
  { if (z1[ct]=='\0') return(ct+1);
    ct++;
  }
  return(ct);
}

int pmatch(z, strings, vals, n, def)
char *z, **strings;
int *vals, n, def;
{ int i, ct, best, best_ct;
  best = -1;
  best_ct = 0;

  for (i=0; i<n; i++)
  { ct = ct_match(z,strings[i]);
    if (ct==strlen(z)+1) return(vals[i]);
    if (ct>best_ct) { best = i; best_ct = ct; }
  }
  if (best==-1) return(def);
  return(vals[best]);
}

static char *famil[17] =
  { "density", "ate",   "hazard",    "gaussian", "binomial",
    "poisson", "gamma", "geometric", "circular", "obust", "huber",
    "weibull", "cauchy","probab",    "logistic", "nbinomial", "vonmises" };
static int   fvals[17] = 
  { TDEN,  TRAT,  THAZ,  TGAUS, TLOGT,
    TPOIS, TGAMM, TGEOM, TCIRC, TROBT, TROBT,
    TWEIB, TCAUC, TPROB, TLOGT, TGEOM, TCIRC };
int lffamily(z)
char *z;
{ int quasi, robu, f;
  quasi = robu = 0;
  while ((z[0]=='q') | (z[0]=='r'))
  { quasi |= (z[0]=='q');
    robu  |= (z[0]=='r');
    z++;
  }
  f = pmatch(z,famil,fvals,16,-1);
  if ((z[0]=='o') | (z[0]=='a')) robu = 0;
  if (f==-1)
  { WARN(("unknown family %s",z));
    f = TGAUS;
  }
  if (quasi) f += 64;
  if (robu)  f += 128;
  return(f);
}

static char *wfuns[13] = {
  "rectangular", "epanechnikov", "bisquare",    "tricube",
  "triweight",   "gaussian",     "triangular",  "ququ",
  "6cub",        "minimax",      "exponential", "maclean", "parametric" };
static int wvals[13] = { WRECT, WEPAN, WBISQ, WTCUB,
  WTRWT, WGAUS, WTRIA, WQUQU, W6CUB, WMINM, WEXPL, WMACL, WPARM };
int lfkernel(char *z)
{ return(pmatch(z, wfuns, wvals, 13, WTCUB));
}

static char *ktype[5] = { "spherical", "product", "center", "lm", "zeon" };
static int   kvals[5] = { KSPH, KPROD, KCE, KLM, KZEON };
int lfketype(char *z)
{ return(pmatch(z, ktype, kvals, 5, KSPH));
}

static char *ltype[8] = { "default", "canonical", "identity", "log",
                          "logi",    "inverse",   "sqrt",     "arcsin" };
static int   lvals[8] = { LDEFAU, LCANON, LIDENT, LLOG,
                          LLOGIT, LINVER, LSQRT,  LASIN };
int lflink(char *z)
{ return(pmatch(z, ltype, lvals, 8, LDEFAU));
}

static char *etype[11]= { "tree",     "phull", "data", "grid", "kdtree",
                          "kdcenter", "cross", "preset", "xbar", "none",
                          "sphere" };
static int   evals[11]= { ETREE, EPHULL, EDATA, EGRID, EKDTR,
                          EKDCE, ECROS,  EPRES, EXBAR, ENONE, ESPHR };
int lfevstr(char *z)
{ return(pmatch(z, etype, evals, 11, ETREE));
}

static char *itype[7] = { "default", "multi", "product", "mlinear",
                          "hazard",  "sphere", "monte" };
static int   ivals[7] = { IDEFA, IMULT, IPROD, IMLIN, IHAZD, ISPHR, IMONT };
int deitype(char *z)
{ return(pmatch(z, itype, ivals, 6, IDEFA));
}

static char *atype[5] = { "none", "cp", "ici", "mindex", "ok" };
static int   avals[5] = { ANONE, ACP, AKAT, AMDI, AOK };
int lfacri(char *z)
{ return(pmatch(z, atype, avals, 5, ANONE));
}

static char *rtype[8] = { "deviance", "d2",    "pearson", "raw",
                          "ldot",     "lddot", "fit",     "mean" };
static int   rvals[8] = { RDEV, RDEV2, RPEAR, RRAW, RLDOT, RLDDT, RFIT, RMEAN};

static char *whtyp[8] = { "coef", "nlx", "infl", "band",
                          "degr", "like", "rdf", "vari" };
static int   whval[8] = { PCOEF, PNLX, PT0, PBAND, PDEGR, PLIK, PRDF, PVARI };

int restyp(z)
char *z;
{ int val;
  
  val = pmatch(z, rtype, rvals, 8, -1);
  if (val==-1) ERROR(("Unknown type = %s",z));
  return(val);
}

int ppwhat(z)
char *z;
{ int val;
  
  val = pmatch(z, whtyp, whval, 8, -1);
  if (val==-1) ERROR(("Unknown what = %s",z));
  return(val);
}
