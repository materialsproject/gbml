/*
 *   Copyright (c) 1998-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Structures, typedefs etc used in Locfit
 */

typedef struct {
  double *wk, *coef, *xbar, *f;
  jacobian xtwx;
  int lwk, haspc;
} paramcomp;
#define haspc(pc) ((pc)->haspc)

typedef struct {
  double *x[MXDIM];
  double *y;
  double *w;
  double *b;
  double *c;
  double sca[MXDIM];
  double xl[2*MXDIM];
  int n, d, ord;
  int sty[MXDIM];
  varname yname, xname[MXDIM], wname, bname, cname;
} lfdata;
#define resp(lfd,i) (((lfd)->y==NULL) ? 0.0 : (lfd)->y[i])
#define base(lfd,i) (((lfd)->b==NULL) ? 0.0 : (lfd)->b[i])
#define prwt(lfd,i) (((lfd)->w==NULL) ? 1.0 : (lfd)->w[i])
#define cens(lfd,i) (((lfd)->c==NULL) ? 0 : (int)(lfd)->c[i])
#define datum(lfd,i,j) ((lfd)->x[i][j])
#define dvari(lfd,i)   ((lfd)->x[i])

typedef struct {
  double nn, fixh, adpen;
  int ker, kt;
  int deg, deg0, p;
  int acri;
  int fam, lin;
  int ubas;
  double (*vb)();
  void (*vbasis)();
} smpar;
#define nn(sp)   ((sp)->nn)
#define fixh(sp) ((sp)->fixh)
#define pen(sp)  ((sp)->adpen)
#define ker(sp)  ((sp)->ker)
#define kt(sp)   ((sp)->kt)
#define deg(sp)  ((sp)->deg)
#define deg0(sp) ((sp)->deg0)
#define npar(sp) ((sp)->p)
#define acri(sp) ((sp)->acri)
#define ubas(sp) ((sp)->ubas)
#define fam(sp)  ((sp)->fam)
#define link(sp) ((sp)->lin)

typedef struct {
  int deriv[MXDEG+2];
  int nd;
} deriv;

typedef struct {
  int ev;
  double *sv;
  double cut;
  double fl[2*MXDIM];
  Sint *iwk, *ce, *s, *lo, *hi;
  int liw, nce, ncm, maxk;
  int mg[MXDIM];
  void (*espec)();
} evstruc;
#define ev(evs)   ((evs)->ev)
#define cut(evs)  ((evs)->cut)
#define mk(evs)   ((evs)->maxk)
#define mg(evs)   ((evs)->mg)

typedef struct {
  double *xev, *coef, *nlx, *t0, *lik, *h, *deg, *L;
  int lev, lwk, ll;
  int d, dcor, geth, hasd;
  int nv, nvm;
  double df0, df1, llk, rv, rsc;
  double kap[10];
} fitpt;
#define evp(fp)     ((fp)->xev)
#define evpt(fp,i)  (&(fp)->xev[(i)*(fp)->d])
#define evptx(fp,i,k) ((fp)->xev[(i)*(fp)->d+(k)])
#define df0(fp) ((fp)->df0)
#define df1(fp) ((fp)->df1)
#define llk(fp) ((fp)->llk)
#define dc(fp)  ((fp)->dcor)
#define geth(fp) ((fp)->geth)
#define rv(fp)  ((fp)->rv)
#define rsc(fp) ((fp)->rsc)

typedef struct {
  int       lf_init_id;
  lfdata    lfd;
  smpar     sp;
  evstruc   evs;
  fitpt     fp;
  deriv     dv;
  paramcomp pc;
  } lfit;
#define LF_INIT_ID 34897239
