/*
 *   Copyright (c) 1998-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Structures and function definitions for the C version interface.
 */

/* typedef char varname[15]; */

/*
 *  Define the vari type for locfit variables and related macros.
 */

typedef struct {
  varname name;
  int n, bytes, mode, stat;
  double *dpr; } vari;

#define checkvarlen(v,n,name,mode) (createvar(name,STSYSTEM,n,mode))
#define vmode(v) ((v)->mode)
#define vlength(v) ((v)->n)

typedef struct {
  char *arg, *val;
  vari *result;
  int used; } carg;

typedef struct {
  void (*AddColor)(), (*SetColor)(), (*ClearScreen)(), (*TextDim)(), (*DoText)();
  void (*DrawPoint)(), (*DrawLine)(), (*DrawPatch)(), (*wrapup)();
  int (*makewin)(), ticklength, defth, deftw;
} device;

typedef struct {
  vari *data[MXDIM], *fit, *se;
  int d, wh, gr;
} pplot;

typedef struct {
  char cmd;
  double x, *v, (*f)();
  int m, nx[3];
  vari *vv; } arstruct;

typedef struct {
  vari *x, *y, *z;
  char type;
  int id, t, n, nx, ny, pch; } plxyz;

typedef struct {
  double theta, phi, xl[2], yl[2], zl[2], sl[10];
  int id, ty, nsl;
  char main[50], xlab[50], ylab[50], zlab[50];
  vari *track, *xyzs; } plots;

#define PLNONE 0
#define PLDATA 1
#define PLFIT  2
#define PLTRK  4

struct lfcol {
  char name[10];
  int n, r, g, b;
};


/* FILES IN THE src-c DIRECTORY */

/* arith.c */
extern int arvect(), intitem();
extern double areval(), arith(), darith(), dareval();
extern vari *varith(), *saveresult(), *arbuild();

/* c_args.c */
#define argused(v,i) (((carg *)viptr(v,i))->used)
#define setused(v,i) { ((carg *)viptr(v,i))->used = 1; }
#define setunused(v,i) { ((carg *)viptr(v,i))->used = 0; }
#define argarg(v,i) (((carg *)viptr(v,i))->arg)
#define argvalis(v,i,z) (strcmp(argval(v,i),z)==0)
extern char *argval(), *getargval();
extern int getarg(), readilist(), getlogic();

/* cmd.c */
extern int dispatch();
extern void setuplf(), recondat(), cmdint();
extern double backtr(), docrit();

/* c_lf.c */
extern vari *vfitted();
extern void cfitted(), cwdiag();

/* c_plot.c */
extern void plotdata(), plotfit(), plottrack(), plotopt(), setplot();

/* help.c */
extern void example();

/* lfd.c */
extern void doreaddata(), dosavedata(), dosavefit();
extern int  setfilename();

/* main.c */
extern void SetWinDev();

/* makecmd.c */
extern vari *getcmd();
extern void makecmd(), del_lines(), inc_forvar(), dec_forvar();

/* post.c */
extern void SetPSDev();

/* pout.c */
extern int pretty();
extern void displayplot();
extern void plotmaple(), plotmathe(), plotmatlb(), plotgnup(), plotxwin();

/* random.c */
extern double rnorm(), rexp(), runif(), rpois();
extern void rseed();

/* readfile.c */
extern void readfile();

/* scbmax.c */
extern void cscbmax();

/* vari.c */
extern int vbytes();
extern vari *createvar(), *findvar(), *growvar();
extern void initdb(), deletevar(), deletename(), deleteifhidden(), setvarname();
extern void *viptr(), vassn();
extern double *vdptr(), vitem();
