/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 * 
 *
 *  Most of the changes formerly needed here are handled through
 *  the Makefiles and #ifdef's.
 */

#ifndef I_LF_H
#define I_LF_H
#include "notR.h"
//#include <R.h>  // RJN

/*
 *   DIRSEP: '/' for unix; '\\' for DOS
 */
#ifdef DOS
#define DIRSEP '\\'
#else
#define DIRSEP '/'
#endif

/*
   Some older math libraries have no lgamma() function, and gamma(arg)
   actually returns log(gamma(arg)). If so, you need to change
   LGAMMA macro below.

   If all else fails, you can also use lflgamma().

   Use the definitions for erf, erfc and daws only if your
   math libraries don't include these functions.
 */
#ifdef DOS
#define LGAMMA(arg) lflgamma(arg)
#define erf(x) lferf(x)
#define erfc(x) lferfc(x)
#else
#define LGAMMA(arg) lgamma(arg)
#endif
#define daws(x) lfdaws(x)

/******** NOTHING BELOW HERE NEEDS CHANGING **********/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define RVERSION

#ifdef SWINVERSION
#define SVERSION
#include "newredef.h"
#endif

#ifdef RVERSION

/* #typedef int Sint is defined in R.h */
typedef int Sint;          // RJN
//#include <R.h>           // RJN
//#include <Rdefines.h>    // RJN
//#include <Rinternals.h>  // RJN
#define list_elt(ev,i) VECTOR_PTR(ev)[i]
#define dval2(ev,i,j) NUMERIC_POINTER(list_elt(ev,i))[j]
#define dvec2(ev,i)   NUMERIC_POINTER(list_elt(ev,i))
#define ivec2(ev,i)   INTEGER_POINTER(list_elt(ev,i))
#undef pmatch
//#define printf Rprintf   // RJN
//#define printe REprintf  // RJN

#else

#ifdef SVERSION
#include <S.h>
typedef long int Sint;
typedef s_object * SEXP;
#define list_elt(ev,i) LIST_POINTER(ev)[i]
#define dval2(ev,i,j) NUMERIC_POINTER(list_elt(ev,i))[j]
#define dvec2(ev,i)   NUMERIC_POINTER(list_elt(ev,i))
#define ivec2(ev,i)   INTEGER_POINTER(list_elt(ev,i))
#else
typedef int Sint;
#endif

#endif

#ifdef RVERSION
//#undef LGAMMA                         // RJN
//#define LGAMMA(arg) Rf_lgammafn(arg)  // RJN
//extern double Rf_lgammafn();          // RJN
#define SVERSION
#endif

#include "mutil.h"
#include "tube.h"

#include "lfcons.h"

typedef char varname[15];

#ifdef CVERSION
#include "cversion.h"
#endif

#include "lfstruc.h"
#include "design.h"
#include "lffuns.h"

#ifdef CVERSION
#undef printf
#define printf lfprintf
extern int lfprintf(const char *format, ...);
extern int printe(const char *format, ...);
/* #else
   #define printe printf */
#endif

#ifdef ERROR
#undef ERROR
#endif

#ifdef WARN
#undef WARN
#endif

/* #define ERROR(args) {printe("Error: "); printe args; printe("\n"); lf_error=
1;} */
//#define ERROR(args) {error args; lf_error=1;}  // RJN
/* #define WARN(args)  {printe("Warning: "); printe args; printe("\n"); } */
//#define WARN(args)  warning args;  // RJN

#define ERROR(args) {printf("Error: ");   printf args; printf("\n"); lf_error=1;}  // RJN
#define WARN(args)  {printf("Warning: "); printf args; printf("\n");}              // RJN

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SGN(x) (((x)>0) ? 1 : -1)
#define SQR(x) ((x)*(x))
#define NOSLN 0.1278433
#define GFACT 2.5
#define EFACT 3.0

#define MAXCOLOR 20
#define MAXWIN 5

#define ISWAP(a,b) { int zz; zz = a; a = b; b = zz; }
extern int lf_error;

#endif /* I_LF_H */
