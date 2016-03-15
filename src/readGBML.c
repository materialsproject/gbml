//
//   Copyright 2015, The Materials Project
//

// Routines for reading a GBM-Locfit model from a file.

#include "local.h"
#include "readGBML.h"

int readGBMLheader(FILE *filePtr, double *learningRate, double *responseMean,
      Sint *nInteractions, Sint *nDescriptors, Sint *nSteps)
{
  if( (fread(learningRate,  sizeof(double),    1, filePtr) != 1) ||
      (fread(responseMean,  sizeof(double),    1, filePtr) != 1) ||
      (fread(nInteractions, sizeof(Sint),      1, filePtr) != 1) ||
      (fread(nDescriptors,  sizeof(Sint),      1, filePtr) != 1) ||
      (fread(nSteps,        sizeof(Sint),      1, filePtr) != 1) ) return(1);
  return(0);
}

int readGBMLdata(FILE *filePtr, Sint nInteractions, Sint **descriptorIndexes, double **xev,
      double **coef, double **sv, Sint **cell, double **pc, double **scale,
      Sint **nvc, Sint **mi, double **dp, Sint **mg, Sint **sty)
{
  Sint nItems, rItems;

  if((*descriptorIndexes = (Sint *) calloc(nInteractions, sizeof(Sint))) == NULL) return(2);
  if((rItems = fread(*descriptorIndexes, sizeof(Sint), nInteractions, filePtr)) != nInteractions) return(3);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(4);
  if((*xev = (double *) calloc(nItems, sizeof(double))) == NULL) return(5);
  if((rItems = fread(*xev, sizeof(double), nItems, filePtr)) != nItems) return(10*rItems+6);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(7);
  if((*coef = (double *) calloc(nItems, sizeof(double))) == NULL) return(8);
  if((rItems = fread(*coef, sizeof(double), nItems, filePtr)) != nItems) return(9);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(10);
  if((*sv = (double *) calloc(nItems, sizeof(double))) == NULL) return(11);
  if((rItems = fread(*sv, sizeof(double), nItems, filePtr)) != nItems) return(12);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(13);
  if((*cell = (Sint *) calloc(nItems, sizeof(Sint))) == NULL) return(14);
  if((rItems = fread(*cell, sizeof(Sint), nItems, filePtr)) != nItems) return(15);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(16);
  if((*pc = (double *) calloc(nItems, sizeof(double))) == NULL) return(17);
  if((rItems = fread(*pc, sizeof(double), nItems, filePtr)) != nItems) return(18);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(19);
  if((*scale = (double *) calloc(nItems, sizeof(double))) == NULL) return(20);
  if((rItems = fread(*scale, sizeof(double), nItems, filePtr)) != nItems) return(21);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(22);
  if((*nvc = (Sint *) calloc(nItems, sizeof(Sint))) == NULL) return(23);
  if((rItems = fread(*nvc, sizeof(Sint), nItems, filePtr)) != nItems) return(24);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(25);
  if((*mi = (Sint *) calloc(nItems, sizeof(Sint))) == NULL) return(26);
  if((rItems = fread(*mi, sizeof(Sint), nItems, filePtr)) != nItems) return(27);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(28);
  if((*dp = (double *) calloc(nItems, sizeof(double))) == NULL) return(29);
  if((rItems = fread(*dp, sizeof(double), nItems, filePtr)) != nItems) return(30);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(31);
  if((*mg = (Sint *) calloc(nItems, sizeof(Sint))) == NULL) return(32);
  if((rItems = fread(*mg, sizeof(Sint), nItems, filePtr)) != nItems) return(33);

  if(fread(&nItems, sizeof(Sint), 1, filePtr) != 1) return(34);
  if((*sty = (Sint *) calloc(nItems, sizeof(Sint))) == NULL) return(35);
  if((rItems = fread(*sty, sizeof(Sint), nItems, filePtr)) != nItems) return(36);

  return(0);
}
