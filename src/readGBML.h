int readGBMLheader(FILE *filePtr, double *learningRate, double *responseMean,
      Sint *nInteractions, Sint *nDescriptors, Sint *nSteps);

int readGBMLdata(FILE *filePtr, Sint nInteractions, Sint **descriptorIndexes, double **xev,
      double **coef, double **sv, Sint **cell, double **pc, double **scale,
      Sint **nvc, Sint **mi, double **dp, Sint **mg, Sint **sty);
