//
//   Copyright 2015, The Materials Project
//

#include <stdio.h>
#include "local.h"
#include "predict.h"
#include "readGBML.h"
#include "spreplot.h"

extern int deitype(char *);  /* in lfstr.c */
static lfit lf;

int predict(char *filename, int nPredictions, double *descriptors, double *predictions)
{
  Sint iFlag, iPrediction, iStep;
  FILE *filePtr;

  double *newData, *incrPredictions, *stdErrors=NULL;

  Sint nInteractions, nDescriptors, nSteps;
  double learningRate, responseMean;

  Sint *descriptorIndexes, *cell, *nvc, *mi, *mg, *sty;
  double *xev, *coef, *sv, *pc, *scale, *dp;

  char what0[5]="coef";
  char what1[5]="none";
  char *what[2]; what[0]=what0; what[1]=what1;
  Sint *dv=NULL;
  Sint nd=0;
  Sint where=1;
  void *bs = NULL;

  // Open input file for reading
  if(!(filePtr = fopen(filename,"rb"))) {
    fprintf(stderr, "Problem opening input file: %s\n", filename);
    return(1);
  }

  // Read header (28 bytes)
  if(readGBMLheader(filePtr, &learningRate, &responseMean, &nInteractions, &nDescriptors, &nSteps)) {
    fprintf(stderr, "Problem reading header from input file: %s\n", filename);
    return(2);
  }

  // Allocate storage
  int nNewData = nPredictions * nInteractions;
  if((newData         = (double *) calloc(nNewData,     sizeof(double))) == NULL) return(3);
  if((incrPredictions = (double *) calloc(nPredictions, sizeof(double))) == NULL) return(4);

  // Initialize predictions to responseMean
  for(iPrediction=0; iPrediction < nPredictions; ++iPrediction)
    predictions[iPrediction] = responseMean;

  // Loop over steps
  for(iStep=0; iStep < nSteps; ++iStep) {

    // Read GBML data for iStep
    if((iFlag = readGBMLdata(filePtr, nInteractions, &descriptorIndexes, &xev, &coef, &sv, &cell, &pc,
        &scale, &nvc, &mi, &dp, &mg, &sty))) {
      fprintf(stderr, "Error %d while reading data for step %d from input file.\n", iFlag, iStep);
      return(5);
    }

    // Fill newData with appropriate descriptors
    int iInteraction, newIndex, oldIndex;
    for(iInteraction=0; iInteraction < nInteractions; ++iInteraction) {
      for(iPrediction=0; iPrediction < nPredictions; ++iPrediction) {
        newIndex = nPredictions * iInteraction + iPrediction;
        oldIndex = nPredictions * (descriptorIndexes[iInteraction] - 1) + iPrediction;
        newData[newIndex] = descriptors[oldIndex];
      }
    }

    // Call spreplot to calculate incremental prediction
    spreplot(xev, coef, sv, cell, newData, incrPredictions, stdErrors, pc, scale, &nPredictions,
        nvc, mi, dp, mg, dv, &nd, sty, &where, what, &bs);

    // Add incremental predictions
    for(iPrediction=0; iPrediction < nPredictions; ++iPrediction)
      predictions[iPrediction] += (learningRate * incrPredictions[iPrediction]);

    // Free memory
    free(descriptorIndexes); free(nvc); free(mi);
    free(xev); free(coef); free(sv); free(cell);
    free(pc); free(scale); free(dp); free(mg);  free(sty);
  }

/*
  printf("\nPredictions:");
  for(iPrediction=0; iPrediction < nPredictions; ++iPrediction)
    printf("  %f", predictions[iPrediction]);
  printf("\n");
*/

  // Close input file
  fclose(filePtr);
  return(0);
}
