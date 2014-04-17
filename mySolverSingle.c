//single neuron
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nr.h"
#include "nrutil.h"
#include "globalVars.h"
#include "NetworkEquations.h"
//void (*derive)(double, double *, double *);
// compile as - gcc nrutil.c rk4.c rkdumb.c -g mysolver.c -lm -o mysolver
// -g for gdb
//-lm reqd. for linking math.h
// 
// GLOBAL VARS
extern double **y, *xx, *input_cur;
double dt; //integration time step
void main(int argc, char **argv) {
    int dim = 4;
    double *vstart;
    double x1 = 0;
    double x2 = 100 ;
    int nSteps;
    FILE *fp;
    int loopIdx=0;
    //    char *saveFormat;

    dt = 0.025;
    nSteps = (int)((x2 - x1) / dt);
    xx = vector(1, nSteps);
    y = matrix(1, dim, 1, nSteps);
    input_cur = vector(1, nSteps);
    vstart = vector(1, dim);
   
      for(loopIdx = 1; loopIdx < nSteps+1;  ++loopIdx) {
	if(loopIdx * dt > 5 && loopIdx *dt <= 50  )
	    input_cur[loopIdx] = 10;
	  else 
	    input_cur[loopIdx] = 0;
      }
    vstart[1] = 0;
    vstart[2] = 0.3176;
    vstart[3] = 0.1;
    vstart[4] = 0.5961;

    
    rkdumb(vstart, dim, x1, x2, nSteps, derivs);
    fp = fopen("outputFile.csv", "w");
    //    for(loopIdx = 1; loopIdx < dim + 1 ; ++loopIdx( {
    for(loopIdx = 1; loopIdx < nSteps+1; ++loopIdx) 
      {
	//[t, V, n, z, h, I_input]
	fprintf(fp, "%f %f %f %f %f %f\n", xx[loopIdx], y[1][loopIdx], y[2][loopIdx], y[3][loopIdx], y[4][loopIdx], input_cur[loopIdx]);
	if(loopIdx%100 == 0) {
	  printf("\r%d", loopIdx);
	}
      }

    fclose(fp);
    //  free_vector(xx, 1, nSteps);
    //  free_vector(input_cur, 1, nSteps);
    //  free_matrix(y, 1, dim, 1, nSteps);
  }   
