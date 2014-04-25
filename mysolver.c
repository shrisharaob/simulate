#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "nr.h"
#include "nrutil.h"
#include "globalVars.h"
#include "varProtos.h"
#include "auxFuncProtos.h"
#include "optmNwEq.c"


//#include "nwEq2.h"
//void (*derive)(double, double *, double *);
// compile as - gcc nrutil.c rk4.c auxFunctions.c rkdumb.c -g mysolver.c -lm -o mysolver
// gcc ran1.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb.c -g mysolver.c -lm -o mysolver
// -g for gdb
//-lm reqd. for linking math.h
// 
// GLOBAL VARS
extern double **y, *xx, *input_cur, *IF_SPK, *expSum, *iSynap, 
  *gaussNoiseE, *gaussNoise, theta, contrast, *gFF, *iFF, *rTotal, muE, muI;
extern FILE *spkTimesF, *outVars;
double dt; //integration time step
void main(int argc, char **argv) {
    // declare
    int dim = 4;
    double *vstart;
    double x1 = 0;
    double x2 = 100;
    int nSteps;
    FILE *fp;
    int loopIdx=0;
    int kNeuron, clmNo;
    double *spkTimes;
    int nSpks = 0;
    // initialize
    spkTimes = vector(1, nSteps);
    dt = DT;
    nSteps = (int)((x2 - x1) / dt);
    xx = vector(1, nSteps);
    y = matrix(1, N_StateVars * N_Neurons, 1, nSteps);
    input_cur = vector(1, nSteps);
    vstart = vector(1, N_StateVars * N_Neurons);
    IF_SPK = vector(1, N_Neurons);   
    spkTimesFp = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/spkTimes","w");
    expSum = vector(1, N_Neurons); // synap input
    iSynap = vector(1, N_Neurons); 
    gEE = vector(1, NE);   
    gEI = vector(1, NE);   
    gIE = vector(1, NI);   
    gII = vector(1, NI);
    conMat = matrix(1, N_Neurons, 1, N_Neurons);
    gaussNoiseE = vector(1, NE);
    gaussNoiseI = vector(1, NI);
    iBg = vector(1, N_Neurons);
    gFF = vector(1, N_Neurons);
    iFF = vector(1, N_Neurons);
    rTotal = vector(1, N_Neurons);
    rTotalPrev = vector(1, N_Neurons);
    tempRandnPrev = vector(1, N_Neurons);
    tempRandnNew = vector(1, N_Neurons);
    tempCurE = vector(1, N_Neurons);
    tempCurI = vector(1, N_Neurons);
    Itgrl = vector(1, N_Neurons);
    ItgrlOld = vector(1, N_Neurons);
    randnXiA = vector(1, N_Neurons);
    randwZiA = matrix(1, N_Neurons, 1, 4);
    randuDelta = vector(1, N_Neurons);
    randuPhi = matrix(1, N_Neurons, 1, 3);
    outVars = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/outvars", "w");
    isynapFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/isynapEI", "w");
    rTotalFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/rTotal", "w");
    gbgrndFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/gBg", "w");
    srand(time(NULL));
    //    genConMat();
    AuxRffTotal();
    //    GenConMat02();
    /* /\********\/ */
    //conMat[1][1] = 0; 
         conMat[1][2] = 1; 
    /* conMat[1][3] = 1; */
    //    conMat[2][1] = 1; 
    /* conMat[2][2] = 0; */
    /* conMat[2][3] = 1; */
    /* conMat[3][1] = 1; */
    /* conMat[3][2] = 0; */
    /* conMat[3][3] = 0; */
    // compute
    printf("\nNE = %d\n", NE);
    printf("\nNI = %d\n", NI);
    printf("\nK = %d\n", (int)K);
    printf("computing...\n");
    //Parse inputs
    if(argc > 1) {
      theta = atof(argv[2]);
      contrast = atof(argv[3]);
      muE = atof(argv[4]);
      muI = atof(argv[5]);
      for(loopIdx = 1; loopIdx < nSteps+1;  ++loopIdx) {
	if(loopIdx * dt > atof(argv[6]) && loopIdx *dt <= atof(argv[7])) {
	  input_cur[loopIdx] = atof(argv[8]); 
	}
	else { 
	  input_cur[loopIdx] = 0;
	}
      }
    }
    else {
      theta = 0;
      contrast = 0.25;
      muE = 0.1;
      muI = 0.1;
      for(loopIdx =1; loopIdx < nSteps+1;  ++loopIdx) {
	input_cur[loopIdx] = 10; 
      }
    }
    // initialize states
    for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
      clmNo =  (kNeuron - 1) * N_StateVars;
      vstart[1 + clmNo] = 0;
      vstart[2 + clmNo] = 0.3176;
      vstart[3 + clmNo] = 0.1;
      vstart[4 + clmNo] = 0.5961;
    }
    rkdumb(vstart, N_StateVars * N_Neurons, x1, x2, nSteps, derivs);
    printf("Done! \n");
    fclose(spkTimesFp);
    fp = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/outputFile.csv", "w");
    for(loopIdx = 1; loopIdx < nSteps+1; ++loopIdx) {
        fprintf(fp, "%f ", xx[loopIdx]);
        for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
          clmNo =  (kNeuron - 1) * N_StateVars;
          //[t, V_m, n, z, h, I_input]
	  //	  fprintf(fp, "%f %f %f %f %f ", y[1 + clmNo][loopIdx], y[2 + clmNo][loopIdx], y[3 + clmNo][loopIdx], y[4 + clmNo][loopIdx], input_cur[loopIdx]);
	  fprintf(fp, "%f ", y[1 + clmNo][loopIdx]);
      }
       fprintf(fp, "\n");
       if(loopIdx%100 == 0) {
         printf("\r%d", loopIdx);
       }
      }
    printf("\n");
    printf("nSteps = %d \n", loopIdx);
    fclose(fp);
    fclose(outVars);
    fclose(isynapFP);
    fclose(rTotalFP);
    fclose(gbgrndFP);
    free_vector(expSum, 1, N_Neurons);
    free_vector(iSynap, 1, N_Neurons); 
    free_vector(gEE, 1, NE);   
    free_vector(gIE, 1, NI);   
    free_vector(gEI, 1, NE);   
    free_vector(gII, 1, NI);   
    free_vector( gaussNoiseE, 1, NE);
    free_vector(gaussNoiseI, 1, NI);
    free_vector(rTotal, 1, N_Neurons);
    free_vector(rTotalPrev, 1, N_Neurons);
    free_vector(tempRandnPrev, 1, N_Neurons);
    free_vector(tempRandnNew, 1, N_Neurons);
    free_vector(iBg, 1, N_Neurons);
    free_vector(gFF, 1, N_Neurons);
    free_vector(iFF, 1, N_Neurons);
    free_vector(tempCurI, 1, N_Neurons);
    free_vector(tempCurE, 1, N_Neurons);
    free_matrix(conMat, 1, N_Neurons, 1, N_Neurons);
    free_matrix(randwZiA, 1, N_Neurons, 1, 4);
    free_matrix(randuPhi, 1, N_Neurons, 1, 3);
}
