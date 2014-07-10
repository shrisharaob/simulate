#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "nr.h"
#include "nrutil.h"
#include "config.h"
#include "globalVars.h"
#include "varProtos.h"
#include "auxFuncProtos.h"
#include "optmNwEq.c"

// compile as - gcc nrutil.c rk4.c auxFunctions.c rkdumb.c -g mysolver.c -lm -o mysolver
// gcc ran1.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb.c -g mysolver.c -lm -o mysolver
// -g for gdb
//-lm reqd. for linking math.h

/* GLOBAL VARS */
extern double **y, *xx, *input_cur, *IF_SPK, *expSum, *iSynap, 
  *gaussNoiseE, *gaussNoise, theta, contrast, *gFF, *iFF, *rTotal, muE, muI;
extern FILE *spkTimesF, *outVars;
sparseMat *sConMat[N_Neurons + 1]; // index staring with 1
void main(int argc, char **argv) {
    // ***** DECLARATION *****//
  int dim = 4;
    double *vstart, *spkTimes;;
    double x1 = 0, // simulation start time
      x2 = 10000.0, // simulation end time
      thetaStep = 0;
    int nSteps, nThetaSteps;
    int kNeuron, clmNo, loopIdx=0;
    long idem;
    FILE *vmFP1;
    clock_t begin, end;
    char fileSuffix[128];
    //***** PARSE INPUT ARGS *****//
    if(argc > 1) {
      theta = atof(argv[1]);
    }
    else {
      theta = 0;
    }
    // ***** INITIALIZATION *****//
    dt = DT;
    nSteps = (int)((x2 - x1) / dt);
    nTotSpks = 0;
    xx = vector(1, STORE_LAST_N_STEPS);
    y = matrix(1, N_Neurons, 1, STORE_LAST_N_STEPS);
    vstart = vector(1, N_StateVars * N_Neurons);
    IF_SPK = vector(1, N_Neurons);   
    expSum = vector(1, N_Neurons); // synap input
    iSynap = vector(1, N_Neurons); 
    gEI_I = vector(1, N_Neurons); // E --> E & E --> I
    gEI_E = vector(1, N_Neurons); // I --> E & I --> I  
    spkTimes = vector(1, nSteps);
    conMat = matrix(1, N_Neurons, 1, N_Neurons);
    //    sConMat = (sparseMat *)malloc((N_Neurons+1) * sizeof(sparseMat));
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
    printf("mem allocation done\n");
    printf("NE = %d\n", NE);
    printf("NI = %d\n", NI);
    printf("K = %d\n", (int)K);
    printf("dt = %f, tStop = %f, nSteps = %d\n", DT, x2, nSteps);  
    printf("computing...\n");
    strcpy(filebase, FILEBASE);
    vmFP = fopen(strcat(filebase, "vm.csv"), "w");
    strcpy(filebase, FILEBASE);
    strcat(filebase, "spkTimes_theta");
    sprintf(fileSuffix, "%03.0f", theta);
    strcat(filebase, fileSuffix);
    spkTimesFp = fopen(strcat(filebase, ".csv"),"w");
    printf("\n%s\n", filebase);
    strcpy(filebase, FILEBASE);
    outVars = fopen(strcat(filebase, "outvars.csv"), "w");
    strcpy(filebase, FILEBASE);
    isynapFP = fopen(strcat(filebase, "isynapEI"), "w");
    strcpy(filebase, FILEBASE);
    rTotalFP = fopen(strcat(filebase, "rTotal.csv"), "w");
    strcpy(filebase, FILEBASE);
    gbgrndFP = fopen(strcat(filebase, "gBg.csv"), "w");
    strcpy(filebase, FILEBASE);
    gEEEIFP = fopen(strcat(filebase, "gEEEI.csv"), "w");
    strcpy(filebase, FILEBASE);
    srand(time(NULL)); // set the seed for random number generator
    fprintf(stdout, "generating conMat..."); fflush(stdout);
    genConMat(); // Generate conection matrix
    printf("done\n");
    //    GenConMat02();
    GenSparseConMat(sConMat);
    //    GenSparseConMatDisp(sConMat);
    AuxRffTotal(); /* auxillary function, generates random variables for the 
                      simulation run; which are used for approximating FF input */
    if(thetaStep > 0) {
      thetaVec = vector(1, 360 / thetaStep);
      LinSpace(0, 360, thetaStep, thetaVec, &nThetaSteps); 
    }
    else {
      thetaVec = vector(1, 1);
      thetaVec[1] = 0;
      nThetaSteps = 1;
    }

    contrast = 0.0;
    muE = 0.1;
    muI = 0.1;
    printf("theta = %f contrast = %f\n", theta, contrast);
    /* INITIALIZE STATE VARIABLES */
    vmFP1 = fopen("vmstart.csv", "w");
    for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
      clmNo =  (kNeuron - 1) * N_StateVars;
      idem = -1 * rand();
      vstart[1 + clmNo] = -70 +  40 * ran1(&idem); // Vm(0) ~ U(-70, -30)
      fprintf(vmFP1, "%f\n", vstart[1+clmNo]);
      vstart[2 + clmNo] = 0.3176;
      vstart[3 + clmNo] = 0.1;
      vstart[4 + clmNo] = 0.5961;
    }
    /* INTEGRATE */
    begin = clock();
    for(loopIdx = 1; loopIdx <= nThetaSteps; ++loopIdx) {
      //      theta = thetaVec[loopIdx];
      rkdumb(vstart, N_StateVars * N_Neurons, x1, x2, nSteps, derivs);
    }
    printf("Done! \n");
    end = clock();
    printf("\n time spent integrating : %.2fs", (double)(end - begin) / CLOCKS_PER_SEC);
    fclose(spkTimesFp);
    /* SAVE TO DISK */
    if(theta == 18.0) {
      for(loopIdx = 1; loopIdx <= STORE_LAST_N_STEPS; ++loopIdx) {
        fprintf(vmFP, "%f ", xx[loopIdx]);
        for(kNeuron = 1; kNeuron <= N_Neurons; ++kNeuron) {
          fprintf(vmFP, "%f ", y[kNeuron][loopIdx]);
        }
        fprintf(vmFP, "\n");
      }
    }
    printf("\nnSteps = %d \n", nSteps);
    printf("nSpks = %d\n", nTotSpks);
    fflush(vmFP);
    fclose(vmFP);
    fclose(outVars);
    fclose(isynapFP);
    fclose(rTotalFP);
    fclose(gbgrndFP);
    fclose(gEEEIFP);
    //***** FREE MEMORY *****//
    free_vector(expSum, 1, N_Neurons);
    free_vector(iSynap, 1, N_Neurons); 
    free_vector(gEI_E, 1, N_Neurons);   
    free_vector(gEI_I, 1, N_Neurons);   
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
    FreeSparseMat(sConMat);
}











