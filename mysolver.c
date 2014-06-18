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


//#include "nwEq2.h"
//void (*derive)(double, double *, double *);
// compile as - gcc nrutil.c rk4.c auxFunctions.c rkdumb.c -g mysolver.c -lm -o mysolver
// gcc ran1.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb.c -g mysolver.c -lm -o mysolver
// -g for gdb
//-lm reqd. for linking math.h


// GLOBAL VARS
extern double **y, *xx, *input_cur, *IF_SPK, *expSum, *iSynap, 
  *gaussNoiseE, *gaussNoise, theta, contrast, *gFF, *iFF, *rTotal, muE, muI;
extern FILE *spkTimesF, *outVars;
sparseMat *sConMat[N_Neurons + 1]; // index staring with 1
void main(int argc, char **argv) {
    // ***** DECLARATION *****//
  int dim = 4;
    double *vstart, *spkTimes;;
    double x1 = 0.0, // simulation start time
      x2 = 100.0, // simulation end time
      thetaStep = 0.0;
    int nSteps, nThetaSteps;
    int kNeuron, clmNo, loopIdx=0;
    long idem;
    FILE *fpElapsedTime;
    clock_t begin, end;
    int IF_SAVE = 1;
    // ***** INITIALIZATION *****//

    dt = DT;
    nSteps = (int)((x2 - x1) / dt);
    printf("\nnSteps = %d \n", nSteps);
    xx = vector(1, nSteps);
    y = matrix(1, N_StateVars * N_Neurons, 1, nSteps);
    input_cur = vector(1, nSteps);
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
    strcpy(filebase, FILEBASE);
    vmFP = fopen(strcat(filebase, "vm.csv"), "w");
    strcpy(filebase, FILEBASE);
    strcpy(filebase, FILEBASE);
    spkTimesFp = fopen(strcat(filebase, "spkTimes.csv"),"w");
    strcpy(filebase, FILEBASE);
    outVars = fopen(strcat(filebase, "outvars.csv"), "w");
    strcpy(filebase, FILEBASE);
    isynapFP = fopen(strcat(filebase, "isynapEI.csv"), "w");
    strcpy(filebase, FILEBASE);
    rTotalFP = fopen(strcat(filebase, "rTotal.csv"), "w");
    strcpy(filebase, FILEBASE);
    gbgrndFP = fopen(strcat(filebase, "gBg.csv"), "w");
    strcpy(filebase, FILEBASE);
    gEEEIFP = fopen(strcat(filebase, "gEEEI.csv"), "w");
    strcpy(filebase, FILEBASE);

    begin = clock();
    srand(time(NULL)); // set the seed for random number generator
    //    genConMat(); // Generate conection matrix
    GenConMat02();
    GenSparseConMat(sConMat);
    //    GenSparseConMatDisp(sConMat);
    AuxRffTotal(); /* auxillary function, generates random variables for the 
                      simulation run; which are used approximating FF input */
    if(thetaStep > 0) {
      thetaVec = vector(1, 360 / thetaStep);
      LinSpace(0, 360, thetaStep, thetaVec, &nThetaSteps); 
    }
    else {
      thetaVec = vector(1, 1);
      thetaVec[1] = 0;
      nThetaSteps = 1;
    }
    printf("theta = %f %d\n", thetaVec[1], nThetaSteps);
    /* /\********\/ */
    //    conMat[1][1] = 0; 
    //    conMat[1][2] = 0; 
    //for(loopIdx = 3; loopIdx <=N_Neurons; ++loopIdx) {
    //conMat[loopIdx][2] = 1;
    //    }
    // conMat[1][3] = 1; */
    //    conMat[2][1] = 1; 
    //    conMat[2][2] = 0; 
    /* conMat[2][3] = 1; */
    /* conMat[3][1] = 1; */
    /* conMat[3][2] = 0; */
    /* conMat[3][3] = 0; */
    // compute
    printf("\nNE = %d\n", NE);
    printf("NI = %d\n", NI);
    printf("K = %d\n", (int)K);
    printf("computing...\n");
    //***** PARSE INPUT ARGS *****//
    if(argc > 1) {
      theta = atof(argv[2]);
      contrast = atof(argv[3]);
      muE = atof(argv[4]);
      muI = atof(argv[5]);
      // current pulse - tStart, tStop, stepSize
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
        input_cur[loopIdx] = 0; 
      }
    }
    //***** INITIALIZE STATE VARIABLES *****//
    idem = -1 * time(NULL);
    for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
      clmNo =  (kNeuron - 1) * N_StateVars;
      vstart[1 + clmNo] = -70 +  40 * ran2(&idem); // Vm(0) ~ U(-70, -30)
      vstart[2 + clmNo] = 0.3176;
      vstart[3 + clmNo] = 0.1;
      vstart[4 + clmNo] = 0.5961;
    }
    //***** INTEGRATE *****//
    for(loopIdx = 1; loopIdx <= nThetaSteps; ++loopIdx) {
      theta = thetaVec[loopIdx];
      fprintf(spkTimesFp, "%f %f\n", theta, theta);
      rkdumb(vstart, N_StateVars * N_Neurons, x1, x2, nSteps, derivs);
      fprintf(spkTimesFp, "%d %d\n", 11, 11); // delimiters for thetas 
      fprintf(spkTimesFp, "%d %d\n", 73, 73);
    }
    printf("Done! \n");
    end = clock();
    printf("\n time spent integrating : %fs", (double)(end - begin) / CLOCKS_PER_SEC);
    fpElapsedTime = fopen("elapsedTime.csv", "a");
    fprintf(fpElapsedTime, "%f %d\n", (double)(end - begin) / CLOCKS_PER_SEC, totalNSpks);
    fclose(spkTimesFp);
    //***** SAVE TO DISK *****//
    if(IF_SAVE) {    
      for(loopIdx = 1; loopIdx <= nSteps; ++loopIdx) {
        //      if(loopIdx <= 2e4) {
        fprintf(vmFP, "%f ", xx[loopIdx]);
        for(kNeuron = 1; kNeuron <= N_Neurons; ++kNeuron) {
          clmNo =  (kNeuron - 1) * N_StateVars;
          // y = [t, V_m, n, z, h, I_input]
          fprintf(vmFP, "%f ", y[1 + clmNo][loopIdx]);
        }
        fprintf(vmFP, "\n");
      }
    }
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
