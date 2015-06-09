/* driver for the algorithm "rk4.c" */
// History:
//     Numerical Recipes Software 9 
//     Modified - Shrisha 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nr.h"
#include "globalVars.h"
#include "auxFuncProtos.h"

extern FILE *spkTimesFp, *vmFP;
double **y, *xx;
void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep, void (*derivs)(double, double [], double [])) { 
  void rk4(double y[], double dydx[], int n,  double x, double h, double yout[], void (*derivs)(double, double [], double []));
  int i, k, mNeuron, clmNo, lastNStepsToStore, logVar;
  double x, h, *vm, *vmOld;
  double *v, *vout, *dv;
  int spksE = 0, spksI = 0;
  FILE *fpIFR = fopen("InstFR.csv", "w");
  v = vector(1,nvar); /*nvar = N_Neurons * N_StateVars  */
  vout = vector(1,nvar);
  dv = vector(1,nvar);
  vm = vector(1, N_Neurons);
  vmOld = vector(1, N_Neurons);
  lastNStepsToStore = nstep - STORE_LAST_N_STEPS;
  /* START */
  for (i=1;i<=nvar;i++) {
    v[i] = vstart[i];
  }
  /* TIME LOOP */
  x = x1;
  h = (x2 - x1) / nstep;
  logVar = nstep / 10;
  for (k = 1; k <= nstep; k++) {
    if(k % logVar == 0) {
      printf("%2.1f %% done \n", (double)k * 100.0  / (double)nstep);
    } 
    for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
        clmNo = (mNeuron - 1) * N_StateVars;
        vmOld[mNeuron] = v[1 + clmNo];   
      }
    //    printf("%f ", vmOld[1]);
    for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
        clmNo = (mNeuron - 1) * N_StateVars;
        if(k >= lastNStepsToStore) {
          y[mNeuron][k - lastNStepsToStore] = v[1+clmNo];
          xx[k - lastNStepsToStore] = x;
        }
    }

    (*derivs)(x,v,dv);
    rk4(v,dv,nvar,x,h,vout,derivs);
    if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
    x += h; 
    for (i=1;i<=nvar;i++) {
      v[i]=vout[i];
    }
      
    // detect spike - zero crossing
    if(k > 2) {
      for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
        clmNo = (mNeuron - 1) * N_StateVars;
        /*        if(k > lastNStepsToStore) {
          y[mNeuron][k - lastNStepsToStore] = v[1+clmNo];
          xx[k - lastNStepsToStore] = x;
          }*/
        IF_SPK[mNeuron] = 0;
        vm[mNeuron] = v[1 + clmNo]; 
        //        printf("%f \n", vm[mNeuron]);
        if(vm[mNeuron] > SPK_THRESH) { 
          if(vmOld[mNeuron] <= SPK_THRESH) {
            IF_SPK[mNeuron] = 1;
            fprintf(spkTimesFp, "%f;%d\n", x, mNeuron);
            nTotSpks += 1;
          }
        }
        if(mNeuron >= 5 && mNeuron < 25) {
          fprintf(isynapFP, "%f %f ", tempCurE[mNeuron], tempCurI[mNeuron]);
        }
      } 
      fprintf(isynapFP, "\n");
    }
    /* /\* else { *\/ */
    /* /\*   for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) { *\/ */
    /* /\*     clmNo = (mNeuron - 1) * N_StateVars; *\/ */
    /* /\*     fprintf(vmFP, "%f ", y[1+clmNo][k]); *\/ */
    /* /\*     //fprintf(outVars, "%f %f %f %f ", x, iSynap[mNeuron], iBg[mNeuron], iFF[mNeuron]); *\/ */
    /* /\*     //          fprintf(isynapFP, "%f %f ", tempCurE[mNeuron], tempCurI[mNeuron]); *\/ */
    /* /\*   } *\/ */
    /* /\*   //fprintf(outVars, "\n"); *\/ */
    /* /\*   //fprintf(isynapFP, "\n");  *\/ */
    /* /\*   fprintf(vmFP, "\n"); *\/ */
    /* /\* }  */
    // compute synaptic current
    Isynap1(vm);
    //compute background current
    IBackGrnd(vm);
    // FF input current
    RffTotal(theta, x);
    Gff(theta, x);
    IFF(vm); 

    for(i = 1; i <= N_NEURONS; ++i) {
      if(IF_SPK[i]) {
        if(i < NE) {
          spksE += 1;
        }
        else{
          spksI += 1;
        }
      }
    }
    if(!(k%(int)(50.0/DT))) {
      fprintf(fpIFR, "%f %f \n", ((double)spksE) / (DT * (double)NE), ((double)spksI) / (0.05 * (double)NI));fflush(fpIFR);
      fprintf(stdout, "%f %f \n", ((double)spksE) / (DT * (double)NE), ((double)spksI) / (0.05 * (double)NI));
      spksE = 0; 
      spksI = 0;
    }
  }
  fclose(fpIFR);
  free_vector(v,1,nvar);
  free_vector(vout,1,nvar);
  free_vector(dv,1,nvar);
  free_vector(vm, 1, N_Neurons);
}
#undef NRANSI
