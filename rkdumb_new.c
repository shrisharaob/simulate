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
  int i, k, mNeuron, clmNo, localLastNSteps;
  double x, h, *vm, *vmOld;
  double *v, *vout, *dv;
  v = vector(1,nvar);
  vout = vector(1,nvar);
  dv = vector(1,nvar);
  vm = vector(1, N_Neurons);
  vmOld = vector(1, N_Neurons);
  //*** START ***//
  for (i=1;i<=nvar;i++) {
    v[i] = vstart[i];
    //    y[i][1] = v[i];
  }
  //*** TIMELOOP ***//
  xx[1] = x1;  
  x = x1;
  h = (x2 - x1) / nstep;
  localLastNSteps = nstep - STORE_LAST_N_STEPS;
  for (k = 1; k <= nstep; k++) {
    for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
      clmNo = (mNeuron - 1) * N_StateVars;
      vmOld[mNeuron] = v[1 + clmNo];
    }
    printf("%f\n", vmOld[1]);
    (*derivs)(x,v,dv);
    rk4(v,dv,nvar,x,h,vout,derivs);
    if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
    x += h; 
    
    /* RENAME */
    for (i=1;i<=nvar;i++) {
        v[i]=vout[i];
        //y[i][k+1] = v[i];
    }
    // detect spike - Zero crossing 
    if(k>2) {
      for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
        clmNo = (mNeuron - 1) * N_StateVars;
        if(k > localLastNSteps) {
          //          printf("%d\n", k);
          y[mNeuron][(k - localLastNSteps)] = v[1+clmNo];
          xx[k - localLastNSteps]  = x;
//          if(theta == 0) {
          //            fprintf(isynapFP, "%f %f ", tempCurE[mNeuron], tempCurI[mNeuron]);
            //          }
        }
        IF_SPK[mNeuron] = 0;
        vm[mNeuron] = v[1 + clmNo]; 
        if(v[1 + clmNo] > SPK_THRESH) { 
          if(vmOld[mNeuron] <= SPK_THRESH) {
            IF_SPK[mNeuron] = 1;
            fprintf(spkTimesFp, "%f %d\n", x, mNeuron);   
          }
        }
        

        //  fprintf(vmFP, "%f ", y[1+clmNo][k]);
        //          fprintf(outVars, "%f %f ", iSynap[mNeuron], iBg[mNeuron], iFF[mNeuron]);
      }
      //        fprintf(outVars, "\n");
      //      if(theta==0) {fprintf(isynapFP, "\n");}
      //        fprintf(vmFP, "\n");
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
  }
  free_vector(v,1,nvar);
  free_vector(vout,1,nvar);
  free_vector(dv,1,nvar);
  free_vector(vm, 1, N_Neurons);
}
#undef NRANSI

