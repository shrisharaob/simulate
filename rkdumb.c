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
  int i, k, mNeuron, clmNo;
  double x, h, *vm;
  double *v, *vout, *dv;
  v = vector(1,nvar);
  vout = vector(1,nvar);
  dv = vector(1,nvar);
  vm = vector(1, N_Neurons);
  //*** START ***//
  for (i=1;i<=nvar;i++) 
    {
      v[i] = vstart[i];
      y[i][1] = v[i];
    }
  //*** TIMELOOP ***//
  xx[1] = x1;  
  x = x1;
  h = (x2 - x1) / nstep;
  /* for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) { */
  /*   fprintf(outVars, "%f %f %f %f ", x, iSynap[mNeuron], iBg[mNeuron], iFF[mNeuron]); */
  /* } */
  //  fprintf(outVars, "\n");
  for (k = 1; k <= nstep; k++) 
    {
      (*derivs)(x,v,dv);
      rk4(v,dv,nvar,x,h,vout,derivs);
      if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
      x += h; 
      xx[k+1] = x;
      /* RENAME */
      for (i=1;i<=nvar;i++) 
        {
          v[i]=vout[i];
          y[i][k+1] = v[i];
        }
      //      fprintf(vmFP, "%f ", xx[k]);
      // detect spike - three consequtive time points are checked -->  ../^\..
      if(k>3) {
        //fprintf(outVars, "%f", x); // first clm  is time 
        for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
          clmNo = (mNeuron - 1) * N_StateVars;
          IF_SPK[mNeuron] = 0;
          vm[mNeuron] = v[1 + clmNo]; 
          if(v[1 + clmNo] > SPK_THRESH) { 
            if(y[1 + clmNo][k] <= SPK_THRESH) {
              IF_SPK[mNeuron] = 1;
              fprintf(spkTimesFp, "%f;%d\n", xx[k+1], mNeuron);
            }
          }
          //  fprintf(vmFP, "%f ", y[1+clmNo][k]);
          //          fprintf(outVars, "%f %f %f ", iSynap[mNeuron], iBg[mNeuron], iFF[mNeuron]);
          //          fprintf(isynapFP, "%f %f ", tempCurE[mNeuron], tempCurI[mNeuron]);
        }
        //        fprintf(outVars, "\n");
        //        fprintf(isynapFP, "\n");
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
/* (C) Copr. 1986-92 Numerical Recipes Software 9,)5. */
