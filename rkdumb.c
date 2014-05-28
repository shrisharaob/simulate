/* driver for the algorithm "rk4.c" */
// History:
//     Modified - Shrisha 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include "nrutil.h"
//#include "nr.h"
#include "globalVars.h"
#include "auxFuncProtos.h"
#include "cudaAuxFuncProtos.h"

extern FILE *spkTimesFp, *vmFP;
__device__ double **y, *tempCurE, *tempCurI;
__device__ void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep, void (*derivs)(double, double [], double [])) { 
  void rk4(double y[], double dydx[], int n,  double x, double h, double yout[], void (*derivs)(double, double [], double []));
  int i, k, mNeuron, clmNo;
  double x, h, vm, xx;
  double v[nvar], vout[nvar], dv[nvar];
  int spkNeuronId;

  //*** START ***//
  for (i = 0; i <= nvar - 1; i++) 
    {
      v[i] = vstart[i];
      y[i][1] = v[i];
    }
  //*** TIMELOOP ***//
  xx = x1;  
  x = x1;
  h = (x2 - x1) / nstep;
  for (k = 0; k <= nstep - 1; k++) 
    {
      nSpks = 0;
      (*derivs)(x,v,dv);
      rk4(v,dv,nvar,x,h,vout,derivs);
      if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
      x += h; 
      xx = x; //xx[k+1] = x;
      /* RENAME */
      for (i=1;i <= nvar - 1; i++) 
        {
          v[i]=vout[i];
          y[i][k+1] = v[i];
        }
      if(k > 2) {
        //        for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {

          clmNo = (mNeuron - 1) * N_StateVars;
          IF_SPK[mNeuron] = 0; // GLOBAL 
          vm = v[1 + clmNo]; 
          spkNeuronId = -1;
          if(v[1 + clmNo] > SPK_THRESH) { 
            if(y[1 + clmNo][k] <= SPK_THRESH) {
              IF_SPK[mNeuron] = 1;
              spkNeuronId = mNeuron;
              //      nSpks += 1;
              // Store SPKTIMES !!!!!!!!!!!!!!!!  fprintf(spkTimesFp, "%f %d\n", xx[k+1], mNeuron);
            }
          }
        }
      //      }
      CudaISynap(spkNeuronId); // allocate memore on device for spkNeuronId vector
      ISynapCudaAux(vm); // returns current 
      IBackGrnd(vm);
      // FF input current
      RffTotal(theta, x);
      Gff(theta, x);
      IFF(vm);
    }
}

