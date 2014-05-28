// synaptic input optimized network equations
#include <cuda.h>
#include "devHostConstants.h"
//extern double *input_cur, *iSynap, *iBg, *outvars;

__device__ float alpha_n(float vm);
__device__ float alpha_m(float vm);
__device__ float alpha_h(float vm);
__device__ float beta_n(float vm);
__device__ float beta_m(float vm);
__device__ float beta_h(float vm);

__device__ float  alpha_n(float vm) {
  float out;
  if(vm != -34) { 
    out = 0.1 * (vm + 34) / (1 - exp(-0.1 * (vm + 34)));
  }
  else {
    out = 0.1;
  }
  return out;
}

__device__ float beta_n(float vm) {
  float out;
  out = 1.25 * exp(- (vm + 44) / 80);
  return out;
}

__device__ float alpha_m(float vm) {
  float out;
  if(vm != -30) { 
    out = 0.1 * (vm + 30) / (1 - exp(-0.1 * (vm + 30)));
  }
  else {
    out = 1;
  }
  return out;
}

__device__ float beta_m(float vm) {
  float out;
  out = 4 * exp(-(vm + 55) / 18);
  return out;
}

__device__ float alpha_h(float vm) {
  float out;
  out = 0.7 * exp(- (vm + 44) / 20);
  return out;
}

__device__ float beta_h(float vm) {
  float out;
  out = 10 / (exp(-0.1 * (vm + 14)) + 1);
  return out;
  }

__device__ float m_inf(float vm) {
  float out, temp;
  temp = alpha_m(vm);
  out = temp / (temp + beta_m(vm));
  return out;
}

//z is the gating varible of the adaptation current
__device__ float z_inf(float(vm)) {
  float out;
  out = 1 / (1 + exp(-0.7 *(vm + 30)));
  return out;
}

extern float dt, *iSynap;
// m_inf 
// stateVar = [vm, n, z, h]
// z - gating variable of the adaptation current
__device__ void derivs(float t, float stateVar[], float dydx[]) {
  int tIdx, kNeuron, colNo;
  double cur = 0;
  tIdx = (int)(t / dt) + 1;
  for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
    colNo = (kNeuron - 1) * N_StateVars;
    /* if(kNeuron == 1 && t >= 30 && t <= 35) {  */
    /*  cur = 10;//input_cur[tIdx];  */
    /*  }  */
    /* else {cur = 0;}  */
    //     cur = 10;
    cur = 0.25 * sqrt(K);
    //       cur=2.8;
    //    printf("\n ICur : %f", cur);
    if (kNeuron <= NE) { 
      dydx[1 + colNo] =  1/Cm * (cur 
                                 - G_Na * pow(m_inf(stateVar[1 + colNo]), 3) * stateVar[4 + colNo] * (stateVar[1 + colNo] - E_Na) 
                                 - G_K * pow(stateVar[2 + colNo], 4) * (stateVar[1 + colNo] - E_K) 
                                 - G_L_E * (stateVar[1 + colNo] - E_L)
                                 - G_adapt * stateVar[3 + colNo] * (stateVar[1 + colNo] - E_K) + iSynap[kNeuron]);// iBg[kNeuron]);//+ iFF[kNeuron]); // N = [NE; NI]
      }
      else {
        dydx[1 + colNo] =  1/Cm * (cur  
                                   - G_Na * pow(m_inf(stateVar[1 + colNo]), 3) * stateVar[4 + colNo] * (stateVar[1 + colNo] - E_Na) 
                                   - G_K * pow(stateVar[2 + colNo], 4) * (stateVar[1 + colNo] - E_K) 
                                   - G_L_I * (stateVar[1 + colNo] - E_L)
                                   - G_adapt * stateVar[3 + colNo] * (stateVar[1 + colNo] - E_K) + iSynap[kNeuron]); // + iBg[kNeuron]);//+ iFF[kNeuron]); // N = [NE; NI]
      }
     
    dydx[2 + colNo] = alpha_n(stateVar[1 + colNo]) * (1 - stateVar[2 + colNo]) 
                      - beta_n(stateVar[1 + colNo]) * stateVar[2 + colNo];
  
    dydx[3 + colNo] = 1 / Tau_adapt * (z_inf(stateVar[1 + colNo]) - stateVar[3 + colNo]);
    
    dydx[4 + colNo] = alpha_h(stateVar[1 + colNo]) * (1 - stateVar[4 + colNo]) 
                      - beta_h(stateVar[1 + colNo]) * stateVar[4 + colNo];
  }

}
