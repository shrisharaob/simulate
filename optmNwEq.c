
// synaptic input optimized network equations
#include <stdio.h>
#include <math.h>

extern double *input_cur, *iSynap, *iBg, *outvars;

double alpha_n(double vm);
double alpha_m(double vm);
double alpha_h(double vm);
double beta_n(double vm);
double beta_m(double vm);
double beta_h(double vm);

double alpha_n(double vm) {
  double out;
  if(vm != -34.0) { 
    out = 0.1 * (vm + 34.0) / (1 - exp(-0.1 * (vm + 34.0)));
  }
  else {
    out = 0.1;
  }
  return out;
}

double beta_n(double vm) {
  double out;
  out = 1.25 * exp(- (vm + 44.0) / 80.0);
  return out;
}

double alpha_m(double vm) {
  double out;
  if(vm != -30.0) { 
    out = 0.1 * (vm + 30.0) / (1 - exp(-0.1 * (vm + 30.0)));
  }
  else {
    out = 1.0;
  }
  return out;
}

double beta_m(double vm) {
  double out;
  out = 4.0 * exp(-(vm + 55.0) / 18.0);
  return out;
}

double alpha_h(double vm) {
  double out;
  out = 0.7 * exp(- (vm + 44.0) / 20.0);
  return out;
}

double beta_h(double vm) {
  double out;
  out = 10.0 / (exp(-0.1 * (vm + 14.0)) + 1.0);
  return out;
  }

double m_inf(double vm) {
  double out, temp;
  temp = alpha_m(vm);
  out = temp / (temp + beta_m(vm));
  return out;
}

//z is the gating varible of the adaptation current
double z_inf(double(vm)) {
  double out;
  out = 1 / (1 + exp(-0.7 *(vm + 30)));
  return out;
}

extern double dt, *iSynap;
// m_inf 
// stateVar = [vm, n, z, h]
// z - gating variable of the adaptation current
void derivs(double t, double stateVar[], double dydx[]) {
  int tIdx, kNeuron, colNo;
  double cur = 0.0, iffScale, iBgScale;
  double *vm;
  vm = vector(1, N_Neurons);
  iBgScale = 1.0;
  iffScale = 0.0;
  tIdx = (int)(t / dt) + 1;
  for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
    colNo = (kNeuron - 1) * N_StateVars;
    vm[kNeuron] = stateVar[1 + colNo];
    if(kNeuron == 2) { fprintf(gbgrndFP, "%f\n", iBg[kNeuron]); }
    if (kNeuron <= NE) { 
      if(kNeuron == -10 && kNeuron == -25 ) { /*Na 2+ current blocked */
        dydx[1 + colNo] =  1/Cm * (cur
                                 - G_K * pow(stateVar[2 + colNo], 4) * (stateVar[1 + colNo] - E_K)
                                 - G_L_E * (stateVar[1 + colNo] - E_L)
                                 - G_adapt * stateVar[3 + colNo] * (stateVar[1 + colNo] - E_K) + iSynap[kNeuron] + iBgScale* iBg[kNeuron] +  iffScale * iFF[kNeuron]);
      }
      else {
        dydx[1 + colNo] =  1/Cm * (cur 
                                 - G_Na * pow(m_inf(stateVar[1 + colNo]), 3) * stateVar[4 + colNo] * (stateVar[1 + colNo] - E_Na) 
                                 - G_K * pow(stateVar[2 + colNo], 4) * (stateVar[1 + colNo] - E_K) 
                                 - G_L_E * (stateVar[1 + colNo] - E_L)
                                 - G_adapt * stateVar[3 + colNo] * (stateVar[1 + colNo] - E_K) + iSynap[kNeuron] + iBgScale * iBg[kNeuron] +  iffScale * iFF[kNeuron]);
      }
    }
    else {
      iBg[kNeuron] = -1 * G_EB * K * RB_E * (RHO * vm[kNeuron] + (1 - RHO) * E_L);
      dydx[1 + colNo] =  1/Cm * (cur  
                                   - G_Na * pow(m_inf(stateVar[1 + colNo]), 3) * stateVar[4 + colNo] * (stateVar[1 + colNo] - E_Na) 
                                   - G_K * pow(stateVar[2 + colNo], 4) * (stateVar[1 + colNo] - E_K) 
                                   - G_L_I * (stateVar[1 + colNo] - E_L)
                                   - 0.0 * G_adapt * stateVar[3 + colNo] * (stateVar[1 + colNo] - E_K) + iSynap[kNeuron] + iBgScale * iBg[kNeuron] + iffScale * iFF[kNeuron]); 
      }
     
    dydx[2 + colNo] = alpha_n(stateVar[1 + colNo]) * (1 - stateVar[2 + colNo]) 
                      - beta_n(stateVar[1 + colNo]) * stateVar[2 + colNo];
    dydx[3 + colNo] = 1 / Tau_adapt * (z_inf(stateVar[1 + colNo]) - stateVar[3 + colNo]);
    dydx[4 + colNo] = alpha_h(stateVar[1 + colNo]) * (1 - stateVar[4 + colNo]) 
                      - beta_h(stateVar[1 + colNo]) * stateVar[4 + colNo];
  }
  
  free_vector(vm, 1, N_Neurons);
}
