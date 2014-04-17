#include <stdio.h>
#include <math.h>

extern double *input_cur;

double alpha_n(double V);
double alpha_m(double V);
double alpha_h(double V);
double beta_n(double V);
double beta_m(double V);
double beta_h(double V);

double alpha_n(double V) {
  double out;
  if(V != -34) { 
    out = 0.1 * (V + 34) / (1 - exp(-0.1 * (V + 34)));
  }
  else {
    out = 0.1;
  }
  return out;
}

double beta_n(double V) {
  double out;
  out = 1.25 * exp(- (V + 44) / 80);
  return out;
}

double alpha_m(double V) {
  double out;
  if(V != -30) { 
    out = 0.1 * (V + 30) / (1 - exp(-0.1 * (V + 30)));
  }
  else {
    out = 1;
  }
  return out;
}

double beta_m(double V) {
  double out;
  out = 4 * exp(-(V + 55) / 18);
  return out;
}

double alpha_h(double V) {
  double out;
  out = 0.7 * exp(- (V + 44) / 20);
  return out;
}

double beta_h(double V) {
  double out;
  out = 10 / (exp(-0.1 * (V + 14)) + 1);
  return out;
  }

double m_inf(double V) {
  double out, temp;
  temp = alpha_m(V);
  out = temp / (temp + beta_m(V));
  return out;
}

//z is the gating varible of the adaptation current
double z_inf(double(V)) {
  double out;
  out = 1 / (1 + exp(-0.7 *(V + 30)));
  return out;
}

double IsynapEE(double t, double spkTimes[], int nSpikes, int weight, double V) {
  int kSpk = 0;
  double out = 0;
  FILE *fp;
  if(weight > 0) {
    for(kSpk = 0; kSpk < nSpikes; ++kSpk) {
      if(spkTimes[kSpk] < t) {
	out += exp(-(t - spkTimes[kSpk]) / Tau_syn);
      }
    }
  }
  fp = fopen("isyynap", "a");
  fprintf(fp, "%f ", out);
  out = out * (1/sqrt(K)) * (1 / Tau_syn) * G_EE * (RHO * (V - V_E) + (1 - RHO) * (E_L - V_E));
   fprintf(fp, "%f\n", out);
  fclose(fp);
  return out;
}

extern double dt;
// m_inf 
// stateVar = [V, n, z, h]
// z - gating variable of the adaptation current
void derivs(double t, double stateVar[], double dydx[]) {
  int tIdx, kNeuron, colNo;
  double cur;
  double nSpkss = 7;
  double spkss[] = {50, 53, 56, 59, 60, 62, 65};
  tIdx = (int)(t / dt) + 1;
  for(kNeuron = 1; kNeuron < N_Neurons + 1; ++kNeuron) {
    colNo = (kNeuron - 1) * N_StateVars;
    if(kNeuron == 1) { cur = input_cur[tIdx];}
    else {cur = 0;}
    dydx[1 + colNo] =  1/Cm * (cur 
		       - G_Na * pow(m_inf(stateVar[1 + colNo]), 3) * stateVar[4 + colNo] * (stateVar[1 + colNo] - E_Na) 
		       - G_K * pow(stateVar[2 + colNo], 4) * (stateVar[1 + colNo] - E_K) 
		       - G_L_e * (stateVar[1 + colNo] - E_L)
		       - G_adapt * stateVar[3 + colNo] * (stateVar[1 + colNo] - E_K)
		       - IsynapEE(t, spkss, nSpkss, 1, stateVar[1 + colNo]));
  
    dydx[2 + colNo] = alpha_n(stateVar[1 + colNo]) * (1 - stateVar[2 + colNo]) - beta_n(stateVar[1 + colNo]) * stateVar[2 + colNo];
  
    dydx[3 + colNo] = 1 / Tau_adapt * (z_inf(stateVar[1 + colNo]) - stateVar[3 + colNo]);

    dydx[4 + colNo] = alpha_h(stateVar[1 + colNo]) * (1 - stateVar[4 + colNo]) - beta_h(stateVar[1 + colNo]) * stateVar[4 + colNo];
  }
}
