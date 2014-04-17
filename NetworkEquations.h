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

extern double dt;
// m_inf 
// stateVar = [V, n, z, h]
// z - gating variable of the adaptation current
void derivs(double t, double stateVar[], double dydx[]) {
  int tIdx;
  tIdx = (int)(t / dt);
  dydx[1] =  1/Cm * (input_cur[tIdx] 
		     - G_Na * pow(m_inf(stateVar[1]), 3) * stateVar[4] * (stateVar[1] - E_Na) 
		     - G_K * pow(stateVar[2], 4) * (stateVar[1] - E_K) 
		     - G_L_e * (stateVar[1] - E_L)
                     - G_adapt * stateVar[3] * (stateVar[1] - E_K));
  
  dydx[2] = alpha_n(stateVar[1]) * (1 - stateVar[2]) - beta_n(stateVar[1]) * stateVar[2];
  
  dydx[3] = 1 / Tau_adapt * (z_inf(stateVar[1]) - stateVar[3]);

  dydx[4] = alpha_h(stateVar[1]) * (1 - stateVar[4]) - beta_h(stateVar[1]) * stateVar[4];
}
