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
  if(V != 10) { 
    out = 0.01 * (-V + 10) / (exp((-V + 10) / 10) - 1);
  }
  else {
    out = 0.1;
  }
  return out;
}

double beta_n(double V) {
  double out;
  out = 0.125 * exp(- V / 80);
  return out;
}

double alpha_m(double V) {
  double out;
  if(V != 25) { 
    out = 0.1 * (-V + 25) / (exp((-V + 25) / 10) - 1);
  }
  else {
    out = 1;
  }
  return out;
}

double beta_m(double V) {
  double out;
  out = 4 * exp(-V / 18);
  return out;
}

double alpha_h(double V) {
  double out;
  out = 0.07 * exp(- V / 20);
  return out;
}

double beta_h(double V) {
  double out;
  out = 1 / (exp((-V + 30) / 10) + 1);
  return out;
  }

extern double dt;
//stateVar = [V, n, m, h]
void derivs(double t, double stateVar[], double dydx[]) {
  int tIdx;
  tIdx = (int)(t / dt);
  dydx[1] =  1/Cm * (input_cur[tIdx] 
		     - G_Na * pow(stateVar[3], 3) * stateVar[4] * (stateVar[1] - E_Na) 
		     - G_K * pow(stateVar[2], 4) * (stateVar[1] - E_K) 
		     - G_L * (stateVar[1] - E_L));
  
  dydx[2] = alpha_n(stateVar[1]) * (1 - stateVar[2]) - beta_n(stateVar[1]) * stateVar[2];
  
  dydx[3] = alpha_m(stateVar[1]) * (1 - stateVar[3]) - beta_m(stateVar[1]) * stateVar[3];

  dydx[4] = alpha_h(stateVar[1]) * (1 - stateVar[4]) - beta_h(stateVar[1]) * stateVar[4];
}
