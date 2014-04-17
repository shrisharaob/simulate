#include<math.h>
#include<stdio.h>

#define Cm 1 //microF / cm^2
#define E_Na 55 //mV
#define E_K -90
#define E_L -65
#define G_Na 100 // mS/cm^2
#define G_K 40
#define G_L_e 0.05 // excitatory
#define G_L_i 0.1 // inhibitory
#define G_adapt 0.5

#define Tau_adapt 60 // in ms

// params network
#define N_StateVars 4
#define N_Neurons 2
#define K 2.0

 // params synapse
#define Tau_syn 3.0
#define V_E 0
#define V_I -80
#define G_EE 0.5

double *input_cur, *IF_SPK;
FILE *spkTimesFp;

#define RHO 0.5 // ratio - smatic / dendritic synapses

//#define J[2][2] {{0, 1},{0, 0}}

#define SPK_THRESH 1
