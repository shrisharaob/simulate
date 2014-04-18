#include<math.h>
#include<stdio.h>

#define PI 3.14159265359

#define DT 0.025
#define SQRT_DT sqrt(DT)

#define Cm 1 //microF / cm^2
#define E_Na 55 //mV
#define E_K -90
#define E_L -65
#define G_Na 100 // mS/cm^2
#define G_K 40
#define G_L_E 0.05 // excitatory
#define G_L_I 0.1 // inhibitory
#define G_adapt 0.5

#define Tau_adapt 60 // in ms

// params network
#define N_StateVars 4
#define NI 5
#define NE 5
#define N_Neurons (NE + NI)
#define K 2
o
// params patch
#define L 1.0
#define CON_SIGMA (L / 5.0)

 // params synapse
#define TAU_SYNAP 3.0
#define INV_TAU_SYNAP (1 / TAU_SYNAP)
#define V_E 0
#define V_I -80
#define G_EE 0.15
#define G_EI 2
#define G_IE 0.45
#define G_II 3
#define EXP_SUM exp(-1 * DT / TAU_SYNAP)

// backgrund input
#define RB_E 20.0
#define RB_I 50.0
double *iBg, *gaussNoiseE, *gaussNoiseI;
#define G_EB (0.3 /sqrt(K))
#define G_IB (0.4 /sqrt(K))

double *input_cur, *IF_SPK, **conMat;
double *iSynap, *expSum, *gEE, *gEI, *gIE, *gII;
FILE *outVars, *spkTimesFp;

// ff input
#define CFF 0.1
#define KFF 1.0
#define GE_FF 0.95
#define GI_FF 1.26
#define R0 2.0
#define R1 20.0
#define INP_FREQ (4 * PI) // 
#define ETA_E 0.4
#define ETA_I 0.4
#define MU_E 0.1
#define MU_I 0.1
#define GFF_E 0.95
#define GFF_I 1.26


double contrast, theta;
double *gFF, *iFF, *rTotal, muE, muI,
  *randnXiA, // norm rand number
  **randwZiA, // weibul rand number
  *randuDelta, // uniform rand (0, PI)
  **randuPhi; // uniform rand (0, 2.PI)


#define RHO 0.5 // ratio - smatic / dendritic synapses

#define SPK_THRESH 1
