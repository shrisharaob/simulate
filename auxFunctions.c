#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "nrutil.h"
#include "nr.h"
#include "globalVars.h"
#include "varProtos.h"
#include "auxFuncProtos.h"

// recurrent synaptic current
void Isynap1(double *vm) {
  int kNeuron, mNeuron;
  double out; 
  FILE *gIIFP;
  //  gIIFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/gII", "a");
  iSynap = vector(1, N_Neurons); 
  for(kNeuron = 1; kNeuron <= N_Neurons; ++kNeuron) {
    expSum[kNeuron] = EXP_SUM * expSum[kNeuron]; // EXP_SUM = exp(- dt / tau_synap)
    if(IF_SPK[kNeuron] >= 1) {
      expSum[kNeuron] += 1;
    }
  }
  for(mNeuron = 1; mNeuron <= NE; ++mNeuron) {
    if(IF_SPK[mNeuron] == 1) {  // if mNeuron fired, then add expSum to all the neurons to which it projects
      for(kNeuron = 1; kNeuron <= NE; ++kNeuron) { // E --> E connections
        if(conMat[mNeuron][kNeuron] == 1) { // kNeuron has a synap input from mNeuron ?
          gEE[kNeuron] += expSum[mNeuron];
        }
      }
      for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) { // E --> I connections
        if(conMat[mNeuron][kNeuron] == 1) { // kNeuron has a synap input from mNeuron ?
          gIE[kNeuron - NE] += expSum[mNeuron];
        }
      }
    }
    else {
      for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
        gEE[kNeuron] = EXP_SUM * gEE[kNeuron];
      }
      for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) {
        gIE[kNeuron - NE] = EXP_SUM * gIE[kNeuron - NE];
      }
    }
  }
  for(mNeuron = NE + 1; mNeuron <= N_Neurons; ++mNeuron) {
    if(IF_SPK[mNeuron] == 1) {  // if mNeuron fired, then add expSum to all the neurons to which it projects
      for(kNeuron = 1; kNeuron <= NE; ++kNeuron) { // I --> E connections
        if(conMat[mNeuron][kNeuron] == 1) { // kNeuron has a synap input from mNeuron ?
          gEI[kNeuron] += expSum[mNeuron];
        }
      }
      for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) { // I --> I connections
        if(conMat[mNeuron][kNeuron] == 1) { // kNeuron has a synap input from mNeuron ?
          gII[kNeuron - NE] += expSum[mNeuron];
          //fprintf(gIIFP, "%f ", gII[kNeuron - NE]);
        }
      }
    }
    else {
      for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
        gEI[kNeuron] = EXP_SUM * gEI[kNeuron];
      }
      for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) {
        gII[kNeuron - NE] = EXP_SUM * gII[kNeuron - NE];
        //        
      }
    }

  }

  /* for(mNeuron = 1; mNeuron <= NE; ++mNeuron) { // ISynap for E neurons */
  /*   tempCurE[mNeuron] = 0; */
  /*   for(kNeuron = 1; kNeuron <= NE; ++kNeuron) { */
  /*     tempCurE[mNeuron] += -1 *  gEE[kNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_EE  */
  /*       * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E)); */
  /*   } */
  /*   tempCurI[mNeuron] = 0; */
  /*   for(kNeuron = 1; kNeuron <= NI; ++kNeuron) { */
  /*       tempCurI[mNeuron] += -1 * gEI[kNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_EI  */
  /*         * (RHO * (vm[mNeuron] - V_I) + (1 - RHO) * (E_L - V_I)); */
  /*   } */
  /*   iSynap[mNeuron] = tempCurE[mNeuron] + tempCurI[mNeuron]; */
  /* } */

  /* for(mNeuron = NE + 1; mNeuron <= N_Neurons; ++mNeuron) { // Isynap for I neurons */
  /*   tempCurI[mNeuron] = 0; */
  /*   for(kNeuron = 1; kNeuron <= NI; ++kNeuron) { */
  /*     tempCurI[mNeuron] += -1 * gII[kNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_II   */
  /*       * (RHO * (vm[mNeuron] - V_I) + (1 - RHO) * (E_L - V_I)); */
  /*   } */
  /*   tempCurE[mNeuron] = 0; */
  /*   for(kNeuron = 1; kNeuron <= NI; ++kNeuron) { */
  /*     tempCurE[mNeuron] += -1 * gIE[kNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_IE */
  /*       * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E)); */
  /*   } */
  /*   iSynap[mNeuron] = tempCurE[mNeuron] + tempCurI[mNeuron];   */
  /* } */

  for(mNeuron = 1; mNeuron <= NE; ++mNeuron) { // ISynap for E neurons
    tempCurE[mNeuron] = -1 *  gEE[mNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_EE
        * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
    tempCurI[mNeuron] = -1 * gEI[mNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_EI
          * (RHO * (vm[mNeuron] - V_I) + (1 - RHO) * (E_L - V_I));
      iSynap[mNeuron] = tempCurE[mNeuron] + tempCurI[mNeuron];
  }
  for(mNeuron = NE + 1; mNeuron <= N_Neurons; ++mNeuron) { // Isynap for I neurons
    tempCurI[mNeuron] = -1 * gII[mNeuron - NE] * (1/sqrt(K)) * INV_TAU_SYNAP * G_II
        * (RHO * (vm[mNeuron] - V_I) + (1 - RHO) * (E_L - V_I));
    tempCurE[mNeuron] = -1 * gIE[mNeuron -NE] * (1/sqrt(K)) * INV_TAU_SYNAP * G_IE
        * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
    iSynap[mNeuron] = tempCurE[mNeuron] + tempCurI[mNeuron];
    //    fprintf(gIIFP, "%f ", gII[mNeuron - NE]);
  }
  //  fprintf(gIIFP, "/n");
}


/* GENERATE CONNECTION MATRIX */
double XCordinate(int neuronIdx, double nA) {
  // nA - number of E or I cells
  double xCord;
  xCord = (1 - fmod((double)neuronIdx, sqrt(nA))) * L  / sqrt(nA);
  /* if(xCord < 0) {
    x = x + L;
  }
  else if(xCord >= L) {
    x = x - L;
  }*/
  return xCord ;
}
double YCordinate(int neuronIdx, double nA) {
  double yCord;
  yCord = (1 -fmod((double)neuronIdx, sqrt(nA))) * L / sqrt(nA);
  /*if(xCord < 0) {
    x = x + L;
  }
  else if(xCord >= L) {
    x = x - L;
    } */
 return yCord; 
}

void genConMat() {
  double **conProb, *zE, *zI, xDiff, yDiff, z1, denom, ranttt, tempZI;
  int clmId, rowId, IF_CONNECT;
  long idum = -1 * rand();
  FILE *conProbFP, *conMatFP;
  conProbFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/conProbMat", "w");
  conMatFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/conMatFp", "w");
  conProb = matrix(1, N_Neurons, 1, N_Neurons);
  zI = vector(1, N_Neurons);
  zE = vector(1, N_Neurons);
  z1 =   (1 / sqrt(2 * PI * CON_SIGMA));
  denom = (2 * CON_SIGMA * CON_SIGMA);
  // connection probablity for E cells
  for(rowId =1; rowId <= NE; ++rowId) {
    for(clmId = 1; clmId <= NE; ++clmId) {
      xDiff = XCordinate(rowId, NE) - XCordinate(clmId, NE);
      yDiff = YCordinate(rowId, NE) - YCordinate(clmId, NE);
      conProb[rowId][clmId] =  z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom)) * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      fprintf(conProbFP ,"%f ", conProb[rowId][clmId]);
    }
    for (clmId = NE + 1; clmId <= N_Neurons; ++clmId) {
      xDiff = XCordinate(rowId, NE) - XCordinate(clmId - NE, NI);
      yDiff = YCordinate(rowId, NE) - YCordinate(clmId - NE, NI);
      conProb[rowId][clmId] =  z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom))
                             * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      fprintf(conProbFP, "%f ", conProb[rowId][clmId]);
    }
    fprintf(conProbFP, "\n");
  }
  // connection probablity for I cells
  for(rowId = 1 + NE; rowId <= N_Neurons; ++rowId) {
    for(clmId = 1; clmId <= NE; ++clmId) {
      xDiff = XCordinate(rowId - NE, NI) - XCordinate(clmId, NE);
      yDiff = YCordinate(rowId - NE, NI) - YCordinate(clmId, NE);
      conProb[rowId][clmId] = z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom)) 
                             * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      fprintf(conProbFP ,"%f ", conProb[rowId][clmId]);
    }
    for (clmId = NE + 1; clmId <= N_Neurons; ++clmId) {
      xDiff = XCordinate(rowId - NE, NI) - XCordinate(clmId - NE, NI);
      yDiff = YCordinate(rowId - NE, NI) - YCordinate(clmId - NE, NI);
      conProb[rowId][clmId] = z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom))
                             * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      fprintf(conProbFP, "%f ", conProb[rowId][clmId]);
    }
    fprintf(conProbFP, "\n");
  }
  fclose(conProbFP);
  // compute pre-factor zB[clm] = K / sum(conProd(:, clm))
  for(clmId = 1; clmId <= N_Neurons; ++clmId) {
    for(rowId = 1; rowId <= NE; ++rowId) {
      zE[clmId] += conProb[rowId][clmId];
    }
    //    printf("b - %f ", zE[clmId]);
    zE[clmId] = (double)K / zE[clmId];
    //printf("a - %f \n", zE[clmId]);
    for(rowId = NE + 1; rowId <= N_Neurons; ++rowId) {
      zI[clmId] += conProb[rowId][clmId];
    }
    //    printf("%f ", zI[clmId]);
    zI[clmId] =(double)K /  zI[clmId];
    //    printf("%f \n", zI[clmId]);
  }
  // randomly connect neurons with probability given by conProb
  for(rowId = 1; rowId <= NE; ++rowId) {
    for(clmId = 1; clmId <= N_Neurons; ++clmId) {
//      srand(time(NULL));
      idum = -1 * rand();
      ranttt = ran1(&idum);
      //    printf("%f %ld \n", ranttt, idum);
      if(ranttt <=  zE[clmId] * conProb[rowId][clmId]) {
        conMat[rowId][clmId] = 1;
      }
      fprintf(conMatFP ,"%f ", conMat[rowId][clmId]);
    }
    fprintf(conMatFP, "\n");
  }
    srand(time(NULL));  
for(rowId = 1 + NE; rowId <= N_Neurons; ++rowId) {
  for(clmId = 1; clmId <= N_Neurons; ++clmId) {
    idum = -1 * rand();
    //  printf("%ld \n", idum);
    if(ran1(&idum) <=  zI[clmId] * conProb[rowId][clmId]) {
      conMat[rowId][clmId] = 1;
    }
    fprintf(conMatFP,"%f ", conMat[rowId][clmId]);
  }
  fprintf(conMatFP,"\n");
 }
 fclose(conMatFP);
}

// Background current
void IBackGrnd(double *vm) {
  double D = 1;
  double gE, gI;
  long idum;
  int kNeuron;
  //  srand(time(NULL));
  for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
    idum = -1 * rand();
    //    printf("%f \n", gasdev(&idum));
    gaussNoiseE[kNeuron]  = gaussNoiseE[kNeuron] + DT * ( D * gasdev(&idum) / SQRT_DT -  gaussNoiseE[kNeuron] * INV_TAU_SYNAP);
    gE = G_EB * K * (RB_E + sqrt(RB_E / K) * gaussNoiseE[kNeuron]);
    iBg[kNeuron] = -1 * gE * (RHO * (vm[kNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
    fprintf(gbgrndFP, "%f ", gE);
  }
  for(kNeuron = 1; kNeuron <= NI; ++kNeuron) {
    idum = -1 * rand();
    gaussNoiseI[kNeuron] = gaussNoiseI[kNeuron] + DT * ( D * gasdev(&idum) / SQRT_DT -  gaussNoiseI[kNeuron] * INV_TAU_SYNAP);
    gI = G_IB * K * (RB_I + sqrt(RB_I / K) * gaussNoiseI[kNeuron]);
    iBg[kNeuron + NE] = -1 * gI * (RHO * (vm[kNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
    fprintf(gbgrndFP,"%f ", gI);
  }
  fprintf(gbgrndFP,"\n");
}

// ff input 
void AuxRffTotal() {
  // draws random numbers from the specified distribution
  int lNeuron, i;
  long idem1, idem2, idem3, idem4; // seeds for rand generator
  for(lNeuron = 1; lNeuron <= N_Neurons; ++lNeuron) {
    idem1 = -1 * rand();
    idem2 = -1 * rand();
    randnXiA[lNeuron] =  gasdev(&idem1);
    randuDelta[lNeuron] = PI * ran1(&idem2);
    for(i = 1; i <= 4; ++i) {
      idem3 = -1 * rand();
      randwZiA[lNeuron][i] = 1.4142135 * sqrt(-1 * log(ran1(&idem3)));
    }
    for(i = 1; i <= 3; ++i) {
      idem4 = -1 * rand();
      randuPhi[lNeuron][i] = 2 * PI * ran1(&idem4);
    }
  }
}
// Eq. 16                                                                                                                                                                                                        
void RffTotal(double theta, double t) {
  double etaE, etaI;
  int lNeuron;
  //  if(t>DT) { 
    for(lNeuron = 1; lNeuron <= NE; ++lNeuron) {
      rTotalPrev[lNeuron] = rTotal[lNeuron]; // rTotal(t - 1)
      rTotal[lNeuron] = CFF * K * (R0 + R1 * log10(1 + contrast)) 
        + sqrt(CFF * K) * R0 * randnXiA[lNeuron]
        + sqrt(CFF * K) * R1 * log10(1 + contrast) * (randnXiA[lNeuron] 
                                                      + etaE * randwZiA[lNeuron][1] * cos(2 * (theta - randuDelta[lNeuron])) 
                                                      + muE * randwZiA[lNeuron][2] * cos(INP_FREQ * t - randuPhi[lNeuron][1])
                                                      + etaE * muE * 0.5 * (randwZiA[lNeuron][3] 
                                                      * cos(2 * theta + INP_FREQ * t - randuPhi[lNeuron][2])
                                                      + randwZiA[lNeuron][4] 
                                                      * cos(2 * theta - INP_FREQ * t + randuPhi[lNeuron][3])));
    }
    for(lNeuron = 1 + NE; lNeuron <= N_Neurons; ++lNeuron) {
      rTotalPrev[lNeuron] = rTotal[lNeuron]; // rTotal(t - 1)
      rTotal[lNeuron] = CFF * K * (R0 + R1 * log10(1 + contrast)) 
        + sqrt(CFF * K) * R0 * randnXiA[lNeuron]
        + sqrt(CFF * K) * R1 * log10(1 + contrast) * (randnXiA[lNeuron] 
                                                      + etaI * randwZiA[lNeuron][1] * cos(2 * (theta - randuDelta[lNeuron])) 
                                                      + muI * randwZiA[lNeuron][2] * cos(INP_FREQ * t - randuPhi[lNeuron][1])
                                                      + etaI * muI * 0.5 * (randwZiA[lNeuron][3] 
                                                      * cos(2 * theta + INP_FREQ * t - randuPhi[lNeuron][2])
                                                      + randwZiA[lNeuron][4] 
                                                      * cos(2 * theta - INP_FREQ * t + randuPhi[lNeuron][3])));

    }
}

void Gff(double theta, double t) {
  int kNeuron;
  double tempGasdev;
  long idem;
  if(t > DT) {
    for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
      idem = -1 * rand();
      ItgrlOld[kNeuron] = Itgrl[kNeuron];
      tempGasdev = gasdev(&idem);
      Itgrl[kNeuron] = rTotal[kNeuron] + sqrt(rTotal[kNeuron]) * tempGasdev; // / sqrt(dt)
      //      Itgrl[kNeuron] = cos(2 * PI * t) + 10 + sqrt(cos(2 * PI * t) + 10) * gasdev(&idem);
      //gFF[kNeuron] += DT * (- 0.5*(1/1000) * ( (1/K) * gFF[kNeuron] - Itgrl[kNeuron]));
      gFF[kNeuron] += DT * (-1 * GFF_E * sqrt(1/K) * INV_TAU_SYNAP 
                            * ( INV_TAU_SYNAP * gFF[kNeuron] - Itgrl[kNeuron]));
      fprintf(rTotalFP, "%f %f ", gFF[kNeuron], Itgrl[kNeuron]);
    }
    for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) {
      idem = -1 * rand();
      ItgrlOld[kNeuron] = Itgrl[kNeuron];
      Itgrl[kNeuron] = rTotal[kNeuron] + sqrt(rTotal[kNeuron]) * gasdev(&idem); // / SQRT_DT;
      gFF[kNeuron] += DT * (-1 * GFF_I * sqrt(1/K) * INV_TAU_SYNAP 
                            * ( INV_TAU_SYNAP * gFF[kNeuron] - Itgrl[kNeuron]));
      fprintf(rTotalFP, "%f %f ", gFF[kNeuron], rTotal[kNeuron]);
    }
    fprintf(rTotalFP, "\n");
  }
  else {
     for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
       idem = -1 * rand();
       Itgrl[kNeuron] = rTotal[kNeuron] + sqrt(rTotal[kNeuron]) * gasdev(&idem); // / SQRT_DT;
     }  
     for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) {
       idem = -1 * rand();
       Itgrl[kNeuron] = rTotal[kNeuron] + sqrt(rTotal[kNeuron]) * gasdev(&idem); // / SQRT_DT;
     }
  }
}
void IFF(double *vm) {
  int mNeuron;
  for (mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
    iFF[mNeuron] = -1 * gFF[mNeuron] * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
  }
}

// minimal conmat for balance
void GenConMat02() {
  int i, j;
  long idem;
  FILE *conMatFP;
  conMatFP = fopen("/home/shrisha/Documents/cnrs/results/network_model_outFiles/conMatFp", "w");
  for(i = 1; i <= NE + NI; ++i) {
    for(j = 1; j <= NE + NI; ++j) {
      if(i <= NE & j <= NE) {conMat[i][j] = 0;} // E --> E
      if (i <= NE & j > NE) {conMat[i][j] = 0;} // E --> I
      if (i > NE & j <= NE) {conMat[i][j] = 1;} // I --> E
      if (i > NE & j > NE) { // I --> I
        idem = -1 * rand();
        if((K / NI) >= ran1(&idem)) {
          conMat[i][j] = 1;
        }
      }
      fprintf(conMatFP,"%f ", conMat[i][j]);
    }
    fprintf(conMatFP, "\n");
    //    printf("%f %f %f %f", conMat[i][1], conMat[i][2], conMat[i][3], conMat[i][4]);
  }
  //  printf("\n");
  fclose(conMatFP);
}

void LinSpace(double startVal, double stopVal, double stepSize, double *outVector, int* nSteps)  {
  // generate equally spaced vector
  int loopId = 0;
  //  nSteps = &nVecSteps;
  if(stepSize > 0) {
    *nSteps = (int)(stopVal - startVal) / stepSize;
    //    outVector = vector(1, *nSteps);
    if(*nSteps > 1) {
      outVector[1] = startVal;
      for(loopId = 2; loopId <= *nSteps; ++loopId) {
        outVector[loopId] = outVector[loopId - 1] + stepSize;
      }
    }
  }
  else {
    //    outVector = vector(1, 2);
    outVector[1] = startVal;
    *nSteps = 1;
  }
}
