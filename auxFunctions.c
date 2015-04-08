#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "config.h"
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
  for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) {
      gEI_E[mNeuron] *= EXP_SUM;
      gEI_I[mNeuron] *= EXP_SUM;
      if(IF_SPK[mNeuron] == 1) {  
        for(kNeuron = 1; kNeuron <= sConMat[mNeuron]->nPostNeurons; ++kNeuron) { 
          if(mNeuron <= NE) {       
            gEI_E[sConMat[mNeuron]->postNeuronIds[kNeuron]] += 1;
          }
          else
            gEI_I[sConMat[mNeuron]->postNeuronIds[kNeuron]] += 1;
        }
      }
  }
  for(mNeuron = 1; mNeuron <= N_Neurons; ++mNeuron) { // ISynap for E neurons
    if(mNeuron <=NE) {
      tempCurE[mNeuron] = -1 *  gEI_E[mNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_EE
                          * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
      tempCurI[mNeuron] = -1 * gEI_I[mNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_EI
                          * (RHO * (vm[mNeuron] - V_I) + (1 - RHO) * (E_L - V_I));
    }
    else {
      tempCurE[mNeuron] = -1 * gEI_E[mNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_IE
                          * (RHO * (vm[mNeuron] - V_E) + (1 - RHO) * (E_L - V_E));
      tempCurI[mNeuron] = -1 * gEI_I[mNeuron] * (1/sqrt(K)) * INV_TAU_SYNAP * G_II
                          * (RHO * (vm[mNeuron] - V_I) + (1 - RHO) * (E_L - V_I));
    }
    iSynap[mNeuron] = tempCurE[mNeuron] + tempCurI[mNeuron];
  }
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
  strcpy(filebase, FILEBASE);
  conProbFP = fopen(strcat(filebase,"conProbMat.csv"), "w");
  strcpy(filebase, FILEBASE);
  conMatFP = fopen(strcat(filebase,"conMat.csv"), "w");
  strcpy(filebase, FILEBASE);
  conProb = matrix(1, N_Neurons, 1, N_Neurons);
  strcpy(filebase, FILEBASE);
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
      //      fprintf(conProbFP ,"%f ", conProb[rowId][clmId]);
    }
    for (clmId = NE + 1; clmId <= N_Neurons; ++clmId) {
      xDiff = XCordinate(rowId, NE) - XCordinate(clmId - NE, NI);
      yDiff = YCordinate(rowId, NE) - YCordinate(clmId - NE, NI);
      conProb[rowId][clmId] =  z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom))
                             * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      //      fprintf(conProbFP, "%f ", conProb[rowId][clmId]);
    }
    //    fprintf(conProbFP, "\n");
  }
  // connection probablity for I cells
  for(rowId = 1 + NE; rowId <= N_Neurons; ++rowId) {
    for(clmId = 1; clmId <= NE; ++clmId) {
      xDiff = XCordinate(rowId - NE, NI) - XCordinate(clmId, NE);
      yDiff = YCordinate(rowId - NE, NI) - YCordinate(clmId, NE);
      conProb[rowId][clmId] = z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom)) 
                             * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      //      fprintf(conProbFP ,"%f ", conProb[rowId][clmId]);
    }
    for (clmId = NE + 1; clmId <= N_Neurons; ++clmId) {
      xDiff = XCordinate(rowId - NE, NI) - XCordinate(clmId - NE, NI);
      yDiff = YCordinate(rowId - NE, NI) - YCordinate(clmId - NE, NI);
      conProb[rowId][clmId] = z1 * z1 * exp(-1 * pow(fmod(xDiff, L), 2) / (denom))
                             * z1 * z1 * exp(-1 * pow(fmod(yDiff, L), 2) / (denom));
      //      fprintf(conProbFP, "%f ", conProb[rowId][clmId]);
    }
    //    fprintf(conProbFP, "\n");
  }
  fclose(conProbFP);
  // compute pre-factor zB[clm] = K / sum(conProd(:, clm))
  for(clmId = 1; clmId <= N_Neurons; ++clmId) {
    for(rowId = 1; rowId <= NE; ++rowId) {
      zE[clmId] += conProb[rowId][clmId];
    }
    zE[clmId] = (double)K / zE[clmId];
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
      ranttt = gsl_rng_uniform(gslRngState);
      //    printf("%f %ld \n", ranttt, idum);
      if(ranttt <=  zE[clmId] * conProb[rowId][clmId]) {
        conMat[rowId][clmId] = 1;
      }
      //      fprintf(conMatFP ,"%f ", conMat[rowId][clmId]);
    }
    //    fprintf(conMatFP, "\n");
  }
    srand(time(NULL));  
for(rowId = 1 + NE; rowId <= N_Neurons; ++rowId) {
  for(clmId = 1; clmId <= N_Neurons; ++clmId) {
    if(gsl_rng_uniform(gslRngState) <=  zI[clmId] * conProb[rowId][clmId]) {
      conMat[rowId][clmId] = 1;
    }
    //    fprintf(conMatFP,"%f ", conMat[rowId][clmId]);
  }
  //  fprintf(conMatFP,"\n");
 }
 fclose(conMatFP);
}

// Background current
void IBackGrnd(double *vm) {
  double D = 1;
  double gE, gI;
  long idum = -1 * time(NULL);
  int kNeuron;
  /*  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  tmp1 = 1 - DT * INV_TAU_SYNAP;
  tmp2 = SQRT_DT * INV_TAU_SYNAP;
  tmp3 = G_EB * K * RB_E ;
  tmp4 = G_IB * K * RB_I ;
  tmp5 = G_EB * K * sqrt(RB_E / K);
  tmp6 = G_IB * K * sqrt(RB_I / K);
*/  
  for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
    gaussNoiseE[kNeuron] = gaussNoiseE[kNeuron] * (1 - DT * INV_TAU_SYNAP) + SQRT_DT  * INV_TAU_SYNAP * gasdev(&idum);
    gaussNoiseE[kNeuron] = 0;
    gE = G_EB * K * (RB_E + sqrt(RB_E / K) * gaussNoiseE[kNeuron]);
    iBg[kNeuron] = -1 * gE * (RHO * (vm[kNeuron] - V_E) + (1 - RHO) * (E_L - V_E));

  }
  for(kNeuron = 1; kNeuron <= NI; ++kNeuron) {
    /*    gaussNoiseI[kNeuron] = gaussNoiseI[kNeuron] + DT * ( D * gasdev(&idum) / SQRT_DT -  gaussNoiseI[kNeuron] * INV_TAU_SYNAP);*/
    /*    gaussNoiseI[kNeuron] = gaussNoiseI[kNeuron] * tmp1 + tmp2 * gasdev(&idum);*
    gaussNoiseI[kNeuron] = 0;
    gI = tmp4 + tmp6 * gaussNoiseI[kNeuron];
    iBg[kNeuron + NE] = -1 * gI * (RHO * (vm[kNeuron] - V_E) + (1 - RHO) * (E_L - V_E));*/
    //    fprintf(gbgrndFP,"%f ", gI);

    gaussNoiseI[kNeuron] = gaussNoiseI[kNeuron] * (1 - DT * INV_TAU_SYNAP) + SQRT_DT  * INV_TAU_SYNAP * gasdev(&idum);
    gaussNoiseI[kNeuron] = 0;
    gI = G_IB * K * (RB_I + sqrt(RB_I / K) * gaussNoiseI[kNeuron]);
    iBg[kNeuron] = -1 * gI * (RHO * (vm[kNeuron + NE] - V_E) + (1 - RHO) * (E_L - V_E));
  }
  //  fprintf(gbgrndFP,"\n");
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
    randuDelta[lNeuron] = PI * gsl_rng_uniform(gslRngState);
    for(i = 1; i <= 4; ++i) {
      randwZiA[lNeuron][i] = 1.4142135 * sqrt(-1 * log(gsl_rng_uniform(gslRngState)));
    }
    for(i = 1; i <= 3; ++i) {
      idem4 = -1 * rand();
      randuPhi[lNeuron][i] = 2 * PI * gsl_rng_uniform(gslRngState);
    }
  }
}
// Eq. 16                                                                                                                                                                                                        
void RffTotal(double theta, double t) {
  int lNeuron;
  //  if(t>DT) {
    for(lNeuron = 1; lNeuron <= NE; ++lNeuron) {
      rTotalPrev[lNeuron] = rTotal[lNeuron]; // rTotal(t - 1)
      rTotal[lNeuron] = CFF * K * (R0 + R1 * log10(1 + contrast)) 
        + sqrt(CFF * K) * R0 * randnXiA[lNeuron]
        + sqrt(CFF * K) * R1 * log10(1 + contrast) * (randnXiA[lNeuron] 
                                                      + ETA_E * randwZiA[lNeuron][1] * cos(2 * (theta - randuDelta[lNeuron])) 
                                                      + MU_E * randwZiA[lNeuron][2] * cos(INP_FREQ * t - randuPhi[lNeuron][1])
                                                      + ETA_E * MU_E * 0.5 * (randwZiA[lNeuron][3] 
                                                      * cos(2 * theta + INP_FREQ * t - randuPhi[lNeuron][2])
                                                      + randwZiA[lNeuron][4] 
                                                      * cos(2 * theta - INP_FREQ * t + randuPhi[lNeuron][3])));
    }
    for(lNeuron = 1 + NE; lNeuron <= N_Neurons; ++lNeuron) {
      rTotalPrev[lNeuron] = rTotal[lNeuron]; // rTotal(t - 1)
      rTotal[lNeuron] = CFF * K * (R0 + R1 * log10(1 + contrast)) 
        + sqrt(CFF * K) * R0 * randnXiA[lNeuron]
        + sqrt(CFF * K) * R1 * log10(1 + contrast) * (randnXiA[lNeuron] 
                                                      + ETA_I * randwZiA[lNeuron][1] * cos(2 * (theta - randuDelta[lNeuron])) 
                                                      + MU_I * randwZiA[lNeuron][2] * cos(INP_FREQ * t - randuPhi[lNeuron][1])
                                                      + ETA_I * MU_I * 0.5 * (randwZiA[lNeuron][3] 
                                                      * cos(2 * theta + INP_FREQ * t - randuPhi[lNeuron][2])
                                                      + randwZiA[lNeuron][4] 
                                                      * cos(2 * theta - INP_FREQ * t + randuPhi[lNeuron][3])));

    }
}

void Gff(double theta, double t) {
  int kNeuron;
  double tempGasdev;
  long idem = -1;
  if(t > DT) {
    for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
      /*      ItgrlOld[kNeuron] = Itgrl[kNeuron];*/
      tempGasdev = gasdev(&idem);
      Itgrl[kNeuron] = Itgrl[kNeuron] * (1 - DT * INV_TAU_SYNAP) + (SQRT_DT * INV_TAU_SYNAP) * tempGasdev;
      gFF[kNeuron] = GFF_E * (rTotal[kNeuron] + sqrt(rTotal[kNeuron] * Itgrl[kNeuron]));
      //      fprintf(rTotalFP, "%f %f ", gFF[kNeuron], Itgrl[kNeuron]);
    }
    for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) {
      Itgrl[kNeuron] = Itgrl[kNeuron] * (1 - DT * INV_TAU_SYNAP) + (SQRT_DT * INV_TAU_SYNAP) * tempGasdev;
      gFF[kNeuron] = GFF_I * (rTotal[kNeuron] + sqrt(rTotal[kNeuron] * Itgrl[kNeuron]));
    }
  }
  else {
    for(kNeuron = 1; kNeuron <= NE; ++kNeuron) {
       Itgrl[kNeuron] = 0.0;
     }  
     for(kNeuron = NE + 1; kNeuron <= N_Neurons; ++kNeuron) {
       Itgrl[kNeuron] = 0.0;
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
void GenConMat02(int IF_GEN) {
  int i, j;
  FILE *conMatFP;

  if(IF_GEN) {
    conMatFP = fopen(strcat(filebase,"conMat.dat"), "wb");
    fprintf(stdout, "generating conMat..."); fflush(stdout);
    for(i = 1; i <= NE + NI; ++i) {
      for(j = 1; j <= NE + NI; ++j) {
        if(i <= NE & j <= NE & NE > 0) { // E --> E
          if((K /(double)NE) >= gsl_rng_uniform(gslRngState)) {
            conMat[i][j] = 1;
          }
        } 
        if (i <= NE & j > NE & NE > 0) { // E --> I
          if((K / (double)NE) >= gsl_rng_uniform(gslRngState)) {
            conMat[i][j] = 1;
          }
        } 
        if (i > NE & j <= NE & NI > 0) { // I --> E
          if((K / (double)NI) >= gsl_rng_uniform(gslRngState)) {
            conMat[i][j] = 1;
          }
        } 
        if (i > NE & j > NE & NI > 0) { // I --> I
          if((K / (double)NI) >= gsl_rng_uniform(gslRngState)) {
            conMat[i][j] = 1;
          }
        }
        /*      fprintf(conMatFP,"%f ", conMat[i][j]);*/
      }
      /*    fprintf(conMatFP, "\n");*/
    }
    fprintf(stdout, "done\n");
    fprintf(stdout, "writing to file..."); fflush(stdout);
    fwrite(&conMat[1][1], sizeof(double), N_Neurons * N_Neurons, conMatFP);
  }
  else {
    conMatFP = fopen(strcat(filebase,"conMat.dat"), "rb");
    fprintf(stdout, "reading file... "); fflush(stdout);
    fread(&conMat[1][1], sizeof(double), N_Neurons * N_Neurons, conMatFP);
  }
  fprintf(stdout, "done\n");
  fclose(conMatFP);
  
  FILE  *fp01 = fopen("countI.csv", "w"),  *fp02 = fopen("countE.csv", "w");
  printf("\nN = %d\n", N_NEURONS);
  int countE = 0, countI = 0;
  for(i = 1; i <= N_NEURONS; ++i) {
    countI = 0;
    countE = 0;
    for(j = 1; j <= N_NEURONS; ++j) {
      if(j < NE) {
        countE += conMat[i][j];
      }
      else {
        countI += conMat[i][j];
      }
    }
    fprintf(fp02, "%d\n", countE); 
    fprintf(fp01, "%d\n", countI);
  }
  fprintf(stdout, " done\n");
  fclose(fp02);   
  fclose(fp01);
}

void GenSparseConMat(sparseMat *sPtr[]) {
  int row, clm, *count, loopid, *idx, tmpCntr;
  count = ivector(1, N_Neurons);
  for(row = 1; row <= N_Neurons; ++row) {
    sPtr[row] = (sparseMat *)malloc(sizeof(sparseMat));
    sPtr[row]->neuronId = row;
    //    printf("neurn %d connects to : ", sPtr[row]->neuronId);
    count[row] = 0;
    for(clm = 1; clm <= N_Neurons; ++clm) {
      if(conMat[row][clm] == 1) {
        count[row] += 1;
      }
    }
    //printf("\n count = %d \n ", count[row]);
    sPtr[row]->nPostNeurons = count[row];
    sPtr[row]->postNeuronIds = ivector(1, count[row]);
    idx = ivector(1, count[row]);
    tmpCntr = 0;
    for(clm = 1; clm <= N_Neurons; ++clm) {
      if(conMat[row][clm] == 1) {
        tmpCntr+=1;
        idx[tmpCntr] = clm;
        //        printf("idx = %d ", idx[tmpCntr]);
      }
    }
    for(loopid = 1; loopid <= count[row]; ++loopid) { 
      sPtr[row]->postNeuronIds[loopid] = idx[loopid];
      //printf(" %d ", sPtr[row]->postNeuronIds[loopid]);
      //      printf(" idx = %d", idx[loopid]);
    }
    //        printf("\n");
    free_ivector(idx, 1, count[row]);
  }
  free_ivector(count, 1, N_Neurons);
}

void FreeSparseMat(sparseMat *sPtr[]) {
  int kNeuron;
  for(kNeuron = 1;  kNeuron<= N_Neurons; ++kNeuron) {
      free(sPtr[kNeuron]->postNeuronIds);
      free(sPtr[kNeuron]);
  }
}
/*
void ReadSparseConMat() {
  FILE *fpSparseConVec, *fpIdxVec, *fpNpostNeurons;
  fpSparseConVec = fopen("sparseConVec.dat", "rb");
  fpIdxVec = fopen("idxVec.dat", "rb");
  fpNpostNeurons = fopen("nPostNeurons.dat", "rb");
  fread(sparseConVec, sizeof(*sparseConVec), N_NEURONS * (2 * K + 1), fpSparseConVec);
  fread(idxVec, sizeof(*idxVec), N_NEURONS, fpIdxVec);
  fread(nPostNeurons, sizeof(*nPostNeurons), N_NEURONS, fpNpostNeurons);
  fclose(fpSparseConVec);
  fclose(fpIdxVec);
  fclose(fpNpostNeurons);

  

}
*/
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

void GenSparseConMatDisp(sparseMat *sPtr[]) {
  int row, clm, count, kNeuron;
  printf("\n****** DISP FUNC *******\n");
  for(row = 1; row <= N_Neurons; ++row) {
    printf(" neurn %d is connected to :", row);
    for(kNeuron = 1; kNeuron <= sPtr[row]->nPostNeurons; ++kNeuron) { 
      printf(" %d", sPtr[row]->postNeuronIds[kNeuron]);
    }   
    printf("\n");
  }
}

/* read connection matrix into global variable conMat from text file */
void ReadConMatFromFile(const char* filename, int nNeurons) {
  FILE* fp;
  int i, j;
  char buffer[nNeurons * nNeurons];
  size_t result;
  
  /*  long fileSize;*/
  fp = fopen(filename, "r");

  /*  fseek(fp, 0, SEEK_END);
  fileSize = ftell(fp);
  rewind(fp);
  buffer = (char *)malloc(sizeof(char) * fileSize);
  if(buffer == NULL) {
    fputs("buffer not allocated", stderr);
    exit(1);
    }
  */
  if(fp != NULL) {
    result = fread(buffer, sizeof(int), nNeurons * nNeurons, fp);
    for(i = 0; i < nNeurons; ++i) {
      for(j = 0; j < nNeurons; ++j) {
        conMat[i+1][j+1] = (double)(buffer[i + j * nNeurons] - '0'); // convert ascii decimal character representation to integer 
        /*        printf("int %d\n", (buffer[i + j * nNeurons]) - '0');
                  printf("%f\n" , conMat[i+1][j+1]);*/
      }
    }
    fclose(fp);
  }
  else {
    printf("\n%s does not exist\n", filename);
  }
}
