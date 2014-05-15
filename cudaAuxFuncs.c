#include <stdio.h>
#include <math.h>
#include <time.h>
#include "cuda.h"
#include <curand.h>
#include "cuda_runtime_api.h"
#include "devHostConstants.h"
#include "cudaAuxFuncProtos.h"

extern float *dev_gEI_E, *dev_gEI_I, *dev_conVec, *gEI_E, *gEI_I, *randList;
extern int randListCounter;
extern int *dev_spkNeuronId;

void MyCudaErrorHandler(cudaError_t err);

void MyCudaErrorHandler(cudaError_t err) {
  if (err != cudaSuccess)
    {
      fprintf(stderr, "ERROR = %s\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
    }
}
void CudaInitISynap(float *conVec) {
  int n;
  float *tmp;
  dev_gEI_E = NULL;
  dev_gEI_E = NULL;
  dev_conVec = NULL;

  MyCudaErrorHandler(cudaMallocHost((void **)&gEI_I, (N_NEURONS + 1) * sizeof(float))); 
  MyCudaErrorHandler(cudaMallocHost((void **)&gEI_E, (N_NEURONS + 1) * sizeof(float))); 
  MyCudaErrorHandler(cudaMalloc((void **)&dev_gEI_E, (N_NEURONS + 1) * sizeof(float)));
  MyCudaErrorHandler(cudaMalloc((void **)&dev_gEI_I, (N_NEURONS + 1) * sizeof(float)));
    /* fill up the host array */
  for(n=0; n <= N_NEURONS ;++n) {
    gEI_E[n] = 0.0;
    gEI_I[n] = 0.0;
  }
  // +1 to make it compatible with the stupid numerical recipies indexing from 1
 MyCudaErrorHandler(cudaMemcpy(dev_gEI_E, gEI_E, (N_NEURONS + 1) * sizeof(float), cudaMemcpyHostToDevice));
 MyCudaErrorHandler(cudaMemcpy(dev_gEI_I, gEI_I, (N_NEURONS + 1) * sizeof(float), cudaMemcpyHostToDevice));
 MyCudaErrorHandler(cudaMalloc((void **)&dev_conVec, (N_NEURONS + 1) * (N_NEURONS + 1) * sizeof(float)));
 MyCudaErrorHandler(cudaMemcpy(dev_conVec, conVec, (N_NEURONS + 1) * (N_NEURONS + 1) * sizeof(float), cudaMemcpyHostToDevice));

  /* printf(" ... host and dev ptrs on aux ... \n"); */
  /* printf("%p \n", gEI_E); */
  /* printf("%p \n", gEI_I); */
  /* printf("%p \n", dev_gEI_E); */
  /* printf("%p \n", dev_gEI_I); */
  /* printf("%p \n", dev_conVec); */

 //  printf("\n convec ");
  MyCudaErrorHandler(cudaMallocHost((void **)&tmp, (N_NEURONS + 1) * (N_NEURONS + 1) * sizeof(float))); 
 MyCudaErrorHandler(cudaMemcpy(tmp, dev_conVec, (N_NEURONS + 1) * (N_NEURONS + 1) * sizeof(float), cudaMemcpyDeviceToHost));
 // MyCudaErrorHandler(cudaMalloc((void **)&dev_spkNeuronId, N_NEURONS * sizeof(int)));
 // printf("\n convec 0x on cuda aux: %p \n", conVec);
  /*  for(n = 0; n < (N_NEURONS + 1) * (N_NEURONS + 1); ++n) { */
  /*   printf("%d ", (int)tmp[n]); */
  /* } */
  /* printf("\n"); */
}
void CudaISynap(int nSpks, int *spkNeuronId) {
  int loopIdx;
  //  FILE *fp;
   if(nSpks > 0) {  // optimize - malloc called on every timestep !!!! latency ???
     MyCudaErrorHandler(cudaMalloc((void **)&dev_spkNeuronId, nSpks * sizeof(float)));
     MyCudaErrorHandler(cudaMemcpy(dev_spkNeuronId, spkNeuronId, nSpks * sizeof(int), cudaMemcpyHostToDevice));
  }
  else {
    dev_spkNeuronId = NULL; 
  }
  cudakernel(nSpks, dev_spkNeuronId, dev_conVec, dev_gEI_E, dev_gEI_I);
  // copy from device array to host array 
  MyCudaErrorHandler(cudaMemcpy(gEI_E, dev_gEI_E, (N_NEURONS + 1) * sizeof(float), cudaMemcpyDeviceToHost));
  MyCudaErrorHandler(cudaMemcpy(gEI_I, dev_gEI_I, (N_NEURONS + 1) * sizeof(float), cudaMemcpyDeviceToHost));
  MyCudaErrorHandler(cudaGetLastError());
  
  /* if(nSpks >0) {  */
  /*   // printf("address spkneuronid on cuda aux =  %p \n", spkNeuronId); */
  /*   printf("\n  SPIKE !!!! \n"); */
  /*   for(loopIdx=0; loopIdx < nSpks; ++loopIdx) { */
  /*     printf("nSpks = %d\n", nSpks); */
  /*     printf("spkneuron ids %d , idx = %d \n", spkNeuronId[loopIdx], loopIdx); */
  /*     printf("%f %f %f", gEI_E[0], gEI_E[1], gEI_E[2]); */
  /*   } */
  /*   printf("\n"); */
  /* } */
}

void CudaFreeMem() {
  cudaError_t  err = cudaSuccess;
  cudaFreeHost(gEI_E);
  cudaFreeHost(gEI_I);
  cudaFree(dev_gEI_E);
  cudaFree(dev_gEI_I);
  cudaFree(dev_conVec);
  //  err = cudaDeviceReset(); necessary ???
  err = cudaGetLastError();
  if (err != cudaSuccess)
    {
      printf("\nFailed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
    }
}

// uniform random number generator
void CudaRandGen(size_t n) {
  curandGenerator_t gen;  
  float *dev_vec;
  //  randList = (float *)calloc(n, sizeof(float));
  /* Allocate n floats on device */
  cudaMalloc((void **)&dev_vec, n * sizeof(float));
  /* Create pseudo-random number generator */
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  /* Set seed */
  curandSetPseudoRandomGeneratorSeed(gen, time(NULL));
  /* Generate n floats on device */
  curandGenerateUniform(gen, dev_vec, n);
  /* Copy device memory to host */
  cudaMemcpy(randList, dev_vec, n * sizeof(float), cudaMemcpyDeviceToHost);
  curandDestroyGenerator(gen); 
  cudaFree(dev_vec);
}

void CudaAccessURandList(int n, float **randVec) {
  int row = 0;
  int ii = 0;
  float *tmpVec;
  printf("Inside cudagenrandlist cntr = %d, randList ptr = %p \n", randListCounter, randList);
  if(randListCounter >= 0 && randListCounter < MAX_UNI_RANDOM_VEC_LENGTH) {
    //    printf("Accessing randList at %p \n", randList);
    //    printf("ranVec at %p of length %d \n", *randVec, n);
    *randVec= (float *)malloc(n * sizeof(float));
    tmpVec = *randVec;
    //    printf("created ranVec at %p of length %d, tmpVec = %p \n", *randVec, n, tmpVec);
    //    printf("ranVec at %f\n ", tmpVec[10]);
    //    pause(5000);
    //printf("\n Accessing randList at %p  \n", randList);
    for(ii = 0; ii < n; ++ii) {
      tmpVec[ii] = randList[randListCounter];
      //printf("%f ", tmpVec[ii]);
      randListCounter += 1;
      if(randListCounter >= MAX_UNI_RANDOM_VEC_LENGTH) {
      	randListCounter = -1;
      	CudaGenURandList();
      }
    }
  }
  //  printf("\n done with access ! \n");
}


void CudaGenURandList() {
  printf("Malloc randList at %p \n", randList);
  printf("Generating list of %d rand floats on GPU...", (int)MAX_UNI_RANDOM_VEC_LENGTH);
  fflush(stdout);
  CudaRandGen((size_t) MAX_UNI_RANDOM_VEC_LENGTH);
  randListCounter = 0;
  printf("done ! \n");
    //    printf("ranVe at %p %p \n", *randVec, randList);
    //pause(5000);
  /*   if(*randVec != randList) { */
  /*     *randVec = (float *)malloc(n * sizeof(float)); */
  /*     //printf("\n Accessing randList at %p  \n", randList); */
  /*     for(i = 0; i < n; ++i) { */
  /* 	*randVec[i] = randList[randListCounter]; */
  /* 	randListCounter += 1; */
  /*     } */
  /*     printf("returned randVec at %p \n", randVec); */
  /*   } */
  /*   //    else{printf("did not malloc randVec !!! \n");} */
  /* } */
}

float CudaURand() {
  //  float *tmpVec;
  //tmpVec = NULL;
  //  CudaAccessURandList(1, &tmpVec);
  if(randListCounter >= 0) {
    randListCounter += 1;
    if(randListCounter < MAX_UNI_RANDOM_VEC_LENGTH) {
      return randList[randListCounter];
    }
  }
  else {
    randListCounter = -1;
    CudaGenURandList();
    randListCounter += 1;
    return randList[0];
  }
}
