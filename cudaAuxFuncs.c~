#include <stdio.h>
#include <math.h>
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "devHostConstants.h"
#include "cudaAuxFuncProtos.h"

extern float *dev_gEI_E, *dev_gEI_I, *dev_conVec, *gEI_E, *gEI_I;
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
  //  printf("here 0---> \n");
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

  printf(" ... host and dev ptrs on aux ... \n");
  printf("%p \n", gEI_E);
  printf("%p \n", gEI_I);
  printf("%p \n", dev_gEI_E);
  printf("%p \n", dev_gEI_I);
  printf("%p \n", dev_conVec);

  printf("\n convec ");
  MyCudaErrorHandler(cudaMallocHost((void **)&tmp, (N_NEURONS + 1) * (N_NEURONS + 1) * sizeof(float))); 
 MyCudaErrorHandler(cudaMemcpy(tmp, dev_conVec, (N_NEURONS + 1) * (N_NEURONS + 1) * sizeof(float), cudaMemcpyDeviceToHost));
  for(n = 0; n<16; ++n) {
    printf("%d ", (int)tmp[n]);
  }
  printf("\n");
}
void CudaISynap(int nSpks, int *spkNeuronId) {
  int loopIdx;
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

  for(loopIdx=1; loopIdx <= N_NEURONS; ++loopIdx) {
    //    printf("nSpks = %d\n", nSpks);
      printf("%f ", gEI_I[loopIdx]);
    }
  printf("\n");
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

