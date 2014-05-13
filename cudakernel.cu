#include "cuda.h"
#include <math.h>
#include "devHostConstants.h"
//#include "devGlobalVars.h"
__constant__ float DEV_EXP = EXP_SUM;
__constant__ int DEV_N_NEURONS = N_NEURONS, DEV_NE = NE;

/* gpu kernel */
__global__ void kernel(int nSpks, int *dev_spkNeuronId, float *dev_conVec, 
		       float *gE_data, float *gI_data){
  //
  int nNeuron = blockIdx.x;// + blockDim.x * blockIdx.x;
  int i;
  /* compute squares of entries in data array */
  // !!!!! neurons ids start from ZERO  !!!!!! 
  if(nNeuron <= DEV_N_NEURONS) {
    //g_result[nNeuron] = g_data[nNeuron];
    gE_data[nNeuron] *= DEV_EXP;
    gI_data[nNeuron] *= DEV_EXP;
    if(nSpks > 0){
      for(i = 0; i < nSpks; ++i) { //  
	if(dev_spkNeuronId[i] <= DEV_NE) {
	  gE_data[nNeuron] += dev_conVec[dev_spkNeuronId[i] * (DEV_N_NEURONS + 1) + nNeuron]; //optimize !!!! gEI_E
	}
	else {
	  gI_data[nNeuron] += dev_conVec[dev_spkNeuronId[i] * (DEV_N_NEURONS + 1) + nNeuron]; //optimize !!!! gEI_I
	}
      }
    }
  }
}
/* only use extern if calling code is C */
extern "C" 
{
  /* driver for kernel */
  void cudakernel(int nSpks, int *dev_spkNeuronId, float *dev_conVec,
		  float *gE_data, float *gI_data){
    /* choose 256 threads per block for high occupancy */
    //int ThreadsPerBlock = 256;
      /* find number of blocks */
    //    int BlocksPerGrid = (DEV_N_NEURONS+ThreadsPerBlock-1)/ThreadsPerBlock;
      /* invoke device on this block/thread grid */
    // kernel <<< BlocksPerGrid, ThreadsPerBlock >>> (DEV_N_NEURONS, nSpks, dev_spkNeuronId, dev_conVec, 
    // 						   g_data, g_result);
    kernel <<< 8, 1>>> (nSpks, dev_spkNeuronId, dev_conVec, gE_data, gI_data);
  }
}
