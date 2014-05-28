

__device__ void kernel(int nSpks, int *dev_spkNeuronId){
  //
  int nNeuron = threadIdx.x + blockDim.x * blockIdx.x;
  int i;
  /* compute squares of entries in data array */
  // !!!!! neurons ids start from ZERO  !!!!!! 
  if(nNeuron <= DEV_N_NEURONS) {
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
