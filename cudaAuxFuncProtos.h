

void cudakernel(int nSpks, int *dev_spkNeuronId, float *dev_conVec, float *gE_data, float *gI_data);
void CudaInitISynap(float *conVec);
void CudaISynap(int nSpks, int *spkNeuronId);
void CudaFreeMem();
