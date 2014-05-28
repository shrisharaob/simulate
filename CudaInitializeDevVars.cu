__device__ double y[N_NEURONS * 2];
__device__ 


global__ void CudaInitializeDevVars() {
  double *vstart, *spkTimes;
  
  dt = DT;
  nSteps = (int)((x2 - x1) / dt);
  xx = (1, nSteps);
  y = matrix(1, N_StateVars * N_Neurons, 1, nSteps);
  input_cur = vector(1, nSteps);
  vstart = vector(1, N_StateVars * N_Neurons);
  IF_SPK = vector(1, N_Neurons);   
  expSum = vector(1, N_Neurons); // synap input
  iSynap = vector(1, N_Neurons); 
  //    gEI_I = vector(1, N_Neurons); // E --> E & E --> I
  //    gEI_E = vector(1, N_Neurons); // I --> E & I --> I  
  spkTimes = vector(1, nSteps);
  conMat = matrix(1, N_Neurons, 1, N_Neurons);
        

}
