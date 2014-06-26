localLastNSteps = nstep - STORE_LAST_N_STEPS;
//*** TIMELOOP ***//
xx = x1;  
x = x1;
h = DT; //(x2 - x1) / nstep;
isynapNew = 0;
for (k = 0; k < nstep; k++) 
  {
    dev_IF_SPK[mNeuron] = 0;
    vmOld = v[0];
    derivs(x, v, dv, isynapNew);
    rk4(v, dv, N_STATEVARS, x, h, vout, isynapNew);
    //  ERROR HANDLE     if ((float)(x+h) == x) //nrerror("Step size too small in routine rkdumb");
    x += h; 
    xx = x; //xx[k+1] = x;
    /* RENAME */
    for (i = 0; i < N_STATEVARS; ++i) {
      v[i]=vout[i];
    }
    if(k > localLastNSteps) {
      y[mNeuron + N_NEURONS * k] = v[0];
      dev_isynap[mNeuron + N_NEURONS * k] = isynapNew;
    }
    if(k > 2) {
      if(v[0] > SPK_THRESH) { 
        if(vmOld <= SPK_THRESH) {
          dev_IF_SPK[mNeuron] = 1;
          atomicAdd(totNSpks, 1); // atomic add on global introduces memory latency
          localTotNspks = *totNSpks;
          spkNeuronId[localTotNspks] = mNeuron;
          spkTimes[localTotNspks] = xx;
        }
      }
    }
    __syncthreads(); // CRUTIAL step to ensure that dev_IF_spk is updated by all threads 
    isynapNew = isynap(v[0], dev_conVec);
