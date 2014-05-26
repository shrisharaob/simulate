extern double dt, *thetaVec;// integration time step 

extern double  *randnXiA, // norm rand number
  **randwZiA, // weibul rand number
  *randuDelta, // uniform rand (0, PI)
  **randuPhi, // uniform rand (0, 2.PI)
  *rTotalPrev, // rTotal(t - 1)
  *tempRandnPrev,
  *tempRandnNew,
  *Itgrl, *ItgrlOld;
extern FILE *isynapFP, *rTotalFP, *gbgrndFP, *gEEEIFP, *vmFP;


// recurrent input 
extern double *tempCurE, *tempCurI;

extern double **conMat, *iSynap, *expSum, *rTotal;
extern float *gEI_I, *gEI_E,  *randList;
extern double contrast, *gFF, *iFF;

extern sparseMat *sConMat[];
extern float conVec[];
extern int randListCounter;
extern double K; 
