#ifndef _VARPROTOS_
#define _VARPROTOS_
extern double dt, *thetaVec;// integration time step 

extern double  *randnXiA, // norm rand number
  **randwZiA, // weibul rand number
  *randuDelta, // uniform rand (0, PI)
  **randuPhi, // uniform rand (0, 2.PI)
  *rTotalPrev, // rTotal(t - 1)
  *tempRandnPrev,
  *tempRandnNew,
  *Itgrl, *ItgrlOld;
extern char filebase[256];
extern FILE *isynapFP, *rTotalFP, *gbgrndFP, *gEEEIFP, *vmFP;


// recurrent input 
extern double *tempCurE, *tempCurI;

extern double **conMat, *iSynap, *expSum, *gEI_I, *gEI_E, *rTotal;

extern double contrast, *gFF, *iFF;

extern sparseMat *sConMat[];

extern int nTotSpks;

#endif
