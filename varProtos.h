extern double  *randnXiA, // norm rand number
  **randwZiA, // weibul rand number
  *randuDelta, // uniform rand (0, PI)
  **randuPhi, // uniform rand (0, 2.PI)
  *rTotalPrev, // rTotal(t - 1)
  *tempRandnPrev,
  *tempRandnNew,
  *Itgrl, *ItgrlOld;
extern FILE *isynapFP, *rTotalFP, *gbgrndFP;


// recurrent input 
double *tempCurE, *tempCurI;

extern double **conMat, *iSynap, *expSum, *gEE, *gEI, *gIE, *gII, *rTotal;

extern double contrast, *gFF, *iFF;
