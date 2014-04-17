#include <stdio.h>
#include <math.h>

void findPeaks(double a[], int lengthA, double peaks[], int *spkCount, double thresh) {

  double temp1, temp2;
  int loopIdx;
  // peaks = vector(1, lengthA);
  for(loopIdx = 1; loopIdx < lengthA - 1; ++loopIdx) {
    if(a[loopIdx] > thresh) { 
      temp1 = a[loopIdx - 1];
      temp2 = a[loopIdx +1];
      if(temp1 <= a[loopIdx] & temp2 <= a[loopIdx]) {
	(*spkCount) ++;
	//peaks[*spkCount] = loopIdx;
	peaks[loopIdx] = 1;
      }
    }
  }
}


/* void main() { */

/*   int  spkCount = 0; */
/*   double *peaks; */
/*   double v[10] = {0,1, 2, 10, 4, 5, 6, 12, 1, 2 }; */
/*   findPeaks(v, 10, peaks, &spkCount, 1); */
/*   printf("%d\n", spkCount); */
/* } */
