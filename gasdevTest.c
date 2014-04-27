#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "nrutil.h"
#include "nr.h"
//  gcc ran1.c nrutil.c -g gaussrand.c  -lm

//double ga(long int *);

void main () {
  long seed, s2=-13;
  int i;
  FILE *fp;
  srand(time(NULL));
  fp =   fopen("../../results/randnout", "w");
  for(i=0; i <1000000; ++i) {
    seed = -1 * rand();
    fprintf(fp, "%f \n", gasdev(&seed));
  }
  fclose(fp);
}
