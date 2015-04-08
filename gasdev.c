
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include "nrutil.h"
#include "nr.h"
#include "globalVars.h"

double gasdev(long int *idum)

{
  double r1;
  double r2;
  const double r8_pi = 3.141592653589793;
  double x;
    
  r1 = ran1(idum);
  r2 = ran1(idum );

  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return gsl_ran_gaussian(gslRngState001, 1);
}

