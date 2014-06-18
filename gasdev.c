/* #include <math.h> */
/* // generate randn */
/* float gasdev(long *idum) */
/* { */
/* 	float ran1(long *idum); */
/* 	static int iset=0; */
/* 	static float gset; */
/* 	float fac,rsq,v1,v2; */

/* 	if  (iset == 0) { */
/* 		do { */
/* 			v1=2.0*ran1(idum)-1.0; */
/* 			v2=2.0*ran1(idum)-1.0; */
/* 			rsq=v1*v1+v2*v2; */
/* 		} while (rsq >= 1.0 || rsq == 0.0); */
/* 		fac=sqrt(-2.0*log(rsq)/rsq); */
/* 		gset=v1*fac; */
/* 		iset=1; */
/* 		return v2*fac; */
/* 	} else { */
/* 		iset=0; */
/* 		return gset; */
/* 	} */
/* } */
/* /\* (C) Copr. 1986-92 Numerical Recipes Software 9,)5. *\/ */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "nrutil.h"
#include "nr.h"

/******************************************************************************/

double gasdev(long int *idum)

/******************************************************************************/
/*
  Purpose:

    R8_NORMAL_01 returns a unit pseudonormal R8.

  Discussion:

    The standard normal probability distribution function (PDF) has 
    mean 0 and standard deviation 1.

    Because this routine uses the Box Muller method, it requires pairs
    of uniform random values to generate a pair of normal random values.
    This means that on every other call, the code can use the second
    value that it calculated.

    However, if the user has changed the SEED value between calls,
    the routine automatically resets itself and discards the saved data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 August 2013

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8_NORMAL_01, a normally distributed random value.
*/
{
  double r1;
  double r2;
  const double r8_pi = 3.141592653589793;
  double x;
    
  r1 = ran2(idum);
  r2 = ran2(idum );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return x;
}

