
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "nrutil.h"
#include "nr.h"

/******************************************************************************/

double gaussrand( long int *seed )

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

  r1 = ran1(seed);
  r2 = ran1(seed );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return x;
}


