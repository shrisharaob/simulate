#include<math.h>
#include<stdio.h>


void derive(double x, double stateVar[], double dydx[]) {
  int ddim = 2;
  dydx[1] = stateVar[2];
  dydx[2] =  -1 * stateVar[1]; 
 }

// T = 2 * PI
