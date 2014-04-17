/*** routine for the Runge-Kutta Method for the calculation of ODE 
 
     given values for y[1,...,n] and dxdy[1,...,n] known at x 

     returns the increment variables as yout[1,...,n] 

     user suplies routine "derives(x,y,dxdy)" which returns dxdy at x
***/


#define NRANSI
#include "nrutil.h"

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
	 void (*derivs)(double, double [], double []))
{
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;

	dym = vector(1,n);  /* define vector! */
	dyt = vector(1,n);
	yt = vector(1,n);
	hh = h*0.5;
	h6 = h/6.0;
	xh = x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i]; /* 1st step */
	(*derivs)(xh,yt,dyt);                     /* 2nd step */
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);                     /* 3rd step */
	for (i=1;i<=n;i++)
	{
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);                    /* 4th step */
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_vector(yt,1,n);
	free_vector(dyt,1,n);
	free_vector(dym,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 9,)5. */
