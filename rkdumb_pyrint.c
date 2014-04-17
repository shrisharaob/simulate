/* driver for the algorithm "rk4.c" */



/* ************************************************************ */
/* ************************************************************ */
/*        This routine DOES NOT save EVERY TIMESTEP!!!          */  
/* ************************************************************ */
/* ************************************************************ */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nr.h"
#include "xpyrint_var.h"
#include "xpyrint_func.h" 

void rkdumb_pyrint(double vstart[], int nvar, double x1, double x2, int nstep,
	    void (*derivs)(double, double [], double []))
{
  void rk4(double y[], double dydx[], int n,  double x, double h, 
	   double yout[],
	   void (*derivs)(double, double [], double []));
  
  
  /* ***************************************************** */
  /* ***************************************************** */
  /* ****************** VARIABLES ************************ */

  long idum = -13;

  int i, j, k, s, l, m=0;

  double x, h;
  double *v, *a, *b, *as, *bs, *vout, *dv;

  double **randpp, **randii, **randip, **randpi;
  int **delaypp, **delayii, **delayip, **delaypi;
  
  double distance, sigp = 4, sigi = 9, da, db;

  int np=0, ni=0; // ??? def
  int *peakp_bin, *peaki_bin;

  double dum, *xi, *noisei, *xp, *noisep; 
  double currentpp=0.0, currentpi=0.0, currentip=0.0, currentii=0.0;

  double frequencyp, frequencyi;

  double ***peakpp, ***peakpi, ***peakip, ***peakii;

  char count[20], filev[10], filep[10], filey[10], filetest[10], filesynmx[10], filesyndy[10]; 
  FILE *outvolt, *outpeak, *outsyncurr, *outtest, *outsynmx, *outsyndy;

  
 //*** NAME OUTFILES ***//
  sprintf(count,"%d",loop);
  
  strcpy(filev,"z_volt_");
  strcat(filev, count);
  
  strcpy(filep,"z_peak_");
  strcat(filep, count);
  
  strcpy(filey,"z_syncurr_");
  strcat(filey, count);

  strcpy(filetest,"z_test_");
  strcat(filetest, count);

  strcpy(filesynmx,"z_synmatrix_");
  strcat(filesynmx, count);

  strcpy(filesyndy,"z_syndelay_");
  strcat(filesyndy, count);


  //*** OPEN FILES ***//
  outvolt = fopen(filev,"w");
  outsyncurr = fopen(filey,"w");
  outpeak = fopen(filep,"w");
  outtest = fopen(filetest,"w");
  outsynmx = fopen(filesynmx,"w");
  outsyndy = fopen(filesyndy,"w");
  

  v = vector(1,nvar);
  a = vector(1,nvar); // ??? def
  b = vector(1,nvar); // ???
  as = vector(1,nvar); // ???
  bs = vector(1,nvar); // ???
  vout = vector(1,nvar);
  dv = vector(1,nvar);

  pint = vector(1,NEURON);
  iint = vector(1,INTER);

  peakp_bin = ivector(1,BIN);
  peaki_bin = ivector(1,BIN);

  noisep = vector(1,NEURON);
  noisei = vector(1,INTER);
  
  xp = vector(1,NEURON); 
  xi = vector(1,INTER); 

  sp = vector(1,NEURON); 
  si = vector(1,INTER); 

  randpp = matrix(1,NEURON,1,NEURON);
  randpi = matrix(1,NEURON,1,INTER);
  randii = matrix(1,INTER,1,INTER);
  randip = matrix(1,INTER,1,NEURON);
  delaypp = imatrix(1,NEURON,1,NEURON);
  delaypi = imatrix(1,NEURON,1,INTER);
  delayii = imatrix(1,INTER,1,INTER);
  delayip = imatrix(1,INTER,1,NEURON);

  synpp = vector(1,NEURON);
  synpi = vector(1,INTER);
  synip = vector(1,NEURON);
  synii = vector(1,INTER);

  peakpp = f3tensor(1,NEURON,1,NEURON,1,DELAYMAX);
  peakpi = f3tensor(1,NEURON,1,INTER,1,DELAYMAX);
  peakip = f3tensor(1,INTER,1,NEURON,1,DELAYMAX);
  peakii = f3tensor(1,INTER,1,INTER,1,DELAYMAX);

  //*** START ***//
  for (i=1;i<=nvar;i++) 
    {
      v[i]=vstart[i];
      a[i]=vstart[i];
      as[i]=vstart[i];
    }
  
  //*** INITIATE ***//
  for(i=1;i<=NEURON;++i)
    {
      noisep[i]=0.0;
      xp[i]=0.0; 
      pint[i] = 0.0;
      synpp[i] = 0.0;
      synip[i] = 0.0;
    }
  for(i=1;i<=INTER;++i)
    {
      noisei[i]=0.0; 
      xi[i]=0.0;
      iint[i] = 0.0;
      synii[i] = 0.0;
      synpi[i] = 0.0;
    }
  

  //for (i=1; i<=20; ++i)
  // {
  //   k= ((i-5+10)%10);
  //  printf("%d ",k);
  //}

  //*** CONNECTIVITY ***//
  // 2 D array

  // 2D matrix: 
  np=0;
  for (i=1; i<=GRIDP;++i)
    {
      for (j=1; j<=GRIDP;++j)
	{
	  np = np+1; 
	  //PP
	  ni = 0;
	  for (k=1; k<=GRIDP;++k)
	    {
	      for (l=1; l<=GRIDP;++l)
		{
		  ni = ni+1;
		  if (abs(k-i)<GRIDP/2) da = pow(k-i,2);
		  else da = pow(GRIDP-abs(k-i),2);
		  if (abs(l-j)<GRIDP/2) db = pow(l-j,2);
		  else db = pow(GRIDP-abs(l-j),2);
		  distance =  da + db;
		  randpp[np][ni] = exp(-distance/sigp)*rand()/(RAND_MAX+1.0);
		  fprintf(outsynmx,"%d %d %f\n",np,ni,randpp[np][ni]);
		}
	    }
	  //PI
	  ni = 0;
	  for (k=1; k<=GRIDI;++k)
	    {
	      for (l=1; l<=GRIDI;++l)
		{
		  ni = ni+1;
		  if (abs(k*2-0.5-i)<GRIDP/2) da = pow(k*2-0.5-i,2);
		  else da = pow(GRIDP-abs(k*2-0.5-i),2);
		  if (abs(l*2-0.5-j)<GRIDP/2) db = pow(l*2-0.5-j,2);
		  else db = pow(GRIDP-abs(l*2-0.5-j),2);
		  distance =  da + db;
		  //distance = pow((k*2-i)%(GRIDP)+0.5,2) + pow((l*2-j)%(GRIDP)+0.5,2);
		  randpi[np][ni] = exp(-distance/sigi)*rand()/(RAND_MAX+1.0);
		  fprintf(outsynmx,"%d %d %f\n",np,ni+NEURON,randpi[np][ni]);
		}
	    }
	}
    }
 np=0;
  for (i=1; i<=GRIDI;++i)
    {
      for (j=1; j<=GRIDI;++j)
	{
	  np = np+1; 
	  //IP
	  ni = 0;
	  for (k=1; k<=GRIDP;++k)
	    {
	      for (l=1; l<=GRIDP;++l)
		{
		  ni = ni+1;
		  if (abs(k-(i*2-0.5))<GRIDP/2) da = pow(k-(i*2-0.5),2);
		  else da = pow(GRIDP-abs(k-(i*2-0.5)),2);
		  if (abs(l-(j*2-0.5))<GRIDP/2) db = pow(l-(j*2-0.5),2);
		  else db = pow(GRIDP-abs(l-(j*2-0.5)),2);
		  distance =  da + db;
		  //distance = pow((k-i*2)%(GRIDP)+0.5,2) + pow((l-j*2)%(GRIDP)+0.5,2);
		  randip[np][ni] = exp(-distance/sigi)*rand()/(RAND_MAX+1.0);
		  fprintf(outsynmx,"%d %d %f\n",np+NEURON,ni,randip[np][ni]);
		}
	    }
	  //II
	  ni = 0;
	  for (k=1; k<=GRIDI;++k)
	    {
	      for (l=1; l<=GRIDI;++l)
		{
		  ni = ni+1;
		  if (abs(k-i)*2<GRIDP/2) da = pow((k-i)*2,2);
		  else da = pow(GRIDP-abs(k-i)*2,2);
		  if (abs(l-j)*2<GRIDP/2) db = pow((l-j)*2,2);
		  else db = pow(GRIDP-abs(l-j)*2,2);
		  distance =  da + db;
		  //distance = pow((k*2-i*2)%(GRIDP),2) + pow((l*2-j*2)%(GRIDP),2);
		  randii[np][ni] = exp(-distance/sigi)*rand()/(RAND_MAX+1.0);
		  fprintf(outsynmx,"%d %d %f\n",np+NEURON,ni+NEURON,randii[np][ni]);
		}
	    }
	}
    }
  

  // synaptic delay
  for (i=1; i<=NEURON; ++i) 
    { 
      /* pyramidal -> pyramidal */      
      for (k=1;k<=NEURON; ++k)
	  delaypp[i][k] = 50;//((int) (1.5*NSTEP/X2*rand()/(RAND_MAX+1.0)));
		    	
      /* pyramidal -> interneuron */	  
      for (k=1;k<=INTER; ++k)
	delaypi[i][k] = 50;//((int) (1.5*NSTEP/X2*rand()/(RAND_MAX+1.0)));
    }

  for (i=1; i<=INTER; i++)
    {
      /* interneuron -> interneuron */
      for (j=1;j<=INTER; j++)
	  delayii[i][j] = 50;//((int) (1.5*NSTEP/X2*rand()/(RAND_MAX+1.0)));

      /* interneuron -> pyramidal */
      for (j=1;j<=NEURON; j++)
	  delayip[i][j] = 50;//((int) (1.5*NSTEP/X2*rand()/(RAND_MAX+1.0)));
    } 
  
  // *** PRINT FRIRST VALUE *** //
  //fprintf(outvolt,"%f %f %f\n",x1,vstart[1],vstart[1+IVAR]);


  //*** TIMELOOP ***//
  xx[1] = x1;  x = x1;
  h = (x2-x1)/nstep;
  
  for (k=1;k<=nstep;k++) 
    {
      
      (*derivs)(x,v,dv);
      
      rk4(v,dv,nvar,x,h,vout,derivs);
      
         
      if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
      
      x += h;
      xx[k+1]=x;
      

      /* RENAME */  
      for (i=1;i<=nvar;i++) 
	{
	  v[i]=vout[i];
	  if(k==1){ a[i]=vstart[i]; b[i]=0.0;}
	}

      
      // PRINT //
      if(~k%(50))
	fprintf(outvolt,"%f %f %f %f %f\n",
		xx[k+1],v[1],v[PVAR+1],v[NEUVAR+1],v[NEUVAR+IVAR+1]);
      
      if(~k%(100)==0) printf(" %f ",x);
      
      // *** NEURONS *** //
      currentpp = 0.0;
      currentip = 0.0;
      for(i=1;i<=NEURON;++i)
        {
          j=(i-1)*PVAR;

	  synpp[i] = synpp[i]*exp(-h/taupx);
	  synip[i] = synip[i]*exp(-h/taupx);
	  
	   //PEAK// find peaks
	  np = 0;
	  if(k>3 & a[j+1]<=b[j+1] & b[j+1]>=v[j+1] & b[j+1]>=1.0)
	    {
	      //PRINT SPIKES   // ??? spike gen code
	      fprintf(outpeak,"%f %d\n",xx[k],i);
	      np = 1;
	    }
		
	  //UPDATE//
	  if(k>=2)
            {
	      a[j+1] = b[j+1];// find peaks
	      b[j+1] = v[j+1];
	    }
	  
	  //SYNAPTIC UPDATE// 
	  for(l=1;l<=NEURON;++l)
	    {//implementation of synaptic delay
	      peakpp[i][l][(k+delaypp[i][l])%(DELAYMAX)] = np;
	      synpp[l] = synpp[l] + randpp[i][l]*peakpp[i][l][(k)%(DELAYMAX)];
	      //fprintf(outtest,"%f ",randpp[i][l]*peakpp[i][l][(k)%(DELAYMAX)]);
	      peakpp[i][l][(k)%(DELAYMAX)]=0;
	    }
	  for(l=1;l<=INTER;++l)
	    {
	      peakpi[i][l][(k+delaypi[i][l])%(DELAYMAX)] = np;
	      synpi[l] = synpi[l] + randpi[i][l]*peakpi[i][l][(k)%(DELAYMAX)];
	      peakpi[i][l][(k)%(DELAYMAX)]=0;
	    }
	  
          //NOISE//
          xp[i] = xp[i]*exp(-h/taupx);
          while(noisep[i]<=xx[k+1])
            {
              xp[i] = xp[i] + exp(-(xx[k+1]-noisep[i])/taupx);
              
              do
                dum=ran1(&idum);
              while (dum==0.0);
              
              noisep[i] = noisep[i] - log(dum)/pnoise; 
            }
          sp[i] = xp[i];
	  //CALLIBRATE SYNAPSE:	  sp[2]=0;

	  //EXTERNAL CURRENT//
	  pint[i] = 0.0;
	  //if(xx[k+1]>=500.0 & xx[k+1]<=700.0)
	  //pint[i] = IPYR*exp(-(xx[k+1]-NSTEP/2)*(xx[k+1]-NSTEP/2));
	  //pint[i] = IPYR*(cos(2*DEGPI*8*X2/NSTEP*k/1000))*exp(-(k-NSTEP/2)*(k-NSTEP/2)/(NSTEP*NSTEP/16));

	  //AVERAGE CURENT - P//
	  currentpp += Isynpp(v[1+j],v[7+j])/NEURON;
	  currentip += Isynip(v[1+j],v[8+j])/NEURON;
	}
      fprintf(outtest,"%f %f\n",xx[k+1],-Isynpp(v[1],v[7]));
	    
      //INTER//
      currentpi = 0.0;
      currentii = 0.0;
      for(i=1;i<=INTER;++i)
        {
          j=(i-1)*IVAR+NEUVAR;

	  synpi[i] = synpi[i]*exp(-h/taupx);
	  synii[i] = synii[i]*exp(-h/taupx);
	  
	  //PEAK//
	  ni=0;
	  if(k>3 & a[j+1]<=b[j+1] & b[j+1]>=v[j+1] & b[j+1]>=0.0)
	    {
	      ni=1;
	      //PRINT SPIKES
	      fprintf(outpeak,"%f %d\n",xx[k],i+NEURON);
	    }
	  
          //UPDATE/
          if(k>=2)
            {
	      a[j+1] = b[j+1];
              b[j+1] = v[j+1];
            }
	  
	  //SYNAPTIC UPDATE//
	  for(l=1;l<=NEURON;++l)
	    {
	      peakip[i][l][(k+delayip[i][l])%(DELAYMAX)] = ni;
	      synip[l] = synip[l] + randip[i][l]*peakip[i][l][(k)%(DELAYMAX)];
	      peakip[i][l][(k)%(DELAYMAX)]=0;
	    }
	  for(l=1;l<=INTER;++l)
	    {
	      peakii[i][l][(k+delayii[i][l])%(DELAYMAX)] = ni;
	      synii[l] = synii[l] + randii[i][l]*peakii[i][l][(k)%(DELAYMAX)];
	      peakii[i][l][(k)%(DELAYMAX)]=0;
	    }

	  // NOISE //
	  xi[i] = xi[i]*exp(-h/taupx);
	  while(noisei[i]<=xx[k+1])
            {
              xi[i] = xi[i] + exp(-(xx[k+1]-noisei[i])/taupx);
              
              do
                dum=ran1(&idum);
              while (dum==0.0);
              
              noisei[i] = noisei[i] - log(dum)/inoise;
            }
          si[i] = xi[i];
	  //CALLIBRATE SYNAPSE:	  si[2]=0;

	  //EXTERNAL CURRENT
	  pint[i] = 0.0;
	  //if(xx[k+1]>=500.0 & xx[k+1]<=700.0)
	  //iint[i] = IINT*exp(-(xx[k+1]-i*X2/INTER)*(xx[k+1]-i*X2/INTER));
	  //iint[i] = IINT*(cos(2*DEGPI*8*X2/NSTEP*k/1000))*exp(-(k-NSTEP/2)*(k-NSTEP/2)/(NSTEP*NSTEP/16));

	  //AVERAGE CURENT - I//
	  currentpi += Isynpi(v[1+j],v[5+j])/INTER;
	  currentii += Isynii(v[1+j],v[6+j])/INTER; 
	}

      fprintf(outsyncurr,"%f %f %f %f %f\n",xx[k+1],currentpp,currentip,currentpi,currentii);

    } // TIME LOOP CLOSE
  
  
  //FREQUENCY//
  //frequencyp = 1000.0*np/NEURON/x2;
  //frequencyi = 1000.0*ni/INTER/x2;
  //printf("%d, %d\n",np,ni);
  
  //PRINT//
  //fprintf(outfreq,"%f %f %f %f %f %f %f %f %f %f %f\n",pnoise,inoise,
  //	  frequencyp,frequp_max,kappap,
  //	  frequencyi,frequi_max,kappai,
  //	  -currentpp/currentip,-currentpi/currentii,-currentpp/currentii);
  //printf("(%d) PYR:f=%3.2f fpop=%3.2f cv=%3.2f\n",
  //	 loop,frequencyp,frequp_max,kappap); 
  //printf("    INT:f=%3.2f fpop=%3.2f cv=%3.2f\n",
  //	 frequencyi,frequi_max,kappai); 
  
    
  
  // CLOSE //
  fclose(outvolt);
  fclose(outsyncurr);
  fclose(outpeak);
  fclose(outtest);

  //free_vector(v,1,nvar);
  //free_vector(a,1,nvar);
  //free_vector(b,1,nvar);
  //free_vector(as,1,nvar);
  //free_vector(bs,1,nvar);
  //free_vector(vout,1,nvar);
  //free_vector(dv,1,nvar);
  
  //free_matrix(randpp,1,NEURON,1,NEURON);
  //free_matrix(randpi,1,NEURON,1,INTER);
  //free_matrix(randii,1,INTER,1,NEURON);
  //free_matrix(randip,1,INTER,1,INTER);
  
  //free_vector(noisep,1,NEURON);
  //free_vector(noisei,1,INTER);
  
  //free_vector(xp,1,NEURON);
  //free_vector(xi,1,INTER);
  
  //free_vector(sp,1,NEURON);
  //free_vector(si,1,INTER);
  
  //free_f3tensor(peakpp,1,NEURON,1,NEURON,1,DELAYMAX);
  //free_f3tensor(peakip,1,INTER,1,NEURON,1,DELAYMAX);
  //free_f3tensor(peakpi,1,NEURON,1,INTER,1,DELAYMAX);
  //ree_f3tensor(peakii,1,INTER,1,INTER,1,DELAYMAX);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 9,)5. */
