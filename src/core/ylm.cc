#include<ylm.h>

complex ylm(long l, long m, double theta, double phi)
{
  double ylm;  
  double ylm0,ylm1;
  long lindex;

  if (m != 0) {fprintf(stderr, "err: Ylm: not implemented yet for m <> 0 \n"); exit(-71);}
  if (l<0) {fprintf(stderr, "err: Ylm: l<0\n"); exit(-69);}
  if (l<fabs(m)) {fprintf(stderr, "err: Ylm: l<|m|\n"); exit(-70);}
  

  if (l==0) {
    ylm=0.5/sqrt(M_PI);
  }
  else {
    if (l==1) {
      ylm=0.5*sqrt(3.0/M_PI)*cos(theta);
    }
    else
      {
	ylm0=0.5/sqrt(M_PI);
	ylm1=0.5*sqrt(3.0/M_PI)*cos(theta);
	lindex=2;
	do 
	  {
	    ylm=sqrt((2.0*lindex-1.0)*(2.0*lindex+1.0))/(1.0*lindex)*cos(theta)*ylm1 - (lindex-1.0)/(1.0*lindex)*sqrt((2.0*lindex+1.0)/(2.0*lindex-3.0))*ylm0;
	    ylm0=ylm1;
	    ylm1=ylm;
	    lindex++;
	  }
	while (lindex<=l);
      };
  };



  return complex(ylm,0.0);
}
