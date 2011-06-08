#include<factorial.h>

double factorial(long n)
{
  double result=1.0;
  long l;
  static double res[172] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};
  static long factorial_no_res = 10;

  if (n<0) {fprintf(stderr, "err: factorial: n<0\n"); exit(-71);}
  if (n>170) {fprintf(stderr, "err: double factorial(): n>170\n"); exit(-72);}  
  if (n<factorial_no_res) return(res[n]);
  
  for (l=factorial_no_res; l<=n; l++) res[l]=res[l-1]*l;

  factorial_no_res = l;
  return res[factorial_no_res-1];
}
