#include "sphi.h"
#include <stdlib.h>

extern "C" { void sp_bess_i_( double* x, int* n, double* val ); }
extern "C" { void modified_sp0_( double* x, double* val ); }

double modified_sp_bessel_i( int n, double x )
{
  /*  
  if ( abs(x) < 1e-20 ) {
    if ( n == 0 ) {
      return 1.0;
    }
    else {
      return 0.0;
    }
  }
  */

  double x_, val;
  int n_;

  x_ = x;
  n_ = n;

  if ( n == 0 ) {
    modified_sp0_( &x, &val );
    return val;
  }

  sp_bess_i_( &x_, &n_, &val );
  return val;
}
