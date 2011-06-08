#include "wf_type.h"
#include "sphi.h"
#include <cstdio>

namespace T_wf
{

  double gaussian::R = 29.76;
  double gaussian::width = 3.17766;
  double gaussian::sigma = width;
  double gaussian::alpha = 0.5 / ( sigma * sigma );
  double gaussian::energy = 1.5 * alpha;
  
  complex gaussian::init( double r, int ell ) {
    double bessel = modified_sp_bessel_i( ell, 2.0 * alpha * R * r );
    double result = sqrt( 4.0 * M_PI * ( 2.0 * ell + 1.0 ) ) * 
      exp( - alpha * ( r*r + R*R ) ) * bessel * r;
    return result;
  }

  complex h1s::init( double r, int ell ) {
    complex result = complex( 0.0, 0.0 );
    if( ell == 0 ) {
      result = complex( r * exp( - r ), 0.0 );
    }
    return result;
  }
}
