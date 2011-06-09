#include "sp_type.h"
#include <cmath>

namespace T_sclpot
{

  double none::spot_( double x, double y, double z, double time, int me ) {
    return 0.0;
  }

  double coulomb::charge = 1.0;

  double coulomb::spot_( double x, double y, double z, double time, int me ) {
    return - charge / x;     
  }

  double Ar_1s::spot_( double x, double y, double z, double time, int me ) {
    return ( - 1.0 / x - 17.0 * exp( -17.43 * x ) / x );
  }
}
