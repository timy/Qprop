#include "sp_type.h"
#include <cmath>

/*
namespace sp_type
{
  namespace Coulomb
  {
    double charge = 1.0;

    void set_charge( double charge_ )
    {
      charge = charge_;
      return;
    }

    double spot( double x, double y, double z, double time, int me )
    {
      double result = - charge / x;
      return result;
    }
  }

  namespace helium_atom
  {
    double spot( double x, double y, double z, double time, int me )
    {
      return -17.0 * exp( -17.43 * x ) / x;
    }
  }

  namespace none
  {
    double spot( double x, double y, double z, double time, int me )
    {
      return 0.0;
    }
  }

}

*/

namespace T_sclpot
{

  double none::spot_( double x, double y, double z, double time, int me ) {
    return 0.0;
  }

  double coulomb::charge = 1.0;

  double coulomb::spot_( double x, double y, double z, double time, int me ) {
    return - charge / x;     
  }

  double helium::spot_( double x, double y, double z, double time, int me ) {
    return - 17.0 * exp( -17.43 * x ) / x;
  }
}
