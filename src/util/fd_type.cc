#include "fd_type.h"
#include <cmath>

namespace T_field
{

  double none::field_( double time, int me )
  {
    return 0.0;
  }

  double sin2_A::w = 0.0227817;
  double sin2_A::E_0 = -0.0168;
  double sin2_A::n_c = 0.5;
    
  double sin2_A::field_( double time, int me )
  {
    double wt = w * time;
    double a = 0.5 * wt / n_c;
    double result = E_0 * sin(a) * ( cos(a) * sin(wt) + n_c * cos(wt) * sin(a) ) / n_c;
    return result;
  }
}

namespace fd_type
{
  double field(double time, int me)
  {
    double result = 0.0;
    return result;
  }  
}
