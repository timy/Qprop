#include "vp_type.h"
/*
namespace vp_type
{

  namespace sin2
  {
    double w = 0.0227817;
    double E_0 = -0.0168;
    double n_c = 2.0;
    double t_on = 0.0;
    double t_off = 2.0 * M_PI * n_c / w;
    double phi = 0.0;

    double vpot( double time, int me )
    {
      double result = 0.0;
      
      if( ( time > t_on ) && ( time < t_off ) ) {
	double el = 0.5 * ( 1.0 - cos( w * time / n_c ) );
	result = - E_0 * el * sin( w * time + phi ) / w;
      }
      return result;
    }
  }

  // namespace sin
  // {
  //   double vpot( double time )
  //   {
  //     double result = 0.0;
  //     if( ( time > t_field_on ) && ( time < t_field_off ) ) {
  // 	result = - E_0 * sin( w * time ) / w;
  //     }
  //     return result;
  //   }
  // }
  
  // namespace ramp
  // {
  //   double vpot( double time )
  //   {
      
  //     double result, ampl, ramping, constperiod, downramp, T;
  //     result = 0.0;
  //     ampl   = E_0 / w;
      
  //     T = 2 * M_PI / w;
      
  //     constperiod = n_c * T;
  //     ramping     = n_r * T;
  //     downramp    = ramping;
      
  //     if( time > 0.0 ) {
  // 	if( time < ( t_field_on + ramping ) ) {
  // 	  result = - ampl / ramping * time * cos( w * time ) 
  // 	    + ampl / ( w * ramping ) * sin( w * time );
  // 	}
  // 	else if ( ( time >= ramping ) && ( time < ramping + constperiod ) )
  // 	  result = - ampl * cos( w * time );
	
  // 	else if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp)) {  
  // 	  result = ampl / downramp * ( time - downramp - constperiod - ramping )
  // 	    * cos( w * time ) - ampl / ( w * downramp ) * sin( w * time );
  // 	}
  // 	else result=0.0;
  //     }
  
  //     return result;
  //   }  
  // }

  namespace none
  {
    double vpot( double time, int me )
    {
      return 0;
    }
  }

}

*/


namespace T_vecpot
{
  double w = 0.0227817;
  double E_0 = -0.0168;
  double n_c = 0.5;
  double t_on = 0.0;
  double t_off = 2.0 * M_PI * n_c / w;
  double phi = 0.0;

  double none::vpot_( double time, int me ) {
    return 0.0;
  }

  double sin2::vpot_( double time, int me ) {
    double result = 0.0; 
    if( ( time > t_on ) && ( time < t_off ) ) {
      double el = 0.5 * ( 1.0 - cos( w * time / n_c ) );
      result = - E_0 * el * sin( w * time + phi ) / w;
    }
    return result;    
  }

}
