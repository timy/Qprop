#include "ip_type.h"

namespace ip_type
{
  
  // amplitude of imaginary absorbing potential
  double  ampl_im = 100.0;

  double ipot(long xindex, long yindex, long zindex, double time, grid g)
  {
    double x,y,z;

    if (ampl_im>1.0) {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
      
      switch(g.dimens()) {
      case 3 : 
	y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	  *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
	z=((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z())
	  *((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z());
	return ampl_im*(x*x*x*x*x*x*x*x+y*y*y*y*y*y*y*y+z*z*z*z*z*z*z*z);
	break;
	  
      case 2 : case 6 :
	y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	  *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
	return ampl_im*(x*x*x*x*x*x*x*x+y*y*y*y*y*y*y*y);
	break;
	  
      case 1 : case 5 : case 15 :
	return ampl_im*(x*x*x*x*x*x*x*x);
	break;
	  
      case 4 : case 14 : case 24 : case 34 :		
	x=((double) xindex + 0.5  )/(g.ngps_x())
	  *((double) xindex + 0.5  )/(g.ngps_x());
	return ampl_im*(x*x*x*x*x*x*x*x);
	break;
	  
      default : return 0.0;
      }
    }
    else {
      return 0.0;
    };
    
  } 
}

