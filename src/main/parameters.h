#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <bar.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
#include <ctime>

namespace parameters
{

  extern double r_max;

  extern double  delta_r;
  
  extern long  n_x;

  extern long  n_l;

  extern double  I_p;

  extern double charge;

  extern double time_step;

  extern double t0;

  extern double t_end;

  extern double field_period;
  
  extern double field_duration;

  extern double t_field_on;

  extern double t_field_off;

  extern double  duration;    
  
  extern long  n_ts;

  extern clock_t start_time;

}

#endif
