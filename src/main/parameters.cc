#include "parameters.h"
#include "vp_type.h"
#include "sp_type.h"


namespace parameters
{ 
  //  double r_max = 800.0;
  double r_max = 400.0;
  
  double delta_r = 0.2;   // step size of radial grid
  
  long n_x = (long) round( ( r_max / delta_r ) );     // number of radial grid

  long n_l = 60;        // number of partial waves (expansion with SH)
  
  double I_p = 0.5;      // ionization potential

  double charge = 0.0;

  double time_step = ( abs( charge ) > 1.0 ) ? delta_r / charge / 4.0 : delta_r / 4.0;

  double field_period = 2.0 * M_PI / T_vecpot::w;    // period of field

  double field_duration = T_vecpot::n_c * field_period;    // the total time for field switched on
  
  //double t0 = 41.7614;    // the starting time for propagation
  double t0 = 0.0;

  double t_end = field_duration;   // the ending time for propagation

  double t_field_on = T_vecpot::t_on;    // the starting time of field

  double t_field_off = T_vecpot::t_off;      // the ending time of the field

  double duration  = t_end - t0;     // the total propagation time

  long n_ts = (long) duration / time_step + 1;     // the number of steps for real-time propagation

  clock_t start_time = clock();     // the time when the code starts to run, for profiling



}
