#include "toolbox.h"
#include "parameters.h"

double my_vpot( double time, int me )
{
  return ( T_vecpot::sin2::vpot_( time, me ) - 
	   T_vecpot::sin2::vpot_( parameters::t0, me ) );
}

int main( int argc, char* argv[] )
{
  grid g, g_load, g_h1s;
  wavefunction staticpot, wf, wf_h1s;
  int tag, me = 0;
  hamop H;

  toolbox::set_log_tag( 100 );

  // load the initial gaussian wave packet
  tag = 0;
  toolbox::init_wavefunction_from_file( g, g_load, wf, tag );
  
  // load the hydrogen 1s state
  tag = 9;
  toolbox::read_wavefunction( g_h1s, wf_h1s, tag );
 
  // initialize the hamiltonian
  T_vecpot::none vpx, vpy;
  T_vecpot::user vpz( my_vpot );
  T_sclpot::none spx, spy, spz;
  T_field::none field;
  toolbox::init_hamilton( g, H, staticpot, vpx, vpy, vpz, spx, spy, spz, field );

  complex time_step = complex( parameters::time_step, 0.0 );
  double time = parameters::t0;
  long output_interval = 100;

  // output laser pulse
  tag = 100;
  toolbox::export_vecpot( vpz, output_interval, tag );
  toolbox::export_field( vpz, output_interval, tag );
  toolbox::disp_elapsed_time();
  fprintf( stdout, "n_ts_re = %ld\n", parameters::n_ts );

  toolbox::export_log( 100, true );
  // propagate
  for( long ts = 0; ts < parameters::n_ts; ts ++ ) {
    time = time + parameters::time_step;
    //  fprintf(stdout, "time = %lf, vpot = %lf\n", time, vpz.vpot(time, 0) );
    wf.propagate( time_step, time, g, H, me, staticpot, 
		  0, parameters::charge );
    toolbox::export_observable( g, g_h1s, H, wf, wf_h1s, staticpot, 
				ts, output_interval, 100 );

    if( ts % output_interval == 0 ) {
      tag = tag + 1;
      toolbox::export_wf( g, wf, tag );
      //toolbox::export_spatial_wf( g, wf, tag );
      toolbox::disp_elapsed_time();
    }
  }

  toolbox::disp_elapsed_time();  
  toolbox::done();
  return 0;
}
