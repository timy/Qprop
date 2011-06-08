#include "toolbox.h"

#include "potentials.h"

int main( int argc, char* argv[] )
{
  grid g, g_load, g_h1s;
  wavefunction staticpot, wf, wf_h1s;
  int tag, me = 0;
  hamop H;
  
  // load the initial gaussian wave packet
  tag = 139;
  init_grid_from_file( g, g_load, tag );
  init_wavefunction_from_file( g, g_load, wf, tag );
  
  // load the hydrogen 1s state
  tag = 9;
  read_grid_from_file( g_h1s, tag );
  read_wavefunction_from_file( g_h1s, wf_h1s, tag );
 
  // load hamiltonian from "potentials.cc"
  init_hamilton( g, H, staticpot );

  // propagate
  tag = 139;
  complex timestep = complex( parameters::real_timestep, 0.0 );
  double time = parameters::t0;                  // set the starting time
  long output_interval = 100;

  // output laser pulse
  // export_vecpot( output_interval, tag+1 );
  // export_electric_field( output_interval, tag+1 );
  disp_elapsed_time();
  fprintf( stdout, "n_ts_re = %ld\n", parameters::n_ts_re ); 

  for( long ts = 0; ts < parameters::n_ts_re; ts ++ ) {
    time = time + parameters::real_timestep;
    wf.propagate( timestep, time, g, H, me, staticpot, 
		  0, parameters::nuclear_charge );
    export_observable( g, g_h1s, H, wf, wf_h1s, staticpot, 
		       ts, output_interval, 140 );
    
    if( ts % output_interval == 0 ) {
      tag = tag + 1;
      export_final_wf( g, wf, tag );
      export_spatial_wf( g, wf, tag );
      disp_elapsed_time();
    }
  }

  disp_elapsed_time();  
  done();
  return 0;
}
