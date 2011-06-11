#include "toolbox.h"
#include "parameters.h"

double my_vpot( double time, int me )
{
  return ( T_vecpot::sin2::vpot_( time, me ) - 
	   T_vecpot::sin2::vpot_( parameters::t0, me ) );
}

int main( int argc, char* argv[] )
{
  int tag = 217;
  toolbox::set_log_tag( tag );
  long n_E = 400; 
  long n_a = 100;
  double E_0 = 0.0; 
  double E_1 = 0.5;

  // initialize the hamiltonian
  T_vecpot::none vpx, vpy, vpz;
  //  T_vecpot::user vpz( my_vpot );
  T_sclpot::coulomb spx;
  T_sclpot::none spy, spz;
  T_field::none field;
  toolbox::set_charge( 1.0 );
  
  toolbox::winop_spectra( tag, n_E, n_a, E_0, E_1, vpx, vpy, vpz, spx, spy, spz, field );

  toolbox::export_log( tag, true ); 
  toolbox::disp_elapsed_time();  
  toolbox::done();
  return 0;
}
