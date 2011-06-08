#include "toolbox.h"
#include "parameters.h"

int main( int argc, char* argv[] )
{
  grid g;
  wavefunction wf, staticpot;
  hamop H;

  toolbox::disp_copyright();

  int tag = 9;
  toolbox::set_log_tag( tag );

  // put the hydrogen 1s state on the grid
  T_wf::h1s wf_h1s;
  toolbox::init_wavefunction( g, wf, wf_h1s );
  wf.normalize( g );

  // prepare the hamiltonian for calculation of energy
  T_vecpot::none vpx, vpy, vpz;
  T_sclpot::coulomb spx; 
  T_sclpot::none spy, spz;
  T_field::none field;

  toolbox::set_charge( 1.0 );
  toolbox::init_hamilton( g, H, staticpot, vpx, vpy, vpz, spx, spy, spz, field );

  // the norm and the energy
  int me = 0;
  double norm = wf.norm( g );
  double energy = real( wf.energy( 0.0, g, H, me, staticpot, parameters::charge ) );
  fprintf( stdout, "norm = %20.15lf, energy of ini_wf = %20.15lf\n", norm, energy );

  toolbox::export_wf( g, wf, tag );
  toolbox::export_log( tag, true );
  toolbox::done();
  
  return 0;
}
