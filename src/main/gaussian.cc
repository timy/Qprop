#include "toolbox.h"
#include "parameters.h"

int main( int argc, char* argv[] )
{
  grid g;
  wavefunction wf, staticpot;
  hamop H;

  toolbox::disp_copyright();
  
  int tag = 0;
  toolbox::set_log_tag( tag );
  
  T_wf::gaussian wf_gauss;
  toolbox::init_wavefunction( g, wf, wf_gauss );
  wf.normalize( g );

  T_vecpot::none vpx, vpy, vpz; 
  T_sclpot::none spx, spy, spz;
  T_field::none field;
  toolbox::init_hamilton( g, H, staticpot, vpx, vpy, vpz, spx, spy, spz, field );
  std::cout << "charge = " << parameters::charge << std::endl;

  int me = 0;
  double norm = wf.norm( g );
  double energy = real( wf.energy( 0.0, g, H, me, staticpot, parameters::charge ) );
  fprintf( stdout, "norm = %20.15lf, energy of ini_wf = %20.15lf\n", norm, energy );

  toolbox::export_wf( g, wf, tag );
  toolbox::export_log( tag, true );
  toolbox::done();
    
  return 0;
}
