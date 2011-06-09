#include <iostream>
#include <stdio.h>

#include "wavefunction.h"
#include "fluid.h"
#include "grid.h"
#include "hamop.h"
#include "winop.h"
#include "xml_parse.h"
#include "ylm.h"

#include "../main/potentials.cc"

int main(int argc, char **argv)
{
  int  iv;

  long ngps_x, ngps_y, ngps_z, dimens, index, i, no_energ, no_theta,
    ltheta, no_dir, l_index, ylm_index;

  double delt_x, delta_energy, energy_ini, energy_fin, energ, gamma, w;

  complex summingup, result_tot;

  FILE *file_wf, *file_res, *file_info;

  char cstr_fname_res[1000], cstr_dname_wf[1000],
    cstr_fname_info[1000], cstr_fname_wf[1000];

  hamop  hamilton;

  fluid  theta, V_ee_0;
  grid  g_angle, g_load, g;
  wavefunction fullchi, angle_result, result_lsub, ylm_array,
    staticpot, wf, wf_load;

  no_theta            = 2; // no_theta angles between 0 and (1-1/no_theta)*pi
  no_energ            = 1000;    
  energy_ini         =  -0.6;    
  energy_fin         =  1.6;
  gamma               = (energy_fin - energy_ini) / ( 2.0*(no_energ-1.0));
  std::cout << "gamma = " << gamma << std::endl;
  //exit(0); //gamma = 0.000500626

  sprintf(cstr_dname_wf, "../main/res");
  sprintf(cstr_fname_wf, "%s/wf_fin-115.dat", cstr_dname_wf);
  sprintf(cstr_fname_info, "%s/info-101.dat", cstr_dname_wf);

  no_dir              = no_theta;
  delta_energy        = 2*gamma; // width of the bin is 2*gamma



  file_info = fopen(cstr_fname_info, "r");  
  if (file_info==NULL) {
    fprintf(stdout, "%s can't be opened!\n", cstr_fname_info); exit(-10);}  

  xml_parse_fgetd( &parameters::nuclear_charge, file_info,"<nuclear_charge>",
		  "</nuclear_charge>");
  xml_parse_fgetl(&ngps_x, file_info,"<ngps_x>", "</ngps_x>");
  xml_parse_fgetl(&ngps_y, file_info,"<ngps_y>", "</ngps_y>");
  xml_parse_fgetl(&ngps_z, file_info,"<ngps_z>", "</ngps_z>");
  xml_parse_fgetd(&delt_x,file_info,"<delt_x>","</delt_x>");
  xml_parse_fgetl(&dimens,file_info,"<dimens>","</dimens>");
  xml_parse_fgetd(&w,file_info,"<w>","</w>");
  fclose(file_info);

  sprintf(cstr_fname_res,  "./res/res-%lf-%ld-%ld-%ld-%.3le-%.3le.dat",
	  w, dimens, ngps_x, ngps_y, delta_energy, gamma);

  iv = 1;
  fprintf(stdout, " WINOP: Calculation of one-electron spectra\n");
  fprintf(stdout, " (C) Copyright by Bauer D and Koval P, Heidelberg (2005)\n");
  fprintf(stdout, " -------------------------------------------------------\n");

  file_wf  = fopen(cstr_fname_wf,"r");
  file_res = fopen(cstr_fname_res,"w");
 
  fprintf(stdout, "%s will be read\n", cstr_fname_wf);
  fprintf(stdout, "%s will be written\n", cstr_fname_res);

  if( file_wf == NULL ) {
    fprintf( stderr, "%s can't be opened!\n", cstr_fname_wf ); 
    exit(-10);
  }

  // specify angles \theta, from 0 to PI
  theta.init( no_theta );   
  for( ltheta=0; ltheta < no_theta; ltheta ++ )
    theta[ltheta] = ltheta * M_PI / no_theta;

  // this grid is used to organize the Y_lm coefficient
  g_angle.set_dim( 44 );
  g_angle.set_ngps( no_dir, ngps_y, 0 ); 
  g_angle.set_delt( 0.0 );
  g_angle.set_offs( 0, 0, 0 );

  // this is the grid for the wavefunction that is to be analyzed
  // (data is taken from the info file)
  g_load.set_dim( dimens );
  g_load.set_ngps( ngps_x, ngps_y, ngps_z );
  g_load.set_delt( delt_x );
  g_load.set_offs( 0, 0, 0);

  // this is the grid on which the wavefunction is analyzed
  // delt_x must be the same for g and g_load
  // increasing ngps_x improves the resolution in the continuum 
  g.set_dim(dimens);
  g.set_ngps(ngps_x,ngps_y,ngps_z); 
  g.set_delt(delt_x);
  g.set_offs(0,0,0);
  
  // angle-resolved wave function
  angle_result.init( no_dir );

  // radial-angle-resolved wave function
  fullchi.init( g.ngps_x()*g.ngps_y() ); 

  // ell-resolved wave function ( partial waves )
  result_lsub.init( g.ngps_y() );
  
  // hartree potential
  V_ee_0.init( g.ngps_x() );

  // array for Ylm coefficients
  ylm_array.init( g_angle.size() );

  // the Hamiltonian
  hamilton.init( g, parameters::vecpot_x, parameters::vecpot_y, parameters::vecpot_z, 
		 parameters::scalarpotx, parameters::scalarpoty, parameters::scalarpotz, 
		 parameters::imagpot, parameters::field );

  // static potential
  staticpot.init( g.size() );
  staticpot.calculate_staticpot( g, hamilton );

  // the wavefunction array for analysis
  wf.init( g.size() );

  // the wavefunction array for loading
  wf_load.init( g_load.size() );

  // calculate coefficient array for Ylm
  for( l_index = 0; l_index < g_angle.ngps_y(); l_index ++ ) {
    for( ltheta = 0; ltheta < no_theta; ltheta ++ ) {
      ylm_index = g_angle.rlmindex( ltheta, l_index, 0 );
      ylm_array[ylm_index] = ylm( l_index, 0, theta[ltheta], 0.0 );  
    }
  }

  wf_load.init( g_load, file_wf, 0, iv );
  fclose( file_wf );

  std::cout << "Reading from file: done." << std::endl;
  
  wf.regrid( g, g_load, wf_load );
  
  
  fprintf(stdout, "Grid: \n", g.ngps_x());
  fprintf(stdout, "g.ngps_x() = %ld\n", g.ngps_x());
  fprintf(stdout, "g.ngps_y() = %ld\n", g.ngps_y());
  fprintf(stdout, "g.ngps_z() = %ld\n", g.ngps_z());
  fprintf(stdout, "g.dimens() = %d\n", g.dimens());
  fprintf(stdout, "g.delt_x() = %20.15le\n", g.delt_x());
  fprintf(stdout, "nuclear_charge    = %20.15le\n", parameters::nuclear_charge);
  fprintf(stdout, "cstr_fname_wf = %s\n", cstr_fname_wf);
  fflush(stdout);

  for( long lenergy = 0; lenergy < no_energ; lenergy ++ ) {
    energ = energy_ini + lenergy * delta_energy;

    // write the current energy
    fprintf( file_res, "%.15le ", energ );

    winop_fullchi( fullchi, result_lsub, &result_tot, energ, gamma,
		   staticpot, V_ee_0, parameters::nuclear_charge, g, wf, iv );

    // write the partial waves ( for individual l ) to file
    for( l_index = 0; l_index < g.ngps_y(); l_index ++ ) {
      fprintf( file_res, "%.15le ", real(result_lsub[l_index]) );
    }
    // write the total result to file
    fprintf( file_res, "%.15le ", real(result_tot) );

    angle_result.nullify();
    for( ltheta = 0; ltheta < no_theta; ltheta ++ ) {
      for( long rindex = 0; rindex < g.ngps_x(); rindex ++ ) {
	summingup = complex( 0.0, 0.0 );
	for( l_index = 0; l_index < g.ngps_y(); l_index ++ ) {
	  index = g.index( rindex, l_index, 0, 0 );
	  summingup += fullchi[index] * ylm_array[g_angle.rlmindex( ltheta, l_index, 0 )];
	}
	angle_result[ltheta] += summingup * conj(summingup);
      }
      angle_result[ltheta] *= g.delt_x();
      // write the energy-angle-resolved probability to file
      fprintf( file_res, "%.15le ", real(angle_result[ltheta]) );
    }
    fprintf(file_res, "\n");
    fflush(file_res);
  } // end of energy loop
  fclose(file_res);
  printf("%s is written\n", cstr_fname_res);
  std::cout << "Hasta la vista... " << std::endl;
}
