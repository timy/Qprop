#include "toolbox.h"
#include "parameters.h"
#include <ctime>
#include "ip_type.h"

namespace toolbox
{
  double ( *fp_imgpot ) ( long, long, long, double, grid ) = NULL;
  int file_log_tag = 0;
  FILE* file_log = NULL;

  int init_grid( grid &g )
  {
    g.set_dim( 34 ); 
    g.set_ngps( parameters::n_x, parameters::n_l, 1 );  // N_r, L, N
    g.set_delt( parameters::delta_r, 1.0, 0.0 );  // delta r
    g.set_offs( 0, 0, 0 );
    std::cout << "n_x = " << parameters::n_x << ", "
	      << "n_l = " << parameters::n_l << std::endl;
    std::cout << "delta_r = " << parameters::delta_r 
	      << std::endl;
    std::cout << std::string( 45, '<' ) << " init_grid" 
	      << std::endl << std::endl;
    return 0;
  }

  int read_grid( grid &g, int tag = 0)
  {
    long ngps_x, ngps_y, ngps_z, dimens;
    double delt_x;
    char cstr_fname_info[400];
    sprintf( cstr_fname_info, "./res/info-%d.dat", tag );
    fprintf( stdout, "start to read info from \"%s\"\n", cstr_fname_info );
    export_log( file_log_tag, false );
    fprintf( file_log, "\nread grid info from file \"%s\"\n", cstr_fname_info );
    FILE* file_info = fopen( cstr_fname_info, "r" );
    if ( file_info != NULL ) {
      xml_parse_fgetl( &ngps_x, file_info, "<ngps_x>", "</ngps_x>" );
      xml_parse_fgetl( &ngps_y, file_info, "<ngps_y>", "</ngps_y>" );
      xml_parse_fgetl( &ngps_z, file_info, "<ngps_z>", "</ngps_z>" );
      xml_parse_fgetd( &delt_x, file_info, "<delt_x>", "</delt_x>" );
      xml_parse_fgetl( &dimens, file_info, "<dimens>", "</dimens>" );
      fclose( file_info );
    }
    else {
      fprintf( stderr, "%s can't be opened!\n", cstr_fname_info );
      exit( -10 );
    }

    g.set_dim( dimens );
    g.set_ngps( ngps_x, ngps_y, ngps_z ); 
    g.set_delt( delt_x, 1.0, 0.0 );
    g.set_offs( 0, 0, 0 );
    std::cout << std::string( 45, '<' ) << " read_grid" 
	      << std::endl << std::endl;
    return 0;
  }

  /* initialize wavefunction with specified method */
  int init_wavefunction( grid &g, wavefunction &wf, C_wf &wfc )
  {
    fluid ells;

    init_grid( g );

    ells.init( g.ngps_z() );
    ells[0] = 0;

    wf.init( g.size() );
    if( wfc.init != NULL ) {
      fprintf( stdout, "wavefunction is set to %s\n", wfc.get_type() );
      export_log( file_log_tag, false );
      fprintf( file_log, "\ninit_wavefunction: wavefunction is set to %s\n", wfc.get_type() );
      wf.init( g, wfc.init );
    }
    else {
      fprintf( stderr, "initialization of wave function is not specified!\n" );
      exit( 0 );
      // TODO: imaginary time propagation will be used to initialize wavefunction here.
    }
    //wf.normalize( g );
    fprintf( stdout, "norm of ini_wf: %20.15lf\n", wf.norm(g) );
    std::cout << std::string( 45, '<' ) << " init_wavefunction" 
	      << std::endl << std::endl;
    return 0;
  }

  int read_wavefunction( grid &g, wavefunction &wf, int tag = 0 )
  {
    char cstr_fname_wf[400];
    int iv = 0;

    read_grid( g, tag );
    wf.init( g.size() );
    sprintf( cstr_fname_wf, "./res/wf-%d.dat", tag );
    fprintf( stdout, "start to read wavefunction from \"%s\"\n", cstr_fname_wf );
    
    export_log( file_log_tag, false );
    fprintf( file_log, "\nread wavefunction: from file \"%s\"\n", cstr_fname_wf );

    FILE* file_wf = fopen(cstr_fname_wf, "r");
    if ( file_wf == NULL ) {
      fprintf( stderr, "err: %s cant be opened!\n", cstr_fname_wf ); 
      exit(-4);
    }     
    wf.init( g, file_wf, 0, iv );
    fprintf( stdout, "norm of wf-%d: %20.15lf\n", tag, wf.norm(g) );
    fprintf( file_log, "norm of wf-%d: %20.15lf\n", tag, wf.norm(g) );
    std::cout << std::string( 45, '<' ) << " read_wavefunction" 
	      << std::endl << std::endl;
    return 0;
  }

  int init_wavefunction_from_file( grid &g, grid &g_load, wavefunction &wf, int tag = 0 )
  {
    wavefunction wf_load;

    init_grid( g );
    wf.init( g.size() );
    read_wavefunction( g_load, wf_load, tag );
    wf.regrid( g, g_load, wf_load );
    fprintf( stdout, "regrid wf: (%ld, %ld)->(%ld, %ld)\n", 
	     g_load.ngps_x(), g_load.ngps_y(), g.ngps_x(), g.ngps_y() );
    fprintf( stdout, "norm of ini_wf: %20.15lf\n", wf.norm(g) );

    export_log( file_log_tag, false );
    fprintf( file_log, "\ninit_wavefunction_from_file\n" );
    fprintf( file_log, "regrid wf: (%ld, %ld)->(%ld, %ld)\n", 
	     g_load.ngps_x(), g_load.ngps_y(), g.ngps_x(), g.ngps_y() );

    std::cout << std::string( 45, '<' ) << " init_wavefunction_from_file" 
	      << std::endl << std::endl;
    return 0;
  }

  int set_imgpot( double (*fp_ip) (long, long, long, double, grid) )
  {
    fp_imgpot = fp_ip;
    return 0;
  }

  void set_charge( double charge )
  {
    parameters::charge = charge;
    export_log( file_log_tag, false );
    fprintf( file_log, "\ncharge is reset to %lf\n", parameters::charge );
  }

  int init_hamilton( grid &g, hamop &hamilton, wavefunction &staticpot,
		     C_vecpot &vpx, C_vecpot &vpy, C_vecpot &vpz, 
		     C_sclpot &spx, C_sclpot &spy, C_sclpot &spz,
		     C_field &field
		     )
  {
    fprintf( stdout, "vecpot_x is set to %s\n", vpx.get_type() );
    fprintf( stdout, "vecpot_y is set to %s\n", vpy.get_type() );
    fprintf( stdout, "vecpot_z is set to %s\n", vpz.get_type() );

    fprintf( stdout, "sclpot_x is set to %s\n", spx.get_type() );
    fprintf( stdout, "sclpot_y is set to %s\n", spy.get_type() );
    fprintf( stdout, "sclpot_z is set to %s\n", spz.get_type() );

    fprintf( stdout, "field    is set to %s\n", field.get_type() );

    export_log( file_log_tag, false );
    fprintf( file_log, "\nHamiltonian Initialization:\n" );
    fprintf( file_log, "vecpot_x is set to %s\n", vpx.get_type() );
    fprintf( file_log, "vecpot_y is set to %s\n", vpy.get_type() );
    fprintf( file_log, "vecpot_z is set to %s\n", vpz.get_type() );

    fprintf( file_log, "sclpot_x is set to %s\n", spx.get_type() );
    fprintf( file_log, "sclpot_y is set to %s\n", spy.get_type() );
    fprintf( file_log, "sclpot_z is set to %s\n", spz.get_type() );
    
    fprintf( file_log, "field    is set to %s\n", field.get_type() );

    if( fp_imgpot == NULL ) {
      fp_imgpot = ip_type::ipot;
      fprintf( stdout, "imgpot is set to ip_type::ipot by default\n" );
      fprintf( file_log, "imgpot is set to ip_type::ipot by default\n" );

    }

    hamilton.init( g, 
		   vpx.vpot, vpy.vpot, vpz.vpot, 
		   spx.spot, spy.spot, spz.spot, 
		   fp_imgpot, field.field );

    // this is the linear and constant part of the Hamiltonian
    staticpot.init( g.size() ); 
    staticpot.calculate_staticpot( g, hamilton );
    std::cout << std::string( 45, '<' ) << " init_hamiltonian" 
	      << std::endl << std::endl;
    return 0;
  }

  int export_wf( grid &g, wavefunction &wf, int tag = 0 )
  {
    char cstr_fname_wf[400];
    sprintf( cstr_fname_wf, "./res/wf-%d.dat", tag );
    fprintf( stdout, "start to write to %s\n", cstr_fname_wf );
    FILE* file_wf = fopen( cstr_fname_wf, "w" );
    if (file_wf != NULL) {
      wf.dump_to_file_sh( g, file_wf, 1 );
      fprintf( stdout, "norm = %20.15lf\n", wf.norm(g) );
      fclose( file_wf );
      std::cout << std::string( 45, '<' ) << " export_wavefunction" 
		<< std::endl << std::endl;
      export_info( g, tag );
      return 0;
    }
    else {
      fprintf( stderr, "%s cant be opened!\n", cstr_fname_wf );
      return -5;
    }
  }

  int export_spatial_wf( grid &g, wavefunction &wf, int tag=0 )
  {
    fluid degeneracies = 0, ms = 0;
    double factor = 0.4;
    int n_steps = 60;
    int stepwidth = (int) ( factor * g.ngps_x() / n_steps );
    std::cout << "stepwidth for spatial wavepacket output = " << stepwidth << std::endl;
    char cstr_fname_wf_spa[400];
    sprintf( cstr_fname_wf_spa, "./res/wf_spa-%d.dat", tag );

    FILE* file_wf_spa = fopen( cstr_fname_wf_spa, "w" );
    if ( file_wf_spa != NULL ) {
      wf.dump_to_file( g, file_wf_spa, stepwidth, degeneracies, ms, factor);
      fclose( file_wf_spa );
      std::cout << std::string( 45, '<' ) << " export_spatial_wavefunction" 
		<< std::endl << std::endl;
      return 0;
    }
    else {
      fprintf( stderr, "%s cant be opened!\n", cstr_fname_wf_spa );
      return -5;
    }
  }

  int export_observable( grid &g, hamop &hamilton, 
			 wavefunction &wf, wavefunction &staticpot, 
			 long ts, long interval, int tag = 0 )
  {
    char cstr_fname_obser[400];
    double E_tot;
    static FILE* file_obser;
    int me = 0;

    if( ts == 0 ) {
      sprintf( cstr_fname_obser, "./res/obser-%d.dat", tag );
      file_obser = fopen( cstr_fname_obser, "w" );
    }

    if( file_obser != NULL ) {
      if( ts % interval == 0 ) {
	E_tot = real( wf.energy( 0.0, g, hamilton, me, staticpot, 
				 parameters::charge ) );
	fprintf( file_obser, "%6ld   % 20.15le\n", ts, E_tot);
	fprintf( stdout, "%ld: E = %20.15lf\n", ts, E_tot );
      }
      return 0;
    }
    else {
      fprintf( stderr, "%s can\'t be opened!\n", cstr_fname_obser );
      exit(-3);
    }
  
    if( ts == parameters::n_ts - 1 ) {
      fclose( file_obser );
      std::cout << std::string( 45, '<' ) << " export_observable" 
		<< std::endl << std::endl;
    
    }
  
  }

  int export_observable( grid &g, grid &g_proj, hamop &hamilton, wavefunction &wf, wavefunction &wf_proj, 
			 wavefunction &staticpot, long ts, long interval, int tag = 0 )
  {
    char cstr_fname_obser[400];
    double time, E_tot, N, z_expect;
    static FILE* file_obser;
    int me = 0;
    wavefunction P;

    if( ts == 0 ) {
      sprintf( cstr_fname_obser, "./res/obser-%d.dat", tag );
      file_obser = fopen( cstr_fname_obser, "w" );
    }

    if (file_obser != NULL) {
      if ( ts % interval == 0 ) {
	time = parameters::t0 + ts * parameters::time_step;
	E_tot = real( wf.energy( 0.0, g, hamilton, me, staticpot, parameters::charge ) );
	P     = wf.project( g, g_proj, wf_proj );
	N     = wf.norm( g );
	z_expect = real( wf.expect_z( g ) );
 
	fprintf( file_obser, "%15.10le %15.10le %15.10le %15.10le %15.10le %ld\n",
		 time, E_tot, real(conj(P[0])*P[0]), N, z_expect, ts );
	
	fprintf( stdout,     "%15.10le %15.10le %15.10le %15.10le %15.10le %ld\n", 
		 time, E_tot, real(conj(P[0])*P[0]), N, z_expect, ts );
      }
      return 0;
    }
    else {
      fprintf( stderr, "%s cant be opened!\n", cstr_fname_obser );
      exit(-3);
    }

    if( ts == parameters::n_ts-1 ) {
      fclose( file_obser );
      std::cout << std::string( 45, '<' ) << " export_observable" 
		<< std::endl << std::endl;
    }

  }


  void set_log_tag( int tag )
  {
    file_log_tag = tag;
    return;
  }
  
  int export_log( int tag, bool b_close=false )
  {
    if( file_log == NULL ) {
      char cstr_fname_log[400];
      sprintf( cstr_fname_log, "./res/log-%d.dat", tag );
      fprintf( stdout, "start to write to %s\n", cstr_fname_log );
      time_t current_time;
      time( &current_time );
      
      file_log = fopen( cstr_fname_log, "w" );
      if( file_log != NULL ) {
	fprintf( file_log, "Qprop LOG recorded at %s\n\n", ctime(&current_time) );
	
	fprintf( file_log, "Parameters for the Grid\n" );
	fprintf( file_log, "n_r       = %ld\n", parameters::n_x );
	fprintf( file_log, "n_l       = %ld\n", parameters::n_l );
	fprintf( file_log, "n_orbital = %ld\n", (long) 1 );
	fprintf( file_log, "g.dimens  = %ld\n", (long) 34 );
	fprintf( file_log, "delta_r   = %lf\n", parameters::delta_r );
	fprintf( file_log, "r_max     = %lf\n\n", parameters::r_max );
	
	fprintf( file_log, "charge    = %lf\n\n", parameters::charge );
	
	fprintf( file_log, "Parameters for the Propagation Time (if propagate)\n");
	fprintf( file_log, "t_0       = %lf\n", parameters::t0 );
	fprintf( file_log, "t_end     = %lf\n", parameters::t_end );
	fprintf( file_log, "duration  = %lf\n", parameters::duration );
	fprintf( file_log, "time_step = %lf\n", parameters::time_step );
	fprintf( file_log, "n_ts      = %ld\n\n", parameters::n_ts );
	
	fprintf( file_log, "Parameters for the Vector Potential (if exist)\n" );
	fprintf( file_log, "t_field_on  = %lf\n", parameters::t_field_on );
	fprintf( file_log, "t_field_off = %lf\n", parameters::t_field_off );
	fprintf( file_log, "field_period = %lf\n", parameters::field_period );
	fprintf( file_log, "field_duration = %lf\n", parameters::field_duration );
	fprintf( file_log, "w         = %lf\n", T_vecpot::w );
	fprintf( file_log, "E_0       = %lf\n", T_vecpot::E_0 );
	fprintf( file_log, "n_c       = %lf\n", T_vecpot::n_c );
	fprintf( file_log, "phi       = %lf\n", T_vecpot::phi );
      }
      else {
	fprintf( stderr, "Cannot open file %s\n", cstr_fname_log );
	return -1;
      }
      return 0;
    }

    if( b_close ) {
      fflush( file_log );
      fclose( file_log );
      file_log = NULL;
      std::cout << std::string( 45, '<' ) << " export_log" 
		<< std::endl << std::endl;
    }
    return 0;
  }

  int export_info( grid &g, int tag = 0 )
  {
    char cstr_fname_info[400];
    sprintf( cstr_fname_info,  "./res/info-%d.dat", tag );
    FILE* file_info = fopen( cstr_fname_info, "w" );
    fprintf( stdout, "start to write to %s\n", cstr_fname_info );
    if( file_info != NULL ) {
      fprintf( file_info, "<ngps_x>%ld</ngps_x>\n", g.ngps_x() );
      fprintf( file_info, "<ngps_y>%ld</ngps_y>\n", g.ngps_y() );
      fprintf( file_info, "<ngps_z>%ld</ngps_z>\n", g.ngps_z() );
      fprintf( file_info, "<delt_x>%20.15le</delt_x>\n", g.delt_x() );
      fprintf( file_info, "<dimens>%ld</dimens>\n", (long)g.dimens() );
      fflush( file_info );
      fclose( file_info );
      std::cout << std::string( 45, '<' ) << " export_info" 
		<< std::endl << std::endl;
      return 0;
    }
    else {
      fprintf( stderr, "Cannot open file %s\n", cstr_fname_info );
      return -1;
    }
  }


  int export_vecpot( C_vecpot &vp, long interval, int tag = 0 )
  {
    char cstr_fname_vpot[400];
    double time = parameters::t0;
    int me = 0;

    if( vp.vpot == NULL ) {
      fprintf( stderr, "export_vecpot: %s.vpot = NULL!\n", vp.get_type() );
      exit( -1 );
    }

    sprintf( cstr_fname_vpot,  "./res/vpot-%d.dat", tag );
    fprintf( stdout, "start to write to %s\n", cstr_fname_vpot );
    FILE* file_vpot = fopen( cstr_fname_vpot, "w" );
    if( file_vpot != NULL ) {
      for ( long ts = 0; ts < parameters::n_ts; ts ++ ) {
	time = time + parameters::time_step;
	if ( ts % interval == 0 ) {
	  fprintf( file_vpot, "%15.10le %15.10le\n", 
		   time, vp.vpot( time, me ) );
	}
      }
      fclose( file_vpot );
      std::cout << std::string( 45, '<' ) << " export_vecpot" 
		<< std::endl << std::endl;
    }
    else {
      fprintf(stdout, "%s cant be opened!\n", cstr_fname_vpot); 
      return -1;
    }
    return 0;
  }

  int export_field( C_vecpot &vp, long interval, int tag = 0 )
  {
    char cstr_fname_elec[400];
    double time = parameters::t0;
    int me = 0;

    if( vp.vpot == NULL ) {
      fprintf( stderr, "export_electric_field: %s.vpot = NULL!\n", vp.get_type() );
      exit( -1 );
    }

    sprintf( cstr_fname_elec,  "./res/elec-%d.dat", tag );
    fprintf( stdout, "start to write to %s\n", cstr_fname_elec );
    FILE* file_elec = fopen( cstr_fname_elec, "w" );
  
    if( file_elec != NULL ) {
      for ( long ts = 1; ts < parameters::n_ts; ts ++ ) {
	time = time + parameters::time_step;
	if ( ts % interval == 0 ) {
	  double elec = ( vp.vpot( time - 0.5 * parameters::time_step, me ) 
			  - vp.vpot( time + 0.5 * parameters::time_step, me ) )
	    / parameters::time_step;
	  fprintf( file_elec, "%15.10le %15.10le\n", time, elec );
	}
      }
      fclose(file_elec);
      std::cout << std::string( 45, '<' ) << " export_electric_field" 
		<< std::endl << std::endl;
    }
    else {
      fprintf( stdout, "%s cant be opened!\n", cstr_fname_elec ); 
      return -1;
    }
    return 0;
  }

  int export_field( C_field &fd, long interval, int tag )
  {
    char cstr_fname_field[400];
    double time = parameters::t0;
    int me = 0;

    if( fd.field == NULL ) {
      fprintf( stderr, "export_field: %s.field = NULL!\n", fd.get_type() );
      exit( -1 );
    }

    sprintf( cstr_fname_field, "./res/field-%d.dat", tag );
    fprintf( stdout, "start to write to %s\n", cstr_fname_field );
    FILE* file_field = fopen( cstr_fname_field, "w" );
    
    if( file_field != NULL ) {
      for( long ts = 0; ts < parameters::n_ts; ts ++ ) {
	time = time + parameters:: time_step;
	if( ts % interval == 0 ) {
	  double field = fd.field( time, me );
	  fprintf( file_field, "%20.15le %20.15le\n", time, field );
	}
      }
      fclose( file_field );
      std::cout << std::string( 45, '<' ) << "export_field" << std::endl << std::endl;
    }
    else {
      fprintf( stdout, "%s cant be opened!\n", cstr_fname_field ); 
      return -1;
    }
    return 0;    
  }

  void disp_copyright()
  {
    fprintf( stdout, 
	     " ---------------------   QProp   -------------------------\n" );
    fprintf( stdout, 
	     " (C) Copyright by Bauer D and Koval P, Heidelberg (2005)  \n" );
    fprintf( stdout,
	     " ---------------------------------------------------------\n" );
    return;
  }

  void done()
  {
    std::cout << "Hasta la vista" << std::endl;
    exit(0);
  }

  void disp_elapsed_time()
  {
    double elapsed_time = ( clock() - parameters::start_time ) 
      / ( (double) CLOCKS_PER_SEC );
    fprintf( stdout, "Time elapsed: %f\n", elapsed_time );
    return;
  }

}
