#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <stdio.h>
#include <complex>
#define complex std::complex<double>
#include <bar.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
#include <xml_parse.h>

#include "vp_type.h"
#include "sp_type.h"
#include "wf_type.h"
#include "fd_type.h"

namespace toolbox
{
  int init_wavefunction( grid &g, wavefunction &wf, C_wf &wfc );

  int read_wavefunction( grid &g, wavefunction &wf, int tag );

  int init_wavefunction_from_file( grid &g, grid &g_load, wavefunction &wf, int tag );

  int set_imgpot( double (*fp_ip) (long, long, long, double, grid) );

  void set_charge( double charge );

  int init_hamilton( grid &g, hamop &hamilton, wavefunction &staticpot,  
		     C_vecpot &vpx, C_vecpot &vpy, C_vecpot &vpz, 
		     C_sclpot &spx, C_sclpot &spy, C_sclpot &spz, 
		     C_field &field );

  int export_wf( grid &g, wavefunction &wf, int tag );

  int export_spatial_wf( grid &g, wavefunction &wf, int tag ) ;

  int export_observable( grid &g, hamop &hamilton, wavefunction &wf, wavefunction &staticpot, 
			 long ts, long interval, int tag );
  int export_observable( grid &g, grid &g_load, hamop &hamilton, 
			 wavefunction &wf, wavefunction &wf_load, 
			 wavefunction &staticpot, long ts, long interval, int tag );

  void set_log_tag( int tag );

  int export_log( int tag, bool b_close );

  int export_info( grid &g, int tag);

  int export_vecpot( C_vecpot &vp, long interval, int tag );
  int export_field( C_vecpot &vp, long interval, int tag );
  int export_field( C_field &fd, long interval, int tag );

  void disp_copyright();
  void done();
  void disp_elapsed_time();

}

#endif
