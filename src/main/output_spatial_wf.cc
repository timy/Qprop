#include "toolbox.h"
#include "potentials.h"
#include <iostream>

int main(int argc, char **argv)
{
  grid g;
  wavefunction wf;
  int tag;

  disp_copyright();

  for( int i = 11; i < 211; i ++ ) {
    tag = i;
    read_grid_from_file( g, 11 );
    read_wavefunction_from_file( g, wf, tag );
    export_spatial_wf( g, wf, tag );
  }

  done();
  return 0;
}
