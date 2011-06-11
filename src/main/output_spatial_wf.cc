#include "toolbox.h"
#include "parameters.h"

int main(int argc, char **argv)
{
  grid g;
  wavefunction wf;
  int tag;

  disp_copyright();

  for( int i = 11; i < 211; i ++ ) {
    tag = i;
    read_wavefunction( g, wf, tag );
    export_spatial_wf( g, wf, tag );
  }

  done();
  return 0;
}
