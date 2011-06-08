#ifndef winop_h
#define winop_h winop_h

#include<wavefunction.h>
#include<grid.h>
#include<fluid.h>

class grid;
class hamop;
class wavefunction;

//
// Applies a window operator W(energy, gamma) to the wavefunction
//
int winop_fullchi(wavefunction &fullchi, wavefunction &spec_result_lsub,
  complex *spec_result_tot, double energ, double energ_width,
  wavefunction staticpot, fluid hartree_potential_zero,
  double nuclear_charge, grid g, wavefunction wf, int iv);

#endif // winop_h
