#ifndef fluid_h
#define fluid_h fluid_h
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <string>

class grid;
class hamop;
class wavefunction;


class fluid
{
 public:
  fluid(long x=0) : fluid_dim(x), start(new double[x]) { }
  
  fluid(const fluid& v) 
    {
      fluid_dim  = v.fluid_dim;
      start = new double[fluid_dim];
      for (long i = 0; i < fluid_dim; i++)
	start[i] = v.start[i];
    }
  
  virtual ~fluid() { delete [] start;}
  
  long  wf_size() const {return fluid_dim;}
  
  double&  operator[](long index)
    {
      //      assert(index >= 0  &&  index < fluid_dim);
      return start[index];
    }
  
  const double&  operator[](long index) const   
    {   
      //      assert(index >= 0  &&  index < fluid_dim);
      return start[index];
    }
  
  fluid& operator=(const fluid&);
  fluid& operator=(const wavefunction&);
  
  double* begin() { return start;}
  double* end()   { return start + fluid_dim;}
  const double* begin() const { return start;}
  const double* end()   const { return start + fluid_dim;}
  
  fluid& operator *= (double z); 
  
  void nullify();
  
  void dump_to_file_sh(grid g,  FILE* os, int stepwidth);  
  void load(FILE* os, int iv);  

  int  init(long isize);  
  double integrate_along_x(grid g);
  double expect_rr(grid g);
  double expect_oneoverr(grid g);

  double calculate_ionkinenergy(const fluid &Q, const double singlecharge, const double singlemass);
  double calculate_ionpotenergy(const fluid &Q);

 void embed_as_x(grid g, long yindex, long zindex, const fluid &v);
 fluid extract_x(grid g, long other_one, long other_two);
  
    
 private:
  long  fluid_dim;
  double   *start; 
  
};


std::ostream& operator<<( std::ostream& os, const fluid& v);
//istream& operator>>(istream& is, fluid& v);
//istream& operator>>(ifstream& is, fluid& v);
double operator * (const fluid &v, const fluid &w );
fluid operator * (double z, const fluid &v);
fluid operator * (const fluid &v, double z);
fluid operator + (const fluid &v, const fluid &w );
fluid operator - (const fluid &v, const fluid &w );

#endif // fluid_h





