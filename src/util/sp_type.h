#ifndef SP_TYPE_
#define SP_TYPE_

#include <cstdlib>
/*
namespace sp_type
{
  namespace Coulomb
  {
    extern double charge;
    double spot( double, double, double, double, int );
  }

  namespace helium_atom
  {
    double spot( double, double, double, double, int );
  }

  namespace none
  {
    double spot( double, double, double, double, int );
  }
}
*/

class C_sclpot
{
public:
  C_sclpot( double (*spot_) (double, double, double, double, int) = NULL ) : spot(spot_) {}
  virtual ~C_sclpot() {}

  double (*spot) ( double, double, double, double, int );
  virtual const char* get_type() { return "sclpot::base"; }  
};

namespace T_sclpot
{
  class none : public C_sclpot
  {
  public:
    none() : C_sclpot( spot_ ), charge( 0.0 ) {}
    inline const char* get_type() { return "sclpot::none"; }
    static double spot_( double, double, double, double, int );
    double charge;
  };
  
  class coulomb : public C_sclpot
  {
  public:
    coulomb() : C_sclpot( spot_ ) {}
    inline const char* get_type() { return "sclpot::Coulomb"; }
    static double spot_( double, double, double, double, int );
    static double charge;
  };

  class helium : public C_sclpot
  {
  public:
    helium() : C_sclpot( spot_ ), charge( 2.0 ) {}
    inline const char* get_type() { return "sclpot::helium"; }
    static double spot_( double, double, double, double, int );
    double charge;
  };
}

#endif
