#ifndef VP_TYPE_
#define VP_TYPE_

#include <cstdio>
#include <cmath>

class C_vecpot
{
public:
  C_vecpot( double (*vpot_) (double, int) = NULL ) : vpot(vpot_) {}
  virtual ~C_vecpot() {}

  double (*vpot) ( double, int );
  virtual const char* get_type() { return "vecpot::base"; }  
};

namespace T_vecpot
{
  extern double w;
  extern double E_0;
  extern double n_c;
  extern double t_on;
  extern double t_off;
  extern double phi;

  class none : public C_vecpot
  {
  public:
    none() : C_vecpot( vpot_ ) {}
    inline const char* get_type() { return "vecpot::none"; }
    static double vpot_( double time, int me );
  };
  
  class sin2 : public C_vecpot
  {
  public:
    sin2() : C_vecpot( vpot_ ) {}
    inline const char* get_type() { return "vecpot::sin2"; }
    static double vpot_( double time, int me );
  };

  class user : public C_vecpot
  {
  public:
    user( double (*vpot_) ( double, int ) ) : C_vecpot( vpot_ ) {}
    inline const char* get_type() { return "vecpot::user_defined"; }
    //    static double vpot_( double time, int me );
  };
}

#endif
