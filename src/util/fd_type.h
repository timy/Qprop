#ifndef FD_TYPE_
#define FD_TYPE_

#include <cstdlib>

class C_field
{
public:
  C_field( double (*field_) (double, int) = NULL ) : field(field_) {}
  virtual ~C_field() {}

  double (*field) ( double, int );
  virtual const char* get_type() { return "field::base"; };
};

namespace T_field
{
  class none : public C_field
  {
  public:
    none() : C_field( field_ ) {}
    inline const char* get_type() { return "field::none"; }
    static double field_( double time, int me );
  };

  class sin2_A : public C_field
  {
  public:
    sin2_A() : C_field( field_ ) {}
    inline const char* get_type() { return "field::vecpot_sin2"; }
    static double field_( double time, int me );
    static double E_0, w, n_c;
  };
}

namespace fd_type
{
  double field( double time, int me );
}

#endif
