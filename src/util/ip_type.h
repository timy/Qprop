#ifndef IP_TYPE_
#define IP_TYPE_

#include "grid.h"

namespace ip_type
{
  double ipot( long, long, long, double, grid );
}

class C_imgpot
{
public:
  C_imgpot( double (*ipot_) (long, long, long, double, grid) = NULL ) : ipot(ipot_) {}
  virtual ~C_imgpot() {}

  double (*ipot) ( long, long, long, double, grid );
  virtual const char* get_type() { return "imgpot::base"; }
};

namespace T_imgpot
{
  class none : public C_imgpot
  {
  public:
    none() : C_imgpot( ipot_ ) {}
    inline const char* get_type() { return "imgpot::none"; }
    static double ipot_( long, long, long, double, grid );
  };

  class imgpot : public C_imgpot
  {
  public:
    imgpot() : C_imgpot( ipot_ ) {}
    inline const char* get_type() { return "imgpot::imgpot"; }
    static double ipot_( long, long, long, double, grid );
  };

}



#endif
