#ifndef WF_TYPE_
#define WF_TYPE_

#include <complex>
#define complex std::complex<double>

class C_wf
{
public:
  C_wf( complex (*fp_wf_) ( double, int ) = NULL ) : init( fp_wf_ ) {}
  virtual ~C_wf() {}

  complex (*init) ( double, int );
  virtual const char* get_type() { return "wf::base"; }
};

namespace T_wf
{

  class gaussian : public C_wf
  {
  public:
    gaussian() : C_wf( init ) {}
    inline const char* get_type() { return "wf::gaussian"; }
    static complex init( double, int );
  private:
    static double R;
    static double width;
    static double sigma;
    static double alpha;
    static double energy;
  };

  class h1s : public C_wf
  {
  public:
    h1s() : C_wf( init ) {}
    inline const char* get_type() { return "wf::hydrogen_1s"; }
    static complex init( double, int );
  };
}

#endif
