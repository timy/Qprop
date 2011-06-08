#ifndef wavefunction_h
#define wavefunction_h wavefunction_h
#include<stdio.h>
#include<assert.h>
#include<complex>
#include<iostream>
#include<fstream>
#include<math.h>
#include<cstdlib>
#include<string>

void f(double, double[], double[]);
double CG(double a, double b, double c, double alpha,double beta,double gamma);

#define complex std::complex<double>

class cmatrix;
class fluid;
class hamop;
class grid;

class wavefunction
{
  public:
    wavefunction(long x=0) : wf_dim(x), start(new complex[x]) { }

    wavefunction(const wavefunction& v) 
      {
	wf_dim  = v.wf_dim;
	start = new complex[wf_dim];
	for (long i = 0; i < wf_dim; i++)
	  start[i] = v.start[i];
      }

    virtual ~wavefunction() { delete [] start;}

    //
    // Functions
    //
    void apply(const wavefunction &a, const wavefunction &b, 
	       const wavefunction &c, const wavefunction &psi);
    complex* begin() { return start;}
    const complex* begin() const { return start;}
    
    void calculate_staticpot(grid g, hamop hamil);
    fluid calculate_x_zero(grid g, const fluid &Lambdavector);
    fluid calculate_c_zero(grid g, const fluid &Lambdavector);
    wavefunction calculate_hartree_one(grid g, const wavefunction &Thetavector);
    wavefunction calculate_hartree_two(grid g, const wavefunction &Xivector);
    wavefunction calculate_Lambda_ki(grid g, const fluid &ms);


    fluid calculate_hartree_zero(grid g, const fluid &Lambdavector);    
    fluid calculate_Lambda(grid g, const fluid &degeneracies);
    wavefunction calculate_Xi(grid g, const fluid &degeneracies, const fluid &ms);
    wavefunction calculate_Theta(grid g, const fluid &degeneracies, const fluid &ms);
    fluid calculate_radial_density(grid g, const fluid &degeneracies);
    wavefunction calculate_x_one(grid g, const fluid &Lambdavector, const wavefunction &Thetavector);

    fluid calculate_kli_zero(grid g, const fluid &Lambdavector,
      const fluid &ells, const fluid &ms, const fluid &degeneracies,
      int slateronly, int iv);

    fluid calculate_sic_gam_zero(grid g, const fluid &Lambdavector, 
				 const fluid &hartreezero, 
				 const fluid &ms, 
				 const fluid &degeneracies);
    wavefunction calculate_U_ji_l(grid g, long j, long i, long order);
    wavefunction calculate_rruxsigma(grid g, grid q, const fluid &ells,
      const fluid &degeneracies);
    wavefunction calculate_vsxj(grid q, const fluid &Lambdavector, 
      const wavefunction &rruxsigma, int non_inversion);
    wavefunction calculate_uxi(grid g, grid q, int non_inversion,
      const fluid &ells, const fluid &degeneracies);
    wavefunction calculate_kliav(grid g, grid q, int non_inversion,
      const fluid &klipot, const fluid &degeneracies);
    complex calculate_Mji(grid q, const fluid &Lambdavector, const fluid &ells,
      const fluid &ms, int j, int i);
    wavefunction conjugate();

    void dump_to_file(grid g, FILE* os, int stepwidth, fluid &degeneracies, fluid &ms, double fact);
    void dump_xvector_to_file(long ngpsx,  FILE* os, int stepwidth);
    void dump_to_file_sh(grid g, FILE* os, int stepwidth);
    int  dump_to_file_sh(grid g, FILE* os, int stepwidth, int iv);
    void dump_expect_z(grid g, FILE* os, fluid &degeneracies, const fluid &ms);

    void embed_as_x(grid g, long yindex, long zindex, const wavefunction &v);
    void embed_as_x(grid g, long l, long m, long i, const wavefunction &v);
    void embed_as_xy(grid g, long zindex, const wavefunction &v);
    complex energy(double time, grid g, hamop hamil, int me,
      const wavefunction &staticpot, double charge);

    wavefunction extract_x(grid g, long other_one, long other_two);
    wavefunction extract_y(grid g, long other_one, long other_two);
    wavefunction extract_z(grid g, long other_one, long other_two);
    wavefunction extract_xy(grid g, long zindex);
    wavefunction expwf(complex a);
    wavefunction expwf(double a);
    complex expect_z(grid g, fluid &degeneracies, const fluid &ms);
    complex expect_z(grid g);
    complex expect_cycl_pol_plus(grid g);
    complex expect_cycl_pol_minus(grid g);


    complex expect_rr(grid g, fluid &degeneracies);

    int init(long isize);
    void init(grid g, int inittype, double width, fluid &ells);
    void init( grid g, complex (* ptr_init ) ( double, int ) );
    void init_wf_gaussian( grid g, double R, double width );
    void init_rlm(grid g, int inittype, double w, fluid &ells, fluid &ms);
    void init(grid g, int inittype, FILE* file, int ooi);
    int init(grid g, FILE* file, int ooi, int iv);



    void regrid(grid g, grid g_small, const wavefunction &v);
    void regrid_and_rebin(grid g, grid g_small, const wavefunction &v);

    void select_single_orbital(grid g, grid g_small, int orb_of_interest,
      double ell, const wavefunction &v);

    double totalenergy_x_lda(grid g, const fluid &Lambdavector);
    complex totalenergy_exact_x(grid g, const fluid &ells, const fluid &degeneracies);
    double totalenergy_c_lda(grid g, const fluid &Lambdavector);
    double totalenergy_hartree(grid g, const fluid &Lambdavector, const fluid &hartree_zero);
    double totalenergy_single_part(grid g, const wavefunction &orb_energs, const fluid &degeneracies);
    double totalenergy_sic(grid g, const fluid &degeneracies, const fluid &ms);
    double total_ks_norm(grid g, const fluid &degeneracies);
    double total_ks_norm_in_sphere(grid g, const fluid &degeneracies,
           long xlimit);
    double total_ks_norm_in_shell(grid g, const fluid &degeneracies,
           long xlowerlimit, long xupperlimit);
   
    void nullify();
    double norm(grid g);
    void normalize(grid g, const fluid &ms);    
    void normalize(grid g);    
    wavefunction project(grid g, grid gorig, wavefunction &orig);
    complex project(grid g, grid gorig, wavefunction &orig, long zindex);
    
    int load(FILE*, int);

    complex* end()   { return start + wf_dim;}
    const complex* end()   const { return start + wf_dim;}

    void propagate(complex timestep, 
		   double time, 
		   grid g, 
		   hamop hamil, 
		   int me, 
		   const wavefunction &staticpot, 
		   const fluid &wf_one, 
		   const wavefunction &wf_two, 
		   const wavefunction &wf_three, 
		   const fluid &ms, 
		   double charge,
		   int propornot[]);
    
    void propagate(complex timestep, 
		   double time, 
		   grid g, 
		   hamop hamil, 
		   int me, 
		   const wavefunction &staticpot, 
		   int m, 
		   double charge);
      
    wavefunction orb_fieldenergies(double  time, grid g, hamop hamil, int me, const fluid &ms);

    wavefunction orbital_energies(double time, grid g, hamop hamil, int me,  const wavefunction &staticpot, double charge);
    complex orbital_energy(double time, grid g, hamop hamil, int me,  const wavefunction &staticpot, double charge, long orb_no);


    fluid orbital_norms(grid g);

    wavefunction orbital_hartrees(double time, grid g, int me, 
				  const fluid &wf_one);
    complex orbital_hartree(double time, grid g, int me, 
				  const fluid &wf_one, long orb_no);

    void realific();

    void solve(const wavefunction &a, const wavefunction &b, 
	       const wavefunction &c, const wavefunction &psi, long dimens);
    void solve_du(const wavefunction &a, const wavefunction &b, 
               const wavefunction &c, const wavefunction &psi, long dimens);
    void solve_toep(double a, double b, double b_upperleft, double b_lowerright,
               double c, const wavefunction &psi, long dimens);
    wavefunction sqrtwf();
    wavefunction sqrtrealwf();

    wavefunction mult_diag_with_diag(const wavefunction &a);

    wavefunction invert();
    long  wf_size() const {return wf_dim;}
    //
    // Operators
    //

    wavefunction& operator *= (double z); 
    wavefunction& operator *= (complex z); 
    complex&  operator[](long index)
    {
      //      assert(index >= 0  &&  index < wf_dim);
      return start[index];
    }

    const complex&  operator[](long index) const   
    {   
      //      assert(index >= 0  &&  index < wf_dim);
      return start[index];
    }

    wavefunction& operator=(const wavefunction&);

  private:
    long  wf_dim;
    complex   *start; 

    void do_muller_ell(complex timestep, 
		       double time, 
		       grid g, 
		       hamop hamil, 
		       const wavefunction &staticpot, 
		       int me, 
		       double charge,
		       int m);
      
    void do_muller_ellm(complex timestep, 
				  double time, 
				  grid g, 
				  hamop hamil, 
				  const wavefunction &staticpot, 
				  int me, 
				  double charge);

    void do_muller_general_tddft(complex timestep, 
					   double time, 
					   grid g, 
					   hamop hamil, 
					   const wavefunction &staticpot, 
					   int me, 
					   const fluid &wf_one, 
					   const wavefunction &wf_two,
					   const wavefunction &wf_three,
					   double charge,
					   const fluid &ms,
					   int propornot[]);
};

std::ostream& operator<<( std::ostream& os, const wavefunction& v);
//istream& operator>>(istream& is, wavefunction& v);
complex operator * (const wavefunction &v, const wavefunction &w );
wavefunction operator * (double z, const wavefunction &v);
wavefunction operator * (complex z, const wavefunction &v);
wavefunction operator * (const wavefunction &v, double z);
wavefunction operator / (const wavefunction &v, double z);
wavefunction operator * (const wavefunction &v, complex z);
wavefunction operator + (const wavefunction &v, const wavefunction &vv);
wavefunction operator - (const wavefunction &v, const wavefunction &vv);
wavefunction operator + (const wavefunction &v, const fluid &vv);
wavefunction operator + (const fluid &v, const wavefunction &vv);


#endif // wavefunction_h
