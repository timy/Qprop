#ifndef cmatrix_h
#define cmatrix_h cmatrix_h
#include<assert.h>
#include<complex>
#include<iostream>
#include<fstream>
#include<grid.h>
#include<math.h>
#include<cstdlib>

#define complex std::complex<double>

extern "C" {  
  int dstegr_(char *jobz, char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info );
  int dsyev_(char *jobz, char *uplo, long *n, double *a, long *lda, double *w, double *work, long *lwork, int *info );
  int zgeev_(char *jobvl, char *jobvr, int *n, double *AT, int *lda, double *w, double *vl, int *ldvl, double *vr, int *ldvr, double *WORK, int *lwork, double *RWORK, int *ok);
  int zheev_(char *jobz, char *uplo, long *n, double *a, long *lda, double *w, double *work, long *lwork, double *rwork, int *info);
  int dsytrf_(char *uplo, long *n, double *a, long *lda, long *ipiv, double *work, long *lwork, int *info);
  int dsytrs_(char *uplo, long *n, long *nrhs, double *a, long *lda, long *ipiv, double *b, long *ldb, int *info );
  int dgetrf_(long *m, long *n, double *a, long *lda, long *ipiv, int *info);
 int dgetrs_(char *trans, long *n, long *nrhs, double *a, long *lda, long *ipiv, double *b, long *ldb, int *info );
}


class wavefunction;
class fluid;
class hamop;

class cmatrix
{
  public:
    cmatrix(long x=0, long y=0) : cm__dim(x*y), cm__rows(x), cm__cols(y), start(new complex[x*y]) { }

    cmatrix(const cmatrix& v) 
      {
	cm__dim  = v.cm__dim;
	cm__rows = v.cm__rows;
	cm__cols = v.cm__cols;
	start = new complex[cm__dim];
	for (long i = 0; i < cm__dim; i++)
	  start[i] = v.start[i];
      }

    cmatrix(grid g)
      {
	cm__dim=g.ngps_x()*g.ngps_y();
	cm__rows=g.ngps_x();
	cm__cols=g.ngps_y();
	start = new complex[cm__dim];
      }


    virtual ~cmatrix() { delete [] start;}

    long  cm_size() const {return cm__dim;}
    long  cm_rows() const {return cm__rows;}
    long  cm_cols() const {return cm__cols;}


    cmatrix& operator=(const cmatrix&);

    complex* begin() { return start;}
    complex* end()   { return start + cm__dim;}
    const complex* begin() const { return start;}
    const complex* end()   const { return start + cm__dim;}


    complex elem(long row, long col);
    complex elem(long row, long col) const;

    void set_elem(long row, long col, complex value);
    void set_elem(long row, long col, double value);
    void set_unity();
    void set_zero();
    void embedding_diag(const wavefunction &wf);  
    void ini_tridiag(complex a, complex b, complex c); 

    wavefunction expand_in_colvecs(const wavefunction  &c);
    
    cmatrix transp();
    cmatrix adjoint();
    

    void nr_tqli(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex);
    void nr_tred2_followed_by_tqli(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex);
    void lap_zgeev(wavefunction &eval, cmatrix &levec, cmatrix &revec, long firstindex, long lastindex);
    void lap_zheev(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex);
    void lap_dstegr(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex);
    void lap_dsyev(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex);

    void lap_dsytrf_and_dsytrs(fluid &rhs, int info);
    void lap_dgetrf_and_dgetrs(fluid &rhs, int info);

    void build_atomic_hamiltonian_in_sturmian_basis(
						    grid g, 
						    double alpha, 
						    double theta, 
						    double charge, 
						    cmatrix &cmS, 
						    long reduced_size
						    );

    void build_sturmian_overlap_matrix(grid g);
    void build_z_in_sturmian_basis(grid g);
    void build_pz_in_sturmian_basis(grid g);
    cmatrix build_trafo_sturmian_to_orthogonal(grid g,  wavefunction &evals);
    void build_rdr_matrix_in_sturmian_basis(grid g, cmatrix &cmS, long size);

    cmatrix build_orth_vecpot(cmatrix &cmZorth, wavefunction &Aevals, grid g, cmatrix &cmR, long reduced_size);
    cmatrix build_orth_from_symmetric_matrix(cmatrix &cmX, wavefunction &evals,grid g, cmatrix &cmR, long reduced_size);
    cmatrix build_orth_from_hermitian_matrix(cmatrix &cmX, wavefunction &evals,grid g, cmatrix &cmR, long reduced_size);
    cmatrix build_orth_from_general_matrix(cmatrix &cmLX, cmatrix &cmRX, wavefunction &evals,grid g, cmatrix &cmR, 
					   long reduced_size);

    void build_atomic_propagator_from_orthogonal_hamiltonian(
							     grid g, 
							     complex timestep, 
							     cmatrix &cmY, 
							     cmatrix &cmYinv, 
							     wavefunction &Haevals
							     );
    
    void build_rotation_matrix_from_orthogonal_rdr_in_sturmian_representation(
									      grid g,
									      long size,
									      double theta, 
									      cmatrix &cmY, 
									      cmatrix &cmYinv, 
									      wavefunction &evals,
									      cmatrix &cmS,
									      cmatrix &cmR
									      );

    void build_rotation_matrix_from_orthogonal_rdr_in_orthogonal_representation(
										grid g,
										long size,
										double theta, 
										cmatrix &cmY, 
										cmatrix &cmYinv, 
										wavefunction &evals
										);
    
    void build_interaction_propagator_in_length_gauge(
						      grid g, 
						      complex timestep, 
						      cmatrix &cmY, 
						      cmatrix &cmYinv,
						      wavefunction &evals, 
						      double fieldampl,
						      double alpha,
						      double theta
						      );
    
    void build_interaction_propagator_in_velocity_gauge(
							grid g, 
							complex timestep, 
							cmatrix &cmY, 
							cmatrix &cmYinv,
							wavefunction &evals, 
							double fieldampl,
							double alpha,
							double theta,
							double time
							);
    

    void build_gauge_transformation_operator(
					     grid g, 
					     cmatrix &cmY, 
					     cmatrix &cmYinv,
					     wavefunction &evals, 
					     double fieldampl,
					     double alpha,
					     double theta,
					     double time
					     );
    

    void build_gauge_transform(grid g, cmatrix &cmY, wavefunction &evals, hamop hamil, double alpha, double theta, double time);
    cmatrix build_kinetic_energy_operator(cmatrix &cmYorth, cmatrix &cmYorthinv, wavefunction &Haevals, grid g, double alpha, double theta, cmatrix &cmR, long reduced_size);
    void build_kinetic_orthogonal_propagator(grid g, complex timestep, cmatrix &cmY, cmatrix &cmYinv, wavefunction &Haevals);
    void inverse(cmatrix &cmY, cmatrix &cmYinv, wavefunction &evals);

    cmatrix mult_diag_with_cm(const wavefunction &wf);
    cmatrix mult_cm_with_diag(const wavefunction &wf);

    wavefunction cm_pick_col(long i);
    wavefunction cm_pick_row(long i);
    wavefunction cm_pick_diag();

  private:
    long  cm__dim;
    long  cm__rows,cm__cols;
    complex   *start; 



};


std::ostream& operator<<( std::ostream& os, const cmatrix& v);
cmatrix operator + (const cmatrix &a, const cmatrix &b);
cmatrix operator - (const cmatrix &a, const cmatrix &b);
cmatrix operator * (const cmatrix &a, const cmatrix &b);
cmatrix operator * (const cmatrix &b, complex a);
cmatrix operator * (complex a, const cmatrix &b);
cmatrix operator * (const cmatrix &b, double a);
cmatrix operator * (double a, const cmatrix &b);
cmatrix operator + (const cmatrix &a, const wavefunction &wf);
cmatrix operator + (const wavefunction &wf, const cmatrix &a);
cmatrix operator + (const cmatrix &a, complex b);
cmatrix operator + (complex b, const cmatrix &a);


#endif // cmatrix_h
