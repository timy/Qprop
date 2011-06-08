#include<wavefunction.h>
#include<hamop.h>
#include<fluid.h>
#include<cmatrix.h>
#include<grid.h>
#include<ned.h>
#include<ylm.h>

#define THRESH 1e-6
#define OOS 1.0/6.0
#define TOT 2.0/3.0
#define FOT 5.0/3.0
#define SQRTT sqrt(2.0)
#define complex std::complex<double>


double CG(double a, double b, double c, double alpha, double beta, double gamma)
{
  double result;
  ned_(&a,&b,&c,&alpha,&beta,&gamma,&result);
  return result;
}

void wavefunction::nullify()
{
  for (long i = 0; i < wf_dim; i++)
    start[i] = complex(0.0,0.0); 

}


/*! \fn int wavefunction::init(long isize)
  allocates array of \a isize complex numbers.
  
  \param  isize number of elements to allocate.
  
  \return a negative integer if an error occurs, null otherwise.

  \remark The function sets to all the allocated elements the zero
  initial value. If the wavefunction was already allocated, then
  the array will be first reallocated (deleted and then new created).
*/
int wavefunction::init(long lsize)
{
  wf_dim = lsize;
  start = new complex[lsize];
  for (long i = 0; i < lsize; i++) {start[i] = complex(0.0,0.0);}
  return 0;
}


/*! \fn int wavefunction::init(grid g, int inittype, double width,
  fluid &ells)

  initializes the wavefunction with an wf from the file.  
  \param  g.
  \param inittype;
  \param width
  \param ells
  
  \return a negative integer if an error occurs, null otherwise.

  \remark The function sets to all the allocated elements the zero
  initial value. If the wavefunction was already allocated, then
  the array will be first reallocated (deleted and then new created).
*/

void wavefunction::init(grid g, int inittype, double w, fluid &ells)
{
  long xindex, yindex, zindex, index_i,index,i;
  double x,y,z,r;
  double offs = 0.0;
  double realpart, imagpart;
   
  if (g.dimens() != 34 )   
    {    
      fprintf(stderr, "init: you are not in the proper propagation mode to initialize using this function!\n");
      exit(-50);
    };
  
  srand((int)(w)+1);

  // xindex runs over the radial grid points
  // yindex runs over the l quantum numbers
  // zindex runs over the orbitals
   
  for (xindex=0; xindex<g.ngps_x(); xindex++) { 
    r=g.r(xindex);      
    for (yindex=0; yindex<g.ngps_y(); yindex++) {
      y=g.y(yindex);
      for (zindex=0; zindex<g.ngps_z(); zindex++) {
	z=g.z(zindex);
	index_i=g.index(xindex,yindex,zindex);
	start[index_i] = complex(0.0, 0.0);
	switch (inittype) {	
	case 1 :
	  if (ells[zindex]==yindex)
	    //std::cout << xindex << '\t' <<  yindex << '\t' << zindex << std::endl; 
	    start[index_i]=complex(-rand()/(RAND_MAX+1.0)+0.5,0.0);
	  else
	    start[index_i]=complex(0.0,0.0);
	  break;
		  
	  //
	  // hydrogen-like orbitals
	  //
	case 2 :
	      
	  switch(zindex) {	
	  case 0: 
	    if ( ells[zindex] == yindex )
	      //std::cout << xindex << '\t' <<  yindex << '\t' << zindex << std::endl; 
	      start[index_i]=r*exp(- w * r );
	    else
	      start[index_i]=complex(0.0,0.0);
	    break;
		
	  case 1: 
	    if ( ells[zindex] == yindex )
	      //std::cout << xindex << '\t' <<  yindex << '\t' << zindex << std::endl; 
	      start[index_i]=r*exp(-0.5*w*r);
	    else
	      start[index_i]=complex(0.0,0.0);
	    break;
		      
	  case 2: 
	    if (ells[zindex]==yindex)
	      //std::cout << xindex << '\t' <<  yindex << '\t' << zindex << std::endl; 
	      start[index_i]=r*r*exp(-0.5*w*r);
	    else
	      start[index_i]=complex(0.0,0.0);
	    break;
		  
	  case 3: 
	    if (ells[zindex]==yindex)
	      //std::cout << xindex << '\t' <<  yindex << '\t' << zindex << std::endl; 
	      start[index_i]=r*r*exp(-0.5*w*r);
	    else
	      start[index_i]=complex(0.0,0.0);
	    break;
		  
	  case 4: 
	    if (ells[zindex]==yindex)
	      //std::cout << xindex << '\t' <<  yindex << '\t' << zindex << std::endl; 
	      start[index_i]=r*r*exp(-0.5*w*r);
	    else
	      start[index_i]=complex(0.0,0.0);
	    break;

	  }
	  break; 
	default :
	  start[index_i]=complex(0.0,0.0);
	}
      }
    }
  }
  return;
}

// put a specified wave function on the grid
void wavefunction::init( grid g, complex (* ptr_init ) ( double, int ) )
{
  long i_r, i_l, index;
  double r;
  std::cout << "Initialization with user-defined wave function" << std::endl;

  for( index = 0; index < g.size(); index ++ ) {
    start[ index ] = complex( 0.0, 0.0 );
  }
  for( i_r = 0; i_r < g.ngps_x(); i_r ++ ) {
    r = g.r( i_r );
    for( i_l = 0; i_l < g.ngps_y(); i_l ++ ) {
      index = g.index( i_r, i_l, 0 );
      start[ index ] = ptr_init( r, i_l );
    }
  }
  return;
}

/*!

 */


void wavefunction::init_rlm(grid g, int inittype, double w, fluid &ells, fluid &ms)
{
  long xindex, yindex, zindex, index_i,index,i,m;
  double x,y,z;
  double offs=0.0;
  double realpart, imagpart;
  long ctrl_id;

  zindex=0;  
  srand((int)(w)+1);

  if (g.dimens() != 44 )   
    {    
      fprintf(stderr, "init_rlm: you are not in rlm-mode!\n");
      exit(-50);
    };
  
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    { 
      x=g.r(xindex);
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  for (m=-yindex; m<=yindex; m++)
	    {
	      index=g.rlmindex(xindex,yindex,m);
	      switch (inittype)
		{
		case 1 :
		  if ((yindex==(int)(ells[zindex])) && (m==(int)(ms[zindex])))
		    {
		      start[index]=complex(10.0*rand()/(RAND_MAX+1.0)-0.5,0.0);
		    }
		  else
		    {
		      start[index]=complex(0.0,0.0);
		    }
		  break;
		case 2 :
		  if ((yindex==(int)(ells[zindex])) && (m==(int)(ms[zindex])))
		    {
		      start[index]=g.r(xindex)*exp(-g.r(xindex));
		    }
		  else
		    {
		      start[index]=complex(0.0,0.0); 
		    };
		  break;
		default :
		  start[index]=complex(0.0,0.0);
		};
	    };
	};
      
    };
      
}



void wavefunction::init(grid g, int inittype, FILE* file,
			int output_of_interest)

{
  init(g, file,output_of_interest, 1);
}


/*!

 */


int wavefunction::init(grid g, FILE* file, int output_of_interest, int iv)
{
  long index;
  double realpart, imagpart;
  long ctrl_id,i;
  
  for (i=0; i<output_of_interest; i++) 
    {
      if (iv==1) fprintf(stdout, "Ignoring of %ld\n", i);
      for (index=0; index<g.size(); index++)
	{
	  ctrl_id=fscanf(file,"%lf %lf",&realpart,&imagpart);
	}
    }

  if (iv==1) fprintf(stdout, "Now I am storing\n");
  for (index=0; index<g.size(); index++)
    {
      ctrl_id=fscanf(file,"%lf %lf",&realpart,&imagpart);
      start[index]=complex(realpart,imagpart);
    }
  return 0;
}


int wavefunction::load(FILE* os, int iv)
{
  double re, im;
  long ctrl_id;

  for (long i=0; i<wf_size(); i++)
    {
      ctrl_id=fscanf(os,"%lf %lf",&re ,&im);
      start[i]=complex(re, im);
    }
  return 0;
}


wavefunction &wavefunction::operator=(const wavefunction &v)
{
  if (this != &v)
    {
      delete[] start;
      wf_dim  = v.wf_dim;
      start = new complex[wf_dim];
      for (long i = 0; i < wf_dim; i++)
	start[i] = v.start[i];
    }
  return *this;
}

wavefunction& wavefunction::operator *= (double z)
{
  for(long i=0; i<wf_dim; i++)
    start[i]=start[i]*(complex)z;
  return *this;
}

wavefunction& wavefunction::operator *= (complex z)
{
  for(long i=0; i<wf_dim; i++)
    start[i]=start[i]*z;
  return *this;
}

wavefunction operator + (const wavefunction &v, const wavefunction &vv)
{
  wavefunction temp=v;
  for(long i=0; i<v.wf_size(); i++)
    temp[i]=v[i]+vv[i];
  return temp;
}

wavefunction operator - (const wavefunction &v, const wavefunction &vv)
{
  wavefunction temp=v;
  for(long i=0; i<v.wf_size(); i++)
    temp[i]=v[i]-vv[i];
  return temp;
}

wavefunction operator + (const fluid &v, const wavefunction &vv)
{
  wavefunction temp=vv;
  for(long i=0; i<vv.wf_size(); i++)
    temp[i]=complex(v[i],0.0)+vv[i];
  return temp;
}

wavefunction operator + (const wavefunction &v, const fluid &vv)
{
  wavefunction temp=v;
  for(long i=0; i<v.wf_size(); i++)
    temp[i]=v[i]+complex(vv[i],0.0);
  return temp;
}


wavefunction operator * (double z, const wavefunction &v)
{
  wavefunction temp=v;
  return temp *= z;
}

wavefunction operator * (complex z, const wavefunction &v)
{
  wavefunction temp=v;
  return temp *= z;
}

wavefunction operator * (const wavefunction &v, double z)
{
  wavefunction temp=v;
  return temp *= z;
}

wavefunction operator / (const wavefunction &v, double z)
{
  wavefunction temp=v;
  return temp *= 1.0/z;
}

wavefunction operator * (const wavefunction &v, complex z)
{
  wavefunction temp=v;
  return temp *= z;
}

complex operator * (const wavefunction &v, const wavefunction &w )
{
  complex result(0.0,0.0);
  for(long i=0; i<v.wf_size(); i++)
    {
      result+=conj(v[i])*w[i];
    }
  return result;
}




std::ostream& operator<<( std::ostream& os, const wavefunction& v)
{  
  for(long i = 0; i < v.wf_size(); i++)
    {   os << real(v[i]) << " " << imag(v[i]) << std::endl;
    }
  return os;
}

void wavefunction::regrid(grid g, grid g_small, const wavefunction &v)
{
  long xindex,yindex,zindex,index,index_small, 
    xupper, yupper, zupper, m_limit, mindex;

  if (g_small.ngps_x() > g.ngps_x()){xupper=g.ngps_x();}
  else {xupper=g_small.ngps_x();}

  if (g_small.ngps_z() > g.ngps_z()) {zupper=g.ngps_z();}
  else {zupper=g_small.ngps_z();}

  if (g_small.ngps_y() > g.ngps_y()) {yupper=g.ngps_y();}
  else {yupper=g_small.ngps_y();}
  
  for (xindex=0; xindex<xupper; xindex++)
    {
      for (yindex=0; yindex<yupper; yindex++)
	{
	  if (g.dimens()==34) {m_limit=0;} else {m_limit=yindex;}
	  for (mindex=-m_limit; mindex<=m_limit; mindex++)
	    {
	      for (zindex=0; zindex<zupper; zindex++)
		{
		  index=g.index(xindex,yindex,mindex,zindex);
		  index_small=g_small.index(xindex,yindex,mindex,zindex);
		  start[index]=v[index_small];
		}
	    }
	}
    }
}


void wavefunction::regrid_and_rebin(grid g, grid g_small, const wavefunction &v)
{
  long xindex,yindex,zindex,index,index_small,xindex_small;
  double r;
  long yupper, zupper;

  if (g_small.ngps_z() > g.ngps_z())
    {
      zupper=g.ngps_z();
    }
  else
    {
      zupper=g_small.ngps_z();
    };
  if (g_small.ngps_y() > g.ngps_y())
    {
      yupper=g.ngps_y();
    }
  else
    {
      yupper=g_small.ngps_y();
    };
  
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      r=g.r(xindex);
      xindex_small=g_small.rindex(r);
      if (xindex_small<g_small.ngps_x())
	{
	  for (yindex=0; yindex<yupper; yindex++)
	    {
	      for (zindex=0; zindex<zupper; zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  index_small=g_small.index(xindex_small,yindex,zindex);
		  start[index]=v[index_small];
		};
	    };
	};
    };

}

void wavefunction::select_single_orbital(grid g, grid g_small, int orb_of_interest, double ell, const wavefunction &v)
{
  long xindex,yindex,zindex,index,index_small,xindex_small;
  double r;
  
  yindex=(long)(ell);
  zindex=orb_of_interest;

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      r=g.r(xindex);
      xindex_small=g_small.rindex(r);
      
      if (xindex_small<g_small.ngps_x())
	{
	  index=g.index(xindex,yindex,0);
	  index_small=g_small.index(xindex_small,yindex,zindex);
	  start[index]=v[index_small];
	};
    };

}


void wavefunction::embed_as_xy(grid g, long zindex, const wavefunction &v)
{
  long xindex,yindex, index, iindex;
  grid gg;
  gg.set_dim(34);
  gg.set_ngps(g.ngps_x(),g.ngps_y(),1);
  gg.set_delt(g.delt_x(),1,1);
  gg.set_offs(0,0,0);

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  iindex=gg.index(xindex,yindex,0);
	  start[index]=v[iindex];
	};
    };
  
}



void wavefunction::embed_as_x(grid g, long yindex, 
			      long zindex, const wavefunction &v)
{
  long i,index;

  for (i=0; i<v.wf_size(); i++)
    {
      index=g.index(i,yindex,zindex);
      start[index]=v[i];
    }
}


void wavefunction::embed_as_x(grid g,long l,long m,long i,const wavefunction &v)
{
  for (long r=0; r<v.wf_size(); r++) start[g.index(r,l,m,i)]=v[r];
}



wavefunction wavefunction::expwf(complex a)
{
  long i;
  wavefunction result(wf_dim);

  for (i=0; i<wf_dim; i++)
    {
      result[i]=exp(a*start[i]);
    };
  
  return result;
  
}  

wavefunction wavefunction::conjugate()
{
  long i;
  wavefunction result(wf_dim);

  for (i=0; i<wf_dim; i++)
    {
      result[i]=conj(start[i]);
    };
  return result;
} 


wavefunction wavefunction::sqrtwf()
{
  long i;
  wavefunction result(wf_dim);

  for (i=0; i<wf_dim; i++)
    {
      result[i]=sqrt(start[i]);
    };
  
  return result;
  
}  

wavefunction wavefunction::sqrtrealwf()
{
  long i;
  wavefunction result(wf_dim);

  for (i=0; i<wf_dim; i++)
    {
      result[i]=sqrt(real(start[i]));
    };
  
  return result;
  
}  

wavefunction wavefunction::expwf(double a)
{
  long i;
  wavefunction result(wf_dim);

  for (i=0; i<wf_dim; i++)
    {
      result[i]=exp((complex)a*start[i]);
    };
  
  return result;
  
}  

wavefunction wavefunction::mult_diag_with_diag(const wavefunction &a)
{
  long i,smaller;
   
  if (wf_dim<a.wf_size()) 
    smaller=wf_dim;
  else 
    smaller=a.wf_size();
   
  wavefunction result(smaller);
   
  for (i=0; i<smaller; i++)
    {
      result[i]=start[i]*a[i];
    };

  return result;
    
}


wavefunction wavefunction::invert()
{
  long i;
  wavefunction result(wf_dim);

  for (i=0; i<wf_dim; i++)
    {
      result[i]=1.0/start[i];
    };
  
  return result;
  
}  


wavefunction wavefunction::extract_x(grid g, long other_one, long other_two)
{
  long i, index;   
  wavefunction result(g.ngps_x());
   
  for (i=0; i<g.ngps_x(); i++)
    {
      index=g.index(i,other_one,other_two);
      result[i]=start[index];
    };

  return result;
}

wavefunction wavefunction::extract_y(grid g, long other_one, long other_two)
{
  long i, index;
  
  wavefunction result(g.ngps_y());
  for (i=0; i<g.ngps_y(); i++)
    {
      index=g.index(other_one,i,other_two);
      result[i]=start[index];
    };

  return result;
}

wavefunction wavefunction::extract_z(grid g, long other_one, long other_two)
{
  long i, index;
  
  wavefunction result(g.ngps_z());
  for (i=0; i<g.ngps_z(); i++)
    {
      index=g.index(other_one,other_two,i);
      result[i]=start[index];
    };

  return result;
}

wavefunction wavefunction::extract_xy(grid g, long zindex)
{
  long xindex,yindex, index, iindex;

  grid gg;
  gg.set_dim(34);
  gg.set_ngps(g.ngps_x(),g.ngps_y(),1);
  gg.set_delt(g.delt_x(),1,1);
  gg.set_offs(0,0,0);

  wavefunction result(g.ngps_x()*g.ngps_y());
  
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  iindex=gg.index(xindex,yindex,0);
	  result[iindex]=start[index];
	}
    }  
  return result;
}





// this is for several KS orbitals; overloaded below
void wavefunction::propagate(complex timestep, 
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
			     int propornot[])
{

  switch (g.dimens())
    {

    case 34 :

      do_muller_general_tddft(timestep,time,g,hamil,staticpot,me,
			      wf_one,wf_two,wf_three,charge,ms,propornot);
      
      break;


    case 44 :

      do_muller_ellm(timestep,time,g,hamil,staticpot,me,charge);
      
      break;

    default :
      fprintf(stderr, "Unknown propagation mode (g.dimens())\n");
      exit(-50);

    };

}

// this is for a single orbital with V_H and V_xc zero; overloaded above
void wavefunction::propagate(complex timestep, 
			     double time, 
			     grid g, 
			     hamop hamil, 
			     int me, 
			     const wavefunction &staticpot, 
			     int m, 
			     double charge)
{

  switch (g.dimens())
    {

    case 34 :

      do_muller_ell(timestep,time,g,hamil,staticpot,me,charge,m);
      
      break;


    case 44 :

      do_muller_ellm(timestep,time,g,hamil,staticpot,me,charge);
      
      break;

    default :
      fprintf(stderr, "Unknown propagation mode (g.dimens())\n");
      exit(-50);

    };

}




void wavefunction::do_muller_ellm(complex timestep, 
				  double time, 
				  grid g, 
				  hamop hamil, 
				  const wavefunction &staticpot, 
				  int me,
				  double charge)
{

  double Eff_Charge = charge;

  long xindex, yindex, m, ell;
  long index, index_i, index_xp, index_xm, index_lp;
  double r;
  wavefunction rhsone(g.size());
  double wfsq;
  wavefunction aa(g.ngps_x());
  wavefunction bb(g.ngps_x());
  wavefunction cc(g.ngps_x());	
  wavefunction tmp_one(g.ngps_x());	
  wavefunction tmp_two(g.ngps_x());	
  wavefunction tmptmp_one(g.ngps_x());	
  wavefunction tmptmp_two(g.ngps_x());	
  complex imagi(0.0,1.0);
  complex A;
  double phase;
  double Amag;
  complex expp;
  complex expm;
  double blm,btildelm,dlm,dtildelm;
  complex Delta_two_upperleft, M_two_upperleft;
  double x=sqrt(3.0)-2.0;
  double cornerpart=OOS*(4.0+x);
  double oosqrttwo=1.0/sqrt(2.0);
  complex halfimagitimestep=(complex)0.5*timestep*imagi;
  complex halfimagitimestepOOS=halfimagitimestep*(complex)OOS;
  complex halfimagitimestepFOT=halfimagitimestep*(complex)FOT;
  complex aaaa, bbbb, cccc;
  double aaa,bbb,ccc,b_upperleft,b_lowerright,zeta,xi,prefaci;

  complex factor,ur,ll,diagonal;
  



  //========================== the right hand side ======================

  A=hamil.vecpot_x(time,me)+imagi*hamil.vecpot_y(time,me);
  phase=arg(A);
  Amag=abs(A);
  expp=exp(imagi*phase);
  expm=exp(-imagi*phase);
  prefaci=real(timestep)*Amag/8.0;

  zeta=real(timestep)*Amag/(16.0*g.delt_x());

  for (ell=0; ell<g.ngps_y()-1; ell++)
    {
      for (m=-ell; m<=ell; m++)
	{
	  
	  blm=sqrt((ell-m+1.0)/((2.0*ell+1)*(2.0*ell+3)))
	    *(m*sqrt(ell-m+2)-sqrt((ell+m+1)*((ell+1)*(ell+2)-m*(m-1))));
	  
	  // ------------------ R -------------------------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      xi=prefaci/g.r(xindex);
	      factor=1.0/(1.0+xi*xi*blm*blm);
	      diagonal=(1.0-xi*xi*blm*blm)*factor;
	      ur=2.0*xi*expm*blm*factor;
	      ll=-2.0*xi*expp*blm*factor;
	      index=g.	rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m-1);
	      rhsone[index]=diagonal*start[index]+ur*start[index_i];
	      rhsone[index_i]=diagonal*start[index_i]+ll*start[index];
	    };
	  // ------- end of --- R -------------------------

	  // ------- apply B and construct Gamma_1 and Gamma_2 -----------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m-1);
	      tmp_one[xindex]=oosqrttwo*(-expp*rhsone[index]+rhsone[index_i]);
	      tmp_two[xindex]=oosqrttwo*(expp*rhsone[index]+rhsone[index_i]);
	    };
	  // ------- end of applying B etc.  --------------


	  dlm=sqrt((ell-m+1.0)*(ell-m+2.0)/((2.0*ell+1.0)*(2.0*ell+3.0)));

	  // ------- calculate Yoverlminusone times Gamma_1
	  xindex=0;
	  tmptmp_one[xindex]=(cornerpart+x*zeta*dlm)*tmp_one[xindex]
	    +(OOS+zeta*dlm)*tmp_one[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_one[xindex]=(OOS-zeta*dlm)*tmp_one[xindex-1]
		+TOT*tmp_one[xindex]
		+(OOS+zeta*dlm)*tmp_one[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_one[xindex]=(OOS-zeta*dlm)*tmp_one[xindex-1]
	    +(cornerpart-x*zeta*dlm)*tmp_one[xindex];
	  // -------- end of calculating Yoverlminusone times Gamma_1


	  // ------- calculate Yoverlminustwo times Gamma_2
	  xindex=0;
	  tmptmp_two[xindex]=(cornerpart-x*zeta*dlm)*tmp_two[xindex]
	    +(OOS-zeta*dlm)*tmp_two[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_two[xindex]=(OOS+zeta*dlm)*tmp_two[xindex-1]
		+TOT*tmp_two[xindex]
		+(OOS-zeta*dlm)*tmp_two[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_two[xindex]=(OOS+zeta*dlm)*tmp_two[xindex-1]
	    +(cornerpart+x*zeta*dlm)*tmp_two[xindex];
	  // -------- end of calculating Yoverlminustwo times Gamma_2


	  // -------- apply Yoverlplusone -----------------------------
	  b_upperleft=cornerpart-x*zeta*dlm;
	  b_lowerright=cornerpart+x*zeta*dlm;
	  aaa=OOS-zeta*dlm;
	  bbb=TOT;
	  ccc=OOS+zeta*dlm;

	  tmp_one.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_one,g.ngps_x());
	  // -------- end of applying  Yoverlplusone ------------------


	  // -------- apply Yoverlplustwo -----------------------------
	  b_upperleft=cornerpart+x*zeta*dlm;
	  b_lowerright=cornerpart-x*zeta*dlm;
	  aaa=OOS+zeta*dlm;
	  bbb=TOT;
	  ccc=OOS-zeta*dlm;

	  tmp_two.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_two,g.ngps_x());
	  // -------- end of applying  Yoverlplustwo ------------------

	  // ------- apply B^+ and assemble the stuff -----------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m-1);
	      start[index]=oosqrttwo*(-expm*tmp_one[xindex]+expm*tmp_two[xindex]);
	      start[index_i]=oosqrttwo*(tmp_one[xindex]+tmp_two[xindex]);
	    };
	  // ------- end of applying B^+ etc.  --------------


	  // ====================== and now the same stuff with tilde =======================
	  // ================= attention! is a different lm-subblock now! ===================

	  btildelm=sqrt((ell+m+1.0)/((2.0*ell+1)*(2.0*ell+3)))
	    *(m*sqrt(ell+m+2)+sqrt((ell-m+1)*((ell+1)*(ell+2)-m*(m+1))));

	  // ------------------ Rtilde -------------------------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      xi=prefaci/g.r(xindex);
	      factor=1.0/(1.0+xi*xi*btildelm*btildelm);
	      diagonal=(1.0-xi*xi*btildelm*btildelm)*factor;
	      ur=2.0*xi*expp*btildelm*factor;
	      ll=-2.0*xi*expm*btildelm*factor;
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m+1);
	      rhsone[index]=diagonal*start[index]+ur*start[index_i];
	      rhsone[index_i]=diagonal*start[index_i]+ll*start[index];
	    };
	  // ------- end of --- Rtilde  -------------------------


	  // ------- apply Btilde and construct Gamma_1tilde and Gamma_2tilde ---
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m+1);
	      tmp_one[xindex]=oosqrttwo*(-expm*rhsone[index]+rhsone[index_i]);
	      tmp_two[xindex]=oosqrttwo*(expm*rhsone[index]+rhsone[index_i]);
	    };
	  // ------- end of applying Btilde etc.  --------------


	  dtildelm=sqrt((ell+m+1.0)*(ell+m+2.0)/((2.0*ell+1.0)*(2.0*ell+3.0)));

	  // ------- calculate Yoverlminustildeone times Gamma_1tilde
	  xindex=0;
	  tmptmp_one[xindex]=(cornerpart-x*zeta*dtildelm)*tmp_one[xindex]
	    +(OOS-zeta*dtildelm)*tmp_one[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_one[xindex]=(OOS+zeta*dtildelm)*tmp_one[xindex-1]
		+TOT*tmp_one[xindex]
		+(OOS-zeta*dtildelm)*tmp_one[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_one[xindex]=(OOS+zeta*dtildelm)*tmp_one[xindex-1]
	    +(cornerpart+x*zeta*dtildelm)*tmp_one[xindex];
	  // -------- end of calculating Yoverlminustildeone times Gamma_1tilde


	  // ------- calculate Yoverlminustildetwo times Gamma_2tilde
	  xindex=0;
	  tmptmp_two[xindex]=(cornerpart+x*zeta*dtildelm)*tmp_two[xindex]
	    +(OOS+zeta*dtildelm)*tmp_two[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_two[xindex]=(OOS-zeta*dtildelm)*tmp_two[xindex-1]
		+TOT*tmp_two[xindex]
		+(OOS+zeta*dtildelm)*tmp_two[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_two[xindex]=(OOS-zeta*dtildelm)*tmp_two[xindex-1]
	    +(cornerpart-x*zeta*dtildelm)*tmp_two[xindex];
	  // -------- end of calculating Yoverlminustildetwo times Gamma_2tilde


	  // -------- apply Yoverlplustildeone -----------------------------
	  b_upperleft=cornerpart+x*zeta*dtildelm;
	  b_lowerright=cornerpart-x*zeta*dtildelm;
	  aaa=OOS+zeta*dtildelm;
	  bbb=TOT;
	  ccc=OOS-zeta*dtildelm;

	  tmp_one.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_one,g.ngps_x());
	  // -------- end of applying  Yoverlplustildeone ------------------


	  // -------- apply Yoverlplustildetwo -----------------------------
	  b_upperleft=cornerpart-x*zeta*dtildelm;
	  b_lowerright=cornerpart+x*zeta*dtildelm;
	  aaa=OOS-zeta*dtildelm;
	  bbb=TOT;
	  ccc=OOS+zeta*dtildelm;

	  tmp_two.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_two,g.ngps_x());
	  // -------- end of applying  Yoverlplustildetwo ------------------

	  // ------- apply Btilde^+ and assemble the stuff -----------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m+1);
	      start[index]=oosqrttwo*(-expp*tmp_one[xindex]+expp*tmp_two[xindex]);
	      start[index_i]=oosqrttwo*(tmp_one[xindex]+tmp_two[xindex]);
	    };
	  // ------- end of applying Btilde^+ etc.  --------------

	};
    };


  //====== end of ============ the right hand side ======================


  // the spatial part
  // The constant part of the matrix L_-(tau)
  aaaa=-OOS-halfimagitimestep/(g.delt_x()*g.delt_x());
  cccc=aaaa;
  bbbb=-FOT+imagi*timestep/(g.delt_x()*g.delt_x());
  Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
  M_two_upperleft=-(complex)2.0*(1.0+(complex)g.delt_x()*(complex)g.delt_x()/12.0*Delta_two_upperleft);
  
  // Calculate the rhs vector W_-(tau) *this
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    { 
      for (m=-yindex; m<=yindex; m++)
	{
	  xindex=0;
	  index=g.rlmindex(xindex,yindex,m);
	  index_xp=g.rlmindex(xindex+1,yindex,m);
	  if (yindex==0)
	    {
	      rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]))*start[index_xp]
		+((complex)M_two_upperleft*(1.0-halfimagitimestep*(staticpot[index]))-imagi*(complex)Delta_two_upperleft*(complex)0.5*timestep)*start[index];
	    }
	  else
	    {
	      rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]))*start[index_xp]
		+(bbbb+halfimagitimestepFOT*(staticpot[index]))*start[index];
	    };
	      
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    { 
	      index=g.rlmindex(xindex,yindex,m);
	      index_xp=g.rlmindex(xindex+1,yindex,m);
	      index_xm=g.rlmindex(xindex-1,yindex,m);
	      rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]))*start[index_xp]
		+(bbbb+halfimagitimestepFOT*(staticpot[index]))*start[index]
		+(cccc+halfimagitimestepOOS*(staticpot[index_xm]))*start[index_xm];
	    };
	      
	  
	  xindex=g.ngps_x()-1;
	  index=g.rlmindex(xindex,yindex,m);
	  index_xm=g.rlmindex(xindex-1,yindex,m);
	  rhsone[index]=(bbbb+halfimagitimestepFOT*(staticpot[index]))*start[index]
	    +(cccc+halfimagitimestepOOS*(staticpot[index_xm]))*start[index_xm];
	};
    };
	  
  // The matrix W_+(tau)
  aaaa=-OOS+halfimagitimestep/(g.delt_x()*g.delt_x());
  cccc=aaaa;
  bbbb=-FOT-imagi*timestep/(g.delt_x()*g.delt_x());
	  
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    { 
      for (m=-yindex; m<=yindex; m++)
	{
	  xindex=0;
	  index=g.rlmindex(xindex,yindex,m);
	  index_xp=g.rlmindex(xindex+1,yindex,m);
	  aa[xindex]=(aaaa-halfimagitimestep*(complex)OOS*(staticpot[index_xp]));
	  cc[xindex]=1.0; // not used
	      
	  if (yindex==0)
	    {
	      bb[xindex]=((complex)M_two_upperleft*(1.0+halfimagitimestep*(staticpot[index]))+imagi*(complex)Delta_two_upperleft*(complex)0.5*timestep);
	    }
	  else
	    {
	      bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]));
	    };
	      
	  tmp_one[xindex]=rhsone[index];
	  
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    { 
	      index=g.rlmindex(xindex,yindex,m);
	      index_xp=g.rlmindex(xindex+1,yindex,m);
	      index_xm=g.rlmindex(xindex-1,yindex,m);
		  
	      aa[xindex]=(aaaa-halfimagitimestepOOS*(staticpot[index_xp]));
	      cc[xindex]=(cccc-halfimagitimestepOOS*(staticpot[index_xm]));
	      bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]));
	      tmp_one[xindex]=rhsone[index];
	    };  
	      
	  xindex=g.ngps_x()-1;
	  r=g.r(xindex);
	  index=g.rlmindex(xindex,yindex,m);
	  index_xm=g.rlmindex(xindex-1,yindex,m);
	  aa[xindex]=aaaa-halfimagitimestepOOS*(hamil.scalarpotx(r+g.delt_x(),yindex,0.0,time,me)
						-imagi*(complex)hamil.imagpot(xindex,yindex,0,time,g));
	  cc[xindex]=(cccc-halfimagitimestepOOS*(staticpot[index_xm]));
	  bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]));
	  tmp_one[xindex]=rhsone[index];
	      
	  tmp_two.solve(aa,bb,cc,tmp_one,g.ngps_x());
	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      index=g.rlmindex(xindex,yindex,m);
	      start[index]=tmp_two[xindex];
	    };  
	};
    };
	  

  //========================== the left hand side ======================

  time=time+0.5*real(timestep);
  A=hamil.vecpot_x(time,me)+imagi*hamil.vecpot_y(time,me);
  phase=arg(A);
  Amag=abs(A);
  expp=exp(imagi*phase);
  expm=exp(-imagi*phase);
  prefaci=real(timestep)*Amag/8.0;


  zeta=real(timestep)*Amag/(16.0*g.delt_x());

  for (ell=g.ngps_y()-2; ell>=0; ell--)
    {
      for (m=ell; m>=-ell; m--)
	{
	  

	  // ====================== now first the stuff with tilde =======================

	  btildelm=sqrt((ell+m+1.0)/((2.0*ell+1)*(2.0*ell+3)))
	    *(m*sqrt(ell+m+2)+sqrt((ell-m+1)*((ell+1)*(ell+2)-m*(m+1))));


	  // ------- apply Btilde and construct Gamma_1tilde and Gamma_2tilde ---
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m+1);
	      tmp_one[xindex]=oosqrttwo*(-expm*start[index]+start[index_i]);
	      tmp_two[xindex]=oosqrttwo*(expm*start[index]+start[index_i]);
	    };
	  // ------- end of applying Btilde etc.  --------------


	  dtildelm=sqrt((ell+m+1.0)*(ell+m+2.0)/((2.0*ell+1.0)*(2.0*ell+3.0)));

	  // ------- calculate Yoverlminustildeone times Gamma_1tilde
	  xindex=0;
	  tmptmp_one[xindex]=(cornerpart-x*zeta*dtildelm)*tmp_one[xindex]
	    +(OOS-zeta*dtildelm)*tmp_one[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_one[xindex]=(OOS+zeta*dtildelm)*tmp_one[xindex-1]
		+TOT*tmp_one[xindex]
		+(OOS-zeta*dtildelm)*tmp_one[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_one[xindex]=(OOS+zeta*dtildelm)*tmp_one[xindex-1]
	    +(cornerpart+x*zeta*dtildelm)*tmp_one[xindex];
	  // -------- end of calculating Yoverlminustildeone times Gamma_1tilde


	  // ------- calculate Yoverlminustildetwo times Gamma_2tilde
	  xindex=0;
	  tmptmp_two[xindex]=(cornerpart+x*zeta*dtildelm)*tmp_two[xindex]
	    +(OOS+zeta*dtildelm)*tmp_two[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_two[xindex]=(OOS-zeta*dtildelm)*tmp_two[xindex-1]
		+TOT*tmp_two[xindex]
		+(OOS+zeta*dtildelm)*tmp_two[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_two[xindex]=(OOS-zeta*dtildelm)*tmp_two[xindex-1]
	    +(cornerpart-x*zeta*dtildelm)*tmp_two[xindex];
	  // -------- end of calculating Yoverlminustildetwo times Gamma_2tilde


	  // -------- apply Yoverlplustildeone -----------------------------
	  b_upperleft=cornerpart+x*zeta*dtildelm;
	  b_lowerright=cornerpart-x*zeta*dtildelm;
	  aaa=OOS+zeta*dtildelm;
	  bbb=TOT;
	  ccc=OOS-zeta*dtildelm;

	  tmp_one.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_one,g.ngps_x());
	  // -------- end of applying  Yoverlplustildeone ------------------


	  // -------- apply Yoverlplustildetwo -----------------------------
	  b_upperleft=cornerpart-x*zeta*dtildelm;
	  b_lowerright=cornerpart+x*zeta*dtildelm;
	  aaa=OOS-zeta*dtildelm;
	  bbb=TOT;
	  ccc=OOS+zeta*dtildelm;

	  tmp_two.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_two,g.ngps_x());
	  // -------- end of applying  Yoverlplustildetwo ------------------

	  // ------- apply Btilde^+ and assemble the stuff -----------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m+1);
	      rhsone[index]=oosqrttwo*(-expp*tmp_one[xindex]+expp*tmp_two[xindex]);
	      rhsone[index_i]=oosqrttwo*(tmp_one[xindex]+tmp_two[xindex]);
	    };
	  // ------- end of applying Btilde^+ etc.  --------------


	  // ------------------ Rtilde -------------------------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      xi=prefaci/g.r(xindex);
	      factor=1.0/(1.0+xi*xi*btildelm*btildelm);
	      diagonal=(1.0-xi*xi*btildelm*btildelm)*factor;
	      ur=2.0*xi*expp*btildelm*factor;
	      ll=-2.0*xi*expm*btildelm*factor;
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m+1);
	      start[index]=diagonal*rhsone[index]+ur*rhsone[index_i];
	      start[index_i]=diagonal*rhsone[index_i]+ll*rhsone[index];
	    };
	  // ------- end of --- Rtilde  -------------------------


	  // ------- apply B and construct Gamma_1 and Gamma_2 -----------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m-1);
	      tmp_one[xindex]=oosqrttwo*(-expp*start[index]+start[index_i]);
	      tmp_two[xindex]=oosqrttwo*(expp*start[index]+start[index_i]);
	    };
	  // ------- end of applying B etc.  --------------


	  dlm=sqrt((ell-m+1.0)*(ell-m+2.0)/((2.0*ell+1.0)*(2.0*ell+3.0)));

	  // ------- calculate Yoverlminusone times Gamma_1
	  xindex=0;
	  tmptmp_one[xindex]=(cornerpart+x*zeta*dlm)*tmp_one[xindex]
	    +(OOS+zeta*dlm)*tmp_one[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_one[xindex]=(OOS-zeta*dlm)*tmp_one[xindex-1]
		+TOT*tmp_one[xindex]
		+(OOS+zeta*dlm)*tmp_one[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_one[xindex]=(OOS-zeta*dlm)*tmp_one[xindex-1]
	    +(cornerpart-x*zeta*dlm)*tmp_one[xindex];
	  // -------- end of calculating Yoverlminusone times Gamma_1


	  // ------- calculate Yoverlminustwo times Gamma_2
	  xindex=0;
	  tmptmp_two[xindex]=(cornerpart-x*zeta*dlm)*tmp_two[xindex]
	    +(OOS-zeta*dlm)*tmp_two[xindex+1];
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmptmp_two[xindex]=(OOS+zeta*dlm)*tmp_two[xindex-1]
		+TOT*tmp_two[xindex]
		+(OOS-zeta*dlm)*tmp_two[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  tmptmp_two[xindex]=(OOS+zeta*dlm)*tmp_two[xindex-1]
	    +(cornerpart+x*zeta*dlm)*tmp_two[xindex];
	  // -------- end of calculating Yoverlminustwo times Gamma_2


	  // -------- apply Yoverlplusone -----------------------------
	  b_upperleft=cornerpart-x*zeta*dlm;
	  b_lowerright=cornerpart+x*zeta*dlm;
	  aaa=OOS-zeta*dlm;
	  bbb=TOT;
	  ccc=OOS+zeta*dlm;

	  tmp_one.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_one,g.ngps_x());
	  // -------- end of applying  Yoverlplusone ------------------


	  // -------- apply Yoverlplustwo -----------------------------
	  b_upperleft=cornerpart+x*zeta*dlm;
	  b_lowerright=cornerpart-x*zeta*dlm;
	  aaa=OOS+zeta*dlm;
	  bbb=TOT;
	  ccc=OOS-zeta*dlm;

	  tmp_two.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmptmp_two,g.ngps_x());
	  // -------- end of applying  Yoverlplustwo ------------------

	  // ------- apply B^+ and assemble the stuff -----------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m-1);
	      rhsone[index]=oosqrttwo*(-expm*tmp_one[xindex]+expm*tmp_two[xindex]);
	      rhsone[index_i]=oosqrttwo*(tmp_one[xindex]+tmp_two[xindex]);
	    };
	  // ------- end of applying B^+ etc.  --------------

	  blm=sqrt((ell-m+1.0)/((2.0*ell+1)*(2.0*ell+3)))
	    *(m*sqrt(ell-m+2)-sqrt((ell+m+1)*((ell+1)*(ell+2)-m*(m-1))));

	  
	  // ------------------ R -------------------------
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      xi=prefaci/g.r(xindex);
	      factor=1.0/(1.0+xi*xi*blm*blm);
	      diagonal=(1.0-xi*xi*blm*blm)*factor;
	      ur=2.0*xi*expm*blm*factor;
	      ll=-2.0*xi*expp*blm*factor;
	      index=g.rlmindex(xindex,ell,m);
	      index_i=g.rlmindex(xindex,ell+1,m-1);
	      start[index]=diagonal*rhsone[index]+ur*rhsone[index_i];
	      start[index_i]=diagonal*rhsone[index_i]+ll*rhsone[index];
	    };
	  // ------- end of --- R -------------------------

	};
    };


  //====== end of ============ the left hand side ======================




}


void wavefunction::do_muller_general_tddft(complex timestep, 
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
					   int propornot[])
{

  double Eff_Charge = charge;

  long xindex, yindex, zindex;
  long index, index_xp, index_xm, index_lp;
  long max_noofgridpoints=g.ngps_x();
  double r;
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double wfsq;
  wavefunction aa(g.ngps_x());
  wavefunction bb(g.ngps_x());
  wavefunction cc(g.ngps_x());	
  wavefunction tmp_plus(g.ngps_x());	
  wavefunction tmp_minus(g.ngps_x());	
  wavefunction tmp_plus_two(g.ngps_x());	
  wavefunction tmp_minus_two(g.ngps_x());	
  complex imagi(0.0,1.0);
  complex a,wnlm,vnlm,Delta_two_upperleft, M_two_upperleft;
  double c, cl, clm, plm, qlm, tl, ctilde, wnlmwnlm, vnlmvnlm, m, mm;
  double aaa,bbb,ccc,Mtilde;
  double b_upperleft, b_lowerright;
  double factor;
  double vecpotwithprefactor;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  double maxpsi;
  complex ul,ur,ll,lr;
  complex aaaa, bbbb, cccc;

  complex halfimagitimestep=(complex)0.5*timestep*imagi;
  complex halfimagitimestepOOS=halfimagitimestep*(complex)OOS;
  complex halfimagitimestepFOT=halfimagitimestep*(complex)FOT;
  
  complex energ;


  
  // to each Kohn-Sham orbital apply the same propagation (since the Hamiltonian is the same)
  for (zindex=0; zindex<g.ngps_z(); zindex++) 
    {
      if ((propornot[zindex]==1) || (time<=0.0))
	{
	  m=ms[zindex];
	  mm=m*m; 
	  //	  std::cout << "Propagating orbital no. " << zindex << " with m = " << m << std::endl;
	  
	  // calculate one l,l+2-block
	  for (xindex=0; xindex<max_noofgridpoints; xindex++)
	    {      
	      r=g.r(xindex);  
	      for (yindex=(long)(fabs(m)); yindex<g.ngps_y()-2; yindex++)
		{
		  qlm=3.0/(2.0*(2.0*yindex+3.0))*sqrt((((yindex+1.0)*(yindex+1.0)-m*m)*((yindex+2.0)*(yindex+2.0)-m*m))/((2.0*yindex+1)*(2.0*yindex+5)));
		  index=g.index(xindex,yindex,zindex);
		  index_lp=g.index(xindex,yindex+2,zindex);
		  vnlm=imagi*0.25*timestep*qlm*wf_three[xindex];
		  vnlmvnlm=abs(vnlm)*abs(vnlm);
		  factor=1.0/(1.0+vnlmvnlm);
		  ul=1.0-vnlmvnlm;
		  ur=-(complex)2.0*vnlm;
		  ll=(complex)2.0*conj(vnlm);
		  lr=1.0-vnlmvnlm;
		  a=factor*(ul*start[index]+ur*start[index_lp]);
		  start[index_lp]=factor*(ll*start[index]+lr*start[index_lp]);
		  start[index]=a;
		};
	    };


	  // calculate one l,l+1-block
	  for (yindex=(long)(fabs(m)); yindex<g.ngps_y()-1; yindex++)
	    {
	      clm=sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	      c=real((complex)0.25*timestep*(complex)hamil.vecpot_z(time,me)*(complex)(clm/(2.0*g.delt_x())));
	      // calculate G_nl Psi
	      tl=(yindex+1.0)*sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	      for (xindex=0; xindex<max_noofgridpoints; xindex++)
		{      
		  r=g.r(xindex);  
		  index=g.index(xindex,yindex,zindex);
		  index_lp=g.index(xindex,yindex+1,zindex);
		  wnlm=(complex)0.25*timestep*(complex)hamil.vecpot_z(time,me)*(complex)(tl/g.r(xindex))+imagi*(complex)0.25*timestep*(complex)clm*(wf_two[xindex]+r*hamil.field(time,me));
		  wnlmwnlm=abs(wnlm)*abs(wnlm);
		  factor=1.0/(2.0*(1.0+wnlmwnlm));
		  ul=1.0+(complex)2.0*conj(wnlm)-wnlmwnlm;
		  ur=1.0-(complex)2.0*wnlm-wnlmwnlm;
		  ll=-1.0+(complex)2.0*conj(wnlm)+wnlmwnlm;
		  lr=1.0+(complex)2.0*wnlm-wnlmwnlm;
		  tmp_plus[xindex]=((complex)ul*start[index]+(complex)ur*start[index_lp])*(complex)factor;
		  tmp_minus[xindex]=((complex)ll*start[index]+(complex)lr*start[index_lp])*(complex)factor;
		};
	      // calculate Y_l^- tmp
	      xindex=0;
	      Mtilde=TOT+OOS*lambda;
	      ctilde=lambda*c;
	      tmp_plus_two[xindex]=(complex)(Mtilde-ctilde)*tmp_plus[xindex] + (complex)(OOS-c)*tmp_plus[xindex+1];
	      tmp_minus_two[xindex]=(complex)(Mtilde+ctilde)*tmp_minus[xindex] + (complex)(OOS+c)*tmp_minus[xindex+1];
	      ul=OOS+c;
	      ur=OOS-c;
	      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
		{
		  tmp_plus_two[xindex]=(complex)ul*tmp_plus[xindex-1] + (complex)TOT*tmp_plus[xindex] + (complex)ur*tmp_plus[xindex+1];
		  tmp_minus_two[xindex]=(complex)ur*tmp_minus[xindex-1] + (complex)TOT*tmp_minus[xindex] + ul*tmp_minus[xindex+1];
		};
	      xindex=max_noofgridpoints-1;
	      Mtilde=TOT+OOS*lambda;
	      ctilde=llambda*c;
	      tmp_plus_two[xindex]=(complex)(OOS+c)*tmp_plus[xindex-1] + (complex)(Mtilde-ctilde)*tmp_plus[xindex];
	      tmp_minus_two[xindex]=(complex)(OOS-c)*tmp_minus[xindex-1] + (complex)(Mtilde+ctilde)*tmp_minus[xindex];
	      // solve for Y_l^+ tmp_new = tmp_old
	      Mtilde=TOT+OOS*lambda;
	      ctilde=lambda*c;
	      b_upperleft=(Mtilde+ctilde);
	      Mtilde=TOT+OOS*lambda;
	      ctilde=llambda*c;
	      b_lowerright=(Mtilde+ctilde);
	      aaa=OOS+c;
	      bbb=TOT;
	      ccc=OOS-c;
	      tmp_plus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_plus_two,max_noofgridpoints);
	      Mtilde=TOT+OOS*lambda;
	      ctilde=lambda*c;
	      b_upperleft=(Mtilde-ctilde);
	      Mtilde=TOT+OOS*lambda;
	      ctilde=llambda*c;
	      b_lowerright=(Mtilde-ctilde);
	      aaa=OOS-c;
	      bbb=TOT;
	      ccc=OOS+c;
	      tmp_minus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_minus_two,max_noofgridpoints);
	      // update the wavefunction
	      // apply B^T (prefactor 1/sqrt(2) already accounted for above)
	      for (xindex=0; xindex<max_noofgridpoints; xindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  index_lp=g.index(xindex,yindex+1,zindex);
		  start[index]=(tmp_plus[xindex]-tmp_minus[xindex]);
		  start[index_lp]=(tmp_plus[xindex]+tmp_minus[xindex]);
		};
	    };
	  
	  // the spatial part
	  // The constant part of the matrix L_-(tau)
	  aaaa=-OOS-halfimagitimestep/(g.delt_x()*g.delt_x());
	  cccc=aaaa;
	  bbbb=-FOT+imagi*timestep/(g.delt_x()*g.delt_x());
	  Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
	  M_two_upperleft=-(complex)2.0*(1.0+(complex)g.delt_x()*(complex)g.delt_x()/12.0*Delta_two_upperleft);
	  
	  // Calculate the rhs vector W_-(tau) *this
	  for (yindex=(long)(fabs(m)); yindex<g.ngps_y(); yindex++)
	    { 
	      xindex=0;
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      plm=(yindex*(yindex+1.0)-3*mm)/((2.0*yindex-1.0)*(2.0*yindex+3.0));
	      if (yindex==0)
		{
		  rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]+wf_one[xindex+1]+plm*wf_three[xindex+1]))*start[index_xp]
		    +((complex)M_two_upperleft*(1.0-halfimagitimestep*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]))-imagi*(complex)Delta_two_upperleft*(complex)0.5*timestep)*start[index];
		}
	      else
		{
		  rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]+wf_one[xindex+1]+plm*wf_three[xindex+1]))*start[index_xp]
		    +(bbbb+halfimagitimestepFOT*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]))*start[index];
		};
	      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  index_xp=g.index(xindex+1,yindex,zindex);
		  index_xm=g.index(xindex-1,yindex,zindex);
		  rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]+wf_one[xindex+1]+plm*wf_three[xindex+1]))*start[index_xp]
		    +(bbbb+halfimagitimestepFOT*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]))*start[index]
		    +(cccc+halfimagitimestepOOS*(staticpot[index_xm]+wf_one[xindex-1]+plm*wf_three[xindex-1]))*start[index_xm];
		};
	      
	      
	      xindex=max_noofgridpoints-1;
	      index=g.index(xindex,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      rhsone[index]=(bbbb+halfimagitimestepFOT*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]))*start[index]
		+(cccc+halfimagitimestepOOS*(staticpot[index_xm]+wf_one[xindex-1]+plm*wf_three[xindex-1]))*start[index_xm];
	    };
	  
	  // The matrix W_+(tau)
	  aaaa=-OOS+halfimagitimestep/(g.delt_x()*g.delt_x());
	  cccc=aaaa;
	  bbbb=-FOT-imagi*timestep/(g.delt_x()*g.delt_x());
	  
	  for (yindex=(long)(fabs(m)); yindex<g.ngps_y(); yindex++)
	    { 
	      xindex=0;
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      plm=(yindex*(yindex+1.0)-3*mm)/((2.0*yindex-1.0)*(2.0*yindex+3.0));
	      aa[xindex]=(aaaa-halfimagitimestep*(complex)OOS*(staticpot[index_xp]+wf_one[xindex+1]+plm*wf_three[xindex+1]));
	      cc[xindex]=1.0; // not used
	      
	      if (yindex==0)
		{
		  bb[xindex]=((complex)M_two_upperleft*(1.0+halfimagitimestep*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]))+imagi*(complex)Delta_two_upperleft*(complex)0.5*timestep);
		}
	      else
		{
		  bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]));
		};
	      
	      tmp_plus[xindex]=rhsone[index];
	      
	      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  index_xp=g.index(xindex+1,yindex,zindex);
		  index_xm=g.index(xindex-1,yindex,zindex);
		  
		  aa[xindex]=(aaaa-halfimagitimestepOOS*(staticpot[index_xp]+wf_one[xindex+1]+plm*wf_three[xindex+1]));
		  cc[xindex]=(cccc-halfimagitimestepOOS*(staticpot[index_xm]+wf_one[xindex-1]+plm*wf_three[xindex-1]));
		  bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]));
		  tmp_plus[xindex]=rhsone[index];
		};  
	      
	      xindex=max_noofgridpoints-1;
	      r=g.r(xindex);
	      index=g.index(xindex,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      aa[xindex]=aaaa;
	      cc[xindex]=(cccc-halfimagitimestepOOS*(staticpot[index_xm]+wf_one[xindex-1]+plm*wf_three[xindex-1]));
	      bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]+wf_one[xindex]+plm*wf_three[xindex]));
	      tmp_plus[xindex]=rhsone[index];
	      
	      tmp_minus.solve(aa,bb,cc,tmp_plus,max_noofgridpoints);
	      
	      for (xindex=0; xindex<max_noofgridpoints; xindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  start[index]=tmp_minus[xindex];
		};  
	    };
	  
	  
	  // calculate one l,l+1-block
	  for (yindex=g.ngps_y()-2; yindex>=(long)(fabs(m)); yindex--)
	    {
	      clm=sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	      c=real((complex)0.25*timestep*(complex)hamil.vecpot_z(real(time+timestep),me)*(complex)clm/(2.0*g.delt_x()));
	      // calculate B_l Psi (without prefactor 1/sqrt(2) (see below))
	      for (xindex=0; xindex<max_noofgridpoints; xindex++)
		{      
		  index=g.index(xindex,yindex,zindex);
		  index_lp=g.index(xindex,yindex+1,zindex);
		  tmp_plus[xindex]=(start[index]+start[index_lp]);
		  tmp_minus[xindex]=(-start[index]+start[index_lp]);
		};
	      // calculate Y_l^- tmp
	      xindex=0;
	      Mtilde=TOT+OOS*lambda;
	      ctilde=lambda*c;
	      tmp_plus_two[xindex]=(complex)(Mtilde-ctilde)*tmp_plus[xindex] + (complex)(OOS-c)*tmp_plus[xindex+1];
	      tmp_minus_two[xindex]=(complex)(Mtilde+ctilde)*tmp_minus[xindex] + (complex)(OOS+c)*tmp_minus[xindex+1];
	      ul=OOS+c;
	      ur=OOS-c;
	      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
		{
		  tmp_plus_two[xindex]=(complex)ul*tmp_plus[xindex-1] + (complex)TOT*tmp_plus[xindex] + ur*tmp_plus[xindex+1];
		  tmp_minus_two[xindex]=(complex)ur*tmp_minus[xindex-1] + (complex)TOT*tmp_minus[xindex] + ul*tmp_minus[xindex+1];
		};
	      xindex=max_noofgridpoints-1;
	      Mtilde=TOT+OOS*lambda;
	      ctilde=llambda*c;
	      tmp_plus_two[xindex]=(complex)(OOS+c)*tmp_plus[xindex-1] + (complex)(Mtilde-ctilde)*tmp_plus[xindex];
	      tmp_minus_two[xindex]=(complex)(OOS-c)*tmp_minus[xindex-1] + (complex)(Mtilde+ctilde)*tmp_minus[xindex];
	      // solve for Y_l^+ tmp_new = tmp_old
	      Mtilde=TOT+OOS*lambda;
	      ctilde=lambda*c;
	      b_upperleft=(Mtilde+ctilde);
	      Mtilde=TOT+OOS*lambda;
	      ctilde=llambda*c;
	      b_lowerright=(Mtilde+ctilde);
	      aaa=OOS+c;
	      bbb=TOT;
	      ccc=OOS-c;
	      tmp_plus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_plus_two,max_noofgridpoints);
	      Mtilde=TOT+OOS*lambda;
	      ctilde=lambda*c;
	      b_upperleft=(Mtilde-ctilde);
	      Mtilde=TOT+OOS*lambda;
	      ctilde=llambda*c;
	      b_lowerright=(Mtilde-ctilde);
	      aaa=OOS-c;
	      bbb=TOT;
	      ccc=OOS+c;
	      tmp_minus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_minus_two,max_noofgridpoints);
	      // calculate S_nl tmp
	      tl=(yindex+1.0)*sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	      for (xindex=0; xindex<max_noofgridpoints; xindex++)
		{
		  r=g.r(xindex);
		  index=g.index(xindex,yindex,zindex);
		  index_lp=g.index(xindex,yindex+1,zindex);
		  wnlm=(complex)0.25*timestep*(complex)hamil.vecpot_z(real(time+timestep),me)*(complex)tl/g.r(xindex)+imagi*(complex)0.25*timestep*(complex)clm*(wf_two[xindex]+r*hamil.field(real(time+timestep),me));
		  wnlmwnlm=abs(wnlm)*abs(wnlm);
		  factor=1.0/(2.0*(1.0+wnlmwnlm));
		  ll=1.0+(complex)2.0*conj(wnlm)-wnlmwnlm;
		  ul=1.0-(complex)2.0*wnlm-wnlmwnlm;
		  ur=-1.0-(complex)2.0*wnlm+wnlmwnlm;
		  lr=1.0-(complex)2.0*conj(wnlm)-wnlmwnlm;
		  start[index]=((complex)ul*tmp_plus[xindex]+(complex)ur*tmp_minus[xindex])*(complex)factor;
		  start[index_lp]=((complex)ll*tmp_plus[xindex]+(complex)lr*tmp_minus[xindex])*(complex)factor;
		};
	    };
	  // calculate one l,l+2-block
	  for (xindex=0; xindex<max_noofgridpoints; xindex++)
	    {      
	      r=g.r(xindex);  
	      for (yindex=g.ngps_y()-3; yindex>=(long)(fabs(m)); yindex--)
		{
		  qlm=3.0/(2.0*(2.0*yindex+3.0))*sqrt((((yindex+1.0)*(yindex+1.0)-m*m)*((yindex+2.0)*(yindex+2.0)-m*m))/((2.0*yindex+1)*(2.0*yindex+5)));
		  index=g.index(xindex,yindex,zindex);
		  index_lp=g.index(xindex,yindex+2,zindex);
		  vnlm=imagi*0.25*timestep*qlm*wf_three[xindex];
		  vnlmvnlm=abs(vnlm)*abs(vnlm);
		  factor=1.0/(1.0+vnlmvnlm);
		  ul=1.0-vnlmvnlm;
		  ur=-(complex)2.0*vnlm;
		  ll=(complex)2.0*conj(vnlm);
		  lr=1.0-vnlmvnlm;
		  a=factor*(ul*start[index]+ur*start[index_lp]);
		  start[index_lp]=factor*(ll*start[index]+lr*start[index_lp]);
		  start[index]=a;
		};
	    };
	};
    };
}



void wavefunction::do_muller_ell(complex timestep, 
				 double time, 
				 grid g, 
				 hamop hamil, 
				 const wavefunction &staticpot, 
				 int me, 
				 double charge,
				 int m)
{

  double Eff_Charge = charge;

  long xindex, yindex, zindex;
  long index, index_xp, index_xm, index_lp;
  double r;
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double wfsq;
  wavefunction aa(g.ngps_x());
  wavefunction bb(g.ngps_x());
  wavefunction cc(g.ngps_x());	
  wavefunction tmp_plus(g.ngps_x());	
  wavefunction tmp_minus(g.ngps_x());	
  wavefunction tmp_plus_two(g.ngps_x());	
  wavefunction tmp_minus_two(g.ngps_x());	
  complex imagi(0.0,1.0);
  complex wnl,Delta_two_upperleft, M_two_upperleft;
  double c, cl, cnl, tl, ctilde, wnlwnl, mm;
  double aaa,bbb,ccc,Mtilde;
  double b_upperleft, b_lowerright;
  double factor;
  double vecpotwithprefactor;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  double maxpsi;
  complex ul,ur,ll,lr;
  complex aaaa, bbbb, cccc;

  complex halfimagitimestep=(complex)0.5*timestep*imagi;
  complex halfimagitimestepOOS=halfimagitimestep*(complex)OOS;
  complex halfimagitimestepFOT=halfimagitimestep*(complex)FOT;
  
  complex energ;

  mm=m*m; 
  
  for (zindex=0; zindex<g.ngps_z(); zindex++) 
    {
      // calculate one l-block
      for (yindex=m; yindex<g.ngps_y()-1; yindex++)
	{
	  cnl=sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	  c=real((complex)0.25*timestep*(complex)hamil.vecpot_z(time,me)*(complex)(cnl/(2.0*g.delt_x())));
	  // calculate G_nl Psi
	  tl=(yindex+1.0)*sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {      
	      index=g.index(xindex,yindex,zindex);
	      index_lp=g.index(xindex,yindex+1,zindex);
	      wnl=(complex)0.25*timestep*(complex)hamil.vecpot_z(time,me)*(complex)(tl/g.r(xindex));
	      wnlwnl=abs(wnl)*abs(wnl);
	      factor=1.0/(2.0*(1.0+wnlwnl));
	      ul=1.0+(complex)2.0*conj(wnl)-wnlwnl;
	      ur=1.0-(complex)2.0*wnl-wnlwnl;
	      ll=-1.0+(complex)2.0*conj(wnl)+wnlwnl;
	      lr=1.0+(complex)2.0*wnl-wnlwnl;
	      tmp_plus[xindex]=((complex)ul*start[index]+(complex)ur*start[index_lp])*(complex)factor;
	      tmp_minus[xindex]=((complex)ll*start[index]+(complex)lr*start[index_lp])*(complex)factor;
	    };
	  // calculate Y_l^- tmp
	  xindex=0;
	  Mtilde=TOT+OOS*lambda;
	  ctilde=lambda*c;
	  tmp_plus_two[xindex]=(complex)(Mtilde-ctilde)*tmp_plus[xindex] + (complex)(OOS-c)*tmp_plus[xindex+1];
	  tmp_minus_two[xindex]=(complex)(Mtilde+ctilde)*tmp_minus[xindex] + (complex)(OOS+c)*tmp_minus[xindex+1];
	  ul=OOS+c;
	  ur=OOS-c;
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmp_plus_two[xindex]=(complex)ul*tmp_plus[xindex-1] + (complex)TOT*tmp_plus[xindex] + (complex)ur*tmp_plus[xindex+1];
	      tmp_minus_two[xindex]=(complex)ur*tmp_minus[xindex-1] + (complex)TOT*tmp_minus[xindex] + ul*tmp_minus[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  Mtilde=TOT+OOS*lambda;
	  ctilde=llambda*c;
	  tmp_plus_two[xindex]=(complex)(OOS+c)*tmp_plus[xindex-1] + (complex)(Mtilde-ctilde)*tmp_plus[xindex];
	  tmp_minus_two[xindex]=(complex)(OOS-c)*tmp_minus[xindex-1] + (complex)(Mtilde+ctilde)*tmp_minus[xindex];
	  // solve for Y_l^+ tmp_new = tmp_old
	  Mtilde=TOT+OOS*lambda;
	  ctilde=lambda*c;
	  b_upperleft=(Mtilde+ctilde);
	  Mtilde=TOT+OOS*lambda;
	  ctilde=llambda*c;
	  b_lowerright=(Mtilde+ctilde);
	  aaa=OOS+c;
	  bbb=TOT;
	  ccc=OOS-c;
	  tmp_plus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_plus_two,g.ngps_x());
	  Mtilde=TOT+OOS*lambda;
	  ctilde=lambda*c;
	  b_upperleft=(Mtilde-ctilde);
	  Mtilde=TOT+OOS*lambda;
	  ctilde=llambda*c;
	  b_lowerright=(Mtilde-ctilde);
	  aaa=OOS-c;
	  bbb=TOT;
	  ccc=OOS+c;
	  tmp_minus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_minus_two,g.ngps_x());
	  // update the wavefunction
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      index_lp=g.index(xindex,yindex+1,zindex);
	      start[index]=(tmp_plus[xindex]-tmp_minus[xindex]);
	      start[index_lp]=(tmp_plus[xindex]+tmp_minus[xindex]);
	    };
	};
      
      // the spatial part
      // The constant part of the matrix L_-(tau)
      aaaa=-OOS-halfimagitimestep/(g.delt_x()*g.delt_x());
      cccc=aaaa;
      bbbb=-FOT+imagi*timestep/(g.delt_x()*g.delt_x());
      Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
      M_two_upperleft=-(complex)2.0*(1.0+(complex)g.delt_x()*(complex)g.delt_x()/12.0*Delta_two_upperleft);
      
      // Calculate the rhs vector W_-(tau) *this
      for (yindex=m; yindex<g.ngps_y(); yindex++)
	{ 
	  xindex=0;
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  if (yindex==0)
	    {
	      rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]))*start[index_xp]
		+((complex)M_two_upperleft*(1.0-halfimagitimestep*(staticpot[index]))-imagi*(complex)Delta_two_upperleft*(complex)0.5*timestep)*start[index];
	    }
	  else
	    {
	      rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]))*start[index_xp]
		+(bbbb+halfimagitimestepFOT*(staticpot[index]))*start[index];
	    };
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    { 
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      rhsone[index]=(aaaa+halfimagitimestepOOS*(staticpot[index_xp]))*start[index_xp]
		+(bbbb+halfimagitimestepFOT*(staticpot[index]))*start[index]
		+(cccc+halfimagitimestepOOS*(staticpot[index_xm]))*start[index_xm];
	    };
	  
	  
	  xindex=g.ngps_x()-1;
	  index=g.index(xindex,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  rhsone[index]=(bbbb+halfimagitimestepFOT*(staticpot[index]))*start[index]
	    +(cccc+halfimagitimestepOOS*(staticpot[index_xm]))*start[index_xm];
	};
      
      // The matrix W_+(tau)
      aaaa=-OOS+halfimagitimestep/(g.delt_x()*g.delt_x());
      cccc=aaaa;
      bbbb=-FOT-imagi*timestep/(g.delt_x()*g.delt_x());
      
      for (yindex=m; yindex<g.ngps_y(); yindex++)
	{ 
	  xindex=0;
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  aa[xindex]=(aaaa-halfimagitimestep*(complex)OOS*(staticpot[index_xp]));
	  cc[xindex]=1.0; // not used
	  
	  if (yindex==0)
	    {
	      bb[xindex]=((complex)M_two_upperleft*(1.0+halfimagitimestep*(staticpot[index]))+imagi*(complex)Delta_two_upperleft*(complex)0.5*timestep);
	    }
	  else
	    {
	      bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]));
	    };
	  
	  tmp_plus[xindex]=rhsone[index];
	  
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    { 
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      
	      aa[xindex]=(aaaa-halfimagitimestepOOS*(staticpot[index_xp]));
	      cc[xindex]=(cccc-halfimagitimestepOOS*(staticpot[index_xm]));
	      bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]));
	      tmp_plus[xindex]=rhsone[index];
	    };  
	  
	  xindex=g.ngps_x()-1;
	  r=g.r(xindex);
	  index=g.index(xindex,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  aa[xindex]=aaaa-halfimagitimestepOOS*(hamil.scalarpotx(r+g.delt_x(),yindex,0.0,time,me)
						-imagi*(complex)hamil.imagpot(xindex,yindex,0,time,g));
	  cc[xindex]=(cccc-halfimagitimestepOOS*(staticpot[index_xm]));
	  bb[xindex]=(bbbb-halfimagitimestepFOT*(staticpot[index]));
	  tmp_plus[xindex]=rhsone[index];
	  
	  tmp_minus.solve(aa,bb,cc,tmp_plus,g.ngps_x());
	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      index=g.index(xindex,yindex,zindex);
	      start[index]=tmp_minus[xindex];
	    };  
	};
      
      
      // calculate one l-block
      for (yindex=g.ngps_y()-2; yindex>=m; yindex--)
	{
	  cnl=sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	  c=real((complex)0.25*timestep*(complex)hamil.vecpot_z(real(time+timestep),me)*(complex)cnl/(2.0*g.delt_x()));
	  // calculate B_l Psi
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {      
	      index=g.index(xindex,yindex,zindex);
	      index_lp=g.index(xindex,yindex+1,zindex);
	      tmp_plus[xindex]=(start[index]+start[index_lp]);
	      tmp_minus[xindex]=(-start[index]+start[index_lp]);
	    };
	  // calculate Y_l^- tmp
	  xindex=0;
	  Mtilde=TOT+OOS*lambda;
	  ctilde=lambda*c;
	  tmp_plus_two[xindex]=(complex)(Mtilde-ctilde)*tmp_plus[xindex] + (complex)(OOS-c)*tmp_plus[xindex+1];
	  tmp_minus_two[xindex]=(complex)(Mtilde+ctilde)*tmp_minus[xindex] + (complex)(OOS+c)*tmp_minus[xindex+1];
	  ul=OOS+c;
	  ur=OOS-c;
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      tmp_plus_two[xindex]=(complex)ul*tmp_plus[xindex-1] + (complex)TOT*tmp_plus[xindex] + ur*tmp_plus[xindex+1];
	      tmp_minus_two[xindex]=(complex)ur*tmp_minus[xindex-1] + (complex)TOT*tmp_minus[xindex] + ul*tmp_minus[xindex+1];
	    };
	  xindex=g.ngps_x()-1;
	  Mtilde=TOT+OOS*lambda;
	  ctilde=llambda*c;
	  tmp_plus_two[xindex]=(complex)(OOS+c)*tmp_plus[xindex-1] + (complex)(Mtilde-ctilde)*tmp_plus[xindex];
	  tmp_minus_two[xindex]=(complex)(OOS-c)*tmp_minus[xindex-1] + (complex)(Mtilde+ctilde)*tmp_minus[xindex];
	  // solve for Y_l^+ tmp_new = tmp_old
	  Mtilde=TOT+OOS*lambda;
	  ctilde=lambda*c;
	  b_upperleft=(Mtilde+ctilde);
	  Mtilde=TOT+OOS*lambda;
	  ctilde=llambda*c;
	  b_lowerright=(Mtilde+ctilde);
	  aaa=OOS+c;
	  bbb=TOT;
	  ccc=OOS-c;
	  tmp_plus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_plus_two,g.ngps_x());
	  Mtilde=TOT+OOS*lambda;
	  ctilde=lambda*c;
	  b_upperleft=(Mtilde-ctilde);
	  Mtilde=TOT+OOS*lambda;
	  ctilde=llambda*c;
	  b_lowerright=(Mtilde-ctilde);
	  aaa=OOS-c;
	  bbb=TOT;
	  ccc=OOS+c;
	  tmp_minus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_minus_two,g.ngps_x());
	  // calculate S_nl tmp
	  tl=(yindex+1.0)*sqrt((yindex+1.0)*(yindex+1.0)-mm)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      index_lp=g.index(xindex,yindex+1,zindex);
	      wnl=(complex)0.25*timestep*(complex)hamil.vecpot_z(real(time+timestep),me)*(complex)tl/g.r(xindex);
	      wnlwnl=abs(wnl)*abs(wnl);
	      factor=1.0/(2.0*(1.0+wnlwnl));
	      ll=1.0+(complex)2.0*conj(wnl)-wnlwnl;
	      ul=1.0-(complex)2.0*wnl-wnlwnl;
	      ur=-1.0-(complex)2.0*wnl+wnlwnl;
	      lr=1.0-(complex)2.0*conj(wnl)-wnlwnl;
	      start[index]=((complex)ul*tmp_plus[xindex]+(complex)ur*tmp_minus[xindex])*(complex)factor;
	      start[index_lp]=((complex)ll*tmp_plus[xindex]+(complex)lr*tmp_minus[xindex])*(complex)factor;
	    };
	};
    };
}







void wavefunction::dump_xvector_to_file(long ngpsx,  FILE* os, int stepwidth)
{
  for (long i=0; i<ngpsx; i=i+stepwidth)
    {
      fprintf(os,"%16.12e %16.12e\n",real(start[i]),imag(start[i]));
    };

}


void wavefunction::dump_to_file(grid g, FILE* os, int stepwidth, fluid &degeneracies, fluid &ms, double fact)
{

  long xindex,yindex,zindex,i,counter, counter_ii, counter_iii,yyindex,ii;
  double u, r, rho, z, legpolnew, legpolold, legpololder;
  complex summingres;
  
  
  counter=0;
  counter_ii=0;
  counter_iii=0;
  
  switch (g.dimens()) {
  case 34 :
    for (zindex=0; zindex<g.ngps_z(); zindex++) {
      counter=0;
      for (rho=g.delt_x(); rho<fact*g.r(g.ngps_x()); rho+=g.delt_x()*stepwidth) {
	//	      std::cout << rho << std::endl;
	counter++;
	counter_ii=0;
	for (z=-g.r(g.ngps_x())*fact; z<g.r(g.ngps_x())*fact; z+=g.delt_x()*stepwidth) {
	  //		  std::cout << "-- " << z << std::endl;
	  counter_ii++;
	  r = sqrt(rho*rho+z*z);
	  u = z/r;
	  xindex = g.rindex(r);
	  if( xindex < g.ngps_x() ) {
	    summingres = complex( 0.0, 0.0 );
	    for( yindex=(long)ms[zindex]; yindex<g.ngps_y(); yindex++ )  {
	      for( yyindex=(long)ms[zindex]; yyindex<g.ngps_y(); yyindex++ ) {
		i = g.index( xindex, yindex, zindex );
		ii = g.index( xindex, yyindex, zindex );
		if( ms[zindex] == 0.0 ) {
		  summingres = summingres + start[i]*conj(start[ii])
		    * ylm( yindex, 0, acos(u), 0.0) 
		    * ylm( yyindex, 0, acos(u), 0.0);
		  // Here, ylm is always real since phi=0.0 !
		}
		else {
		  summingres=summingres+start[i]*conj(start[ii])*2.0
		    *ylm(yindex,(int)ms[zindex],acos(u),0.0)
		    *ylm(yyindex,(int)ms[zindex],acos(u),0.0);
		}
	      }
	    }
	    fprintf(os, "%e\n", 2.0*real(summingres)/(g.r(xindex)*g.r(xindex)));
	  }
	  else {
	    fprintf(os, "%e\n", 0.0); 
	  }
	}
      }
      std::cout << "transforming orbital no. " << zindex << std::endl;
    };
    break;
  };
  
  std::cout << "dumped " << counter << ", " << counter_ii << ", " << counter_iii << std::endl;
  std::cout << "rho_max=" << fact*g.r(g.ngps_x()) << ", " << "z=[" << -g.r(g.ngps_x())*fact
	    << ", " << g.r(g.ngps_x())*fact << "]" << std::endl;
}


void wavefunction::dump_to_file_sh(grid g, FILE* os, int stepwidth)
{
  dump_to_file_sh(g, os, stepwidth, 1);
}

/*!

 */
int wavefunction::dump_to_file_sh(grid g, FILE* os, int stepwidth, int iv)
{
  long index;

  if (iv==1) {fprintf(stdout, "I'm doing a sh-output... "); fflush(stdout);}

  switch (g.dimens())
    {
    case 34: case 44:
      for (index=0; index<wf_size(); index++)
	fprintf(os,"%20.15le %20.15le\n",real(start[index]), imag(start[index]));
      fflush(os);
      if (iv==1) {fprintf(stdout, "done!\n"); fflush(stdout);}
      break;

    default :
      fprintf(stderr, "err: unknown propagation mode (g.dimens())\n");
      exit(-54);
    }

  return(0);
}


/*! \fn wavefunction::calculate_kli_zero(grid g, const fluid &Lambdavector,
  const fluid &ells, const fluid &ms, long non_inversion,
  const fluid &degeneracies, int slateronly, int iv)

*/

fluid wavefunction::calculate_kli_zero(grid g, const fluid &Lambdavector,
				       const fluid &ells, const fluid &ms, const fluid &degeneracies,
				       int slateronly, int iv)
{
  double clebsch_i, clebsch_ii, summ, y1, y2, x1, x2, k, b, sign;
  long L, non_inversion;

  int info=0, ikb;

  non_inversion = 1;

  if ((g.ngps_z()<2) && (slateronly!=1))
    {
      fprintf(stderr, "err: kli_zero(): g.ngps_z()==%ld i.e. \
slateronly=1 is only possible.\n", g.ngps_z());
      exit(-52);
    }

  fluid result(g.ngps_x());
  result.nullify();

  complex M_tmp,tmp; 

  grid q;
  q.set_dim(34); // all the same
  q.set_ngps(g.ngps_x(),g.ngps_z(),g.ngps_z());
  q.set_delt(g.delt_x(),1.0,1.0);
  q.set_offs(0,0,0);


  wavefunction Lambda(g.ngps_x()*g.ngps_z()*g.ngps_z());  
  wavefunction rruxsigma(g.ngps_x());
  fluid rhs(g.ngps_z()-non_inversion);

  long j,i,rindex,index;

  Lambda=calculate_Lambda_ki(g,ells);

  rruxsigma=Lambda.calculate_rruxsigma(g,q,ells,degeneracies);
  //  std::cout << "rruxsigma" << std::endl << rruxsigma << std::endl;
  wavefunction uxi(g.ngps_z()-non_inversion);

  if (slateronly==0)
    {
      if (iv==1) std::cout << "I'm solving the matrix equation for KLI..." << std::endl;
      wavefunction vsxj(g.ngps_z()-non_inversion);
      cmatrix M(g.ngps_z()-non_inversion,g.ngps_z()-non_inversion);

      vsxj=Lambda.calculate_vsxj(q,Lambdavector,rruxsigma,non_inversion);
      // std::cout << "vsxj" << endl << vsxj << std::endl;
      uxi=Lambda.calculate_uxi(g,q,non_inversion,ells,degeneracies);
      // std::cout << "uxi" << std::endl << uxi << std::endl;

      // calculate symmetric matrix (delta_ij-M_ij) and the right-hand side
      for (j=0; j<g.ngps_z()-non_inversion; j++)
	{
	  M_tmp=1.0/(0.5*degeneracies[j]) -
	    Lambda.calculate_Mji(q,Lambdavector,ells,ms,j,j);
	  M.set_elem(j,j,M_tmp);
	  for (i=0; i<j; i++)
	    {
	      M_tmp=0.0-Lambda.calculate_Mji(q,Lambdavector,ells,ms,j,i);
	      M.set_elem(j,i,M_tmp);
	      M.set_elem(i,j,M_tmp);
	    }
	  rhs[j]=real(vsxj[j]-uxi[j]);
	}

      if (iv==1){
	std::cout << "rhs before" << rhs << std::endl;
	std::cout << "M" << std::endl << M << std::endl;}

      M.lap_dsytrf_and_dsytrs(rhs,info);

      for (j=0; j<g.ngps_z()-non_inversion; j++)
	rhs[j]=rhs[j]/(0.5*degeneracies[j]);

      if (iv==1){
	std::cout << "rhs after" << std::endl << rhs << std::endl;
	std::cout << "... done!" << std::endl;}

    }
  else
    {
      if (iv==1) std::cout << "I'm using the Slater potential only." << std::endl;
    }
  
  ikb = 0;

  // finally calculate the KLI potential
  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      result[rindex]=real(rruxsigma[rindex]);
      if (slateronly==0)
	{
	  tmp=complex(0.0,0.0);
	  for (j=0; j<g.ngps_z()-non_inversion; j++)
	    {
	      index=q.index(rindex,j,j);
	      tmp+=Lambda[index]*rhs[j]*0.5*degeneracies[j];
	    }
	  result[rindex]+=real(tmp);
	}

      if (Lambdavector[rindex]>1e-150) 
	{
	  result[rindex]*=2.0/Lambdavector[rindex];
	}
      else 
	{  
	  if (ikb==0)
	    {
	      y2 = log(fabs(result[rindex-2]));
	      y1 = log(fabs(result[rindex-3]));
	      x2 =    g.r(rindex-2);
	      x1 =    g.r(rindex-3);

	      k = ( y2 - y1 ) / ( x2 - x1 );
	      b =   y2 - k * x2;

	      ikb = 1;
	      if (result[rindex-2]>0) { sign=1; } else { sign = -1; }
	    }

	  result[rindex]=sign*exp(k*g.r(rindex)+b);
	}
    
    }
  return result;
}





fluid wavefunction::calculate_sic_gam_zero(grid g, const fluid &Lambdavector, 
					   const fluid &hartreezero, 
					   const fluid &ms, 
					   const fluid &degeneracies)
{

  fluid result(g.ngps_x());
  result.nullify();

  double current_norm;
  

  // the grid on which Lambda_i(r,t) and v_i lives
  grid q;
  q.set_dim(34); // all the same
  q.set_ngps(g.ngps_x(),1,g.ngps_z());
  q.set_delt(g.delt_x(),1.0,1.0);
  q.set_offs(0,0,0);

  fluid Lambda(g.ngps_x()*g.ngps_z());
  fluid tmp(g.ngps_x());
  fluid tmp_ii(g.ngps_x());


  Lambda.nullify();
  tmp.nullify();
  tmp_ii.nullify();

  long rindex,iindex,index,lindex,wfindex,jindex,iiindex,jjindex;

  // calculate Lambda_i
  for (iindex=0; iindex<g.ngps_z(); iindex++)
    {
      for (rindex=0; rindex<g.ngps_x(); rindex++)
	{
	  index=q.index(rindex,0,iindex);
	  for (lindex=(long)(fabs(ms[iindex])); lindex<g.ngps_y(); lindex++)
	    {
	      wfindex=g.index(rindex,lindex,iindex);
	      Lambda[index]+=real(conj(start[wfindex])*start[wfindex]);
	    };
	};
    };
  
  current_norm=total_ks_norm(g,degeneracies);

  // calculate the sum of the SIC exchange
  for (iindex=0; iindex<g.ngps_z(); iindex++)
    {
      tmp=Lambda.extract_x(q,0,iindex);
      tmp_ii=tmp_ii+calculate_x_zero(g,2.0*tmp);
    };
  tmp_ii*=1.0/(0.5*current_norm);

  result=-1.0/(0.5*current_norm)*0.5*hartreezero-tmp_ii;

  return result;
}




fluid wavefunction::calculate_hartree_zero(grid g, const fluid &Lambdavector) 
{
  fluid result(g.ngps_x());
  long rindex,rrindex;
  double r;
  double rr;
  double rm;
  double upperpart;
  double lowerpart;



  rindex=0;
  lowerpart=0.0;
  upperpart=0.0;
  for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
    {
      rr=g.r(rrindex);
      upperpart+=Lambdavector[rrindex]/rr;
    };
  result[rindex]=upperpart;

  for (rindex=1; rindex<g.ngps_x(); rindex++)
    {      
      r=g.r(rindex);
      rm=g.r(rindex-1);

      lowerpart=lowerpart*rm/r+Lambdavector[rindex-1]/r;
      upperpart=upperpart-Lambdavector[rindex-1]/rm;

      result[rindex]=lowerpart+upperpart;

    };
      



  return result*g.delt_x();
}


wavefunction wavefunction::calculate_Theta(grid g, const fluid &degeneracies, const fluid &ms)
{
  wavefunction result(g.ngps_x());
  long rindex,lindex,index,indexlpo,indexlmo,zindex;
  double r;
  double rr;
  double m;
  complex tmp(0.0,0.0);


  if (g.ngps_y()<2) fprintf(stdout,"warn: calculate_Theta: g.ngps_y()<2\n");

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      result[rindex]=complex(0.0,0.0);
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  m=ms[zindex];
	  if (g.ngps_y()>(long)(fabs(m)))
	    {

	      lindex=(long)(fabs(m));
	      index=g.index(rindex,lindex,zindex);
	      indexlpo=g.index(rindex,lindex+1,zindex);
	      tmp=(complex)sqrt(((lindex+1.0)*(lindex+1.0)-m*m)/((2.0*lindex+1.0)*(2.0*lindex+3.0)))*conj(start[indexlpo])*start[index];
	      for (lindex=(long)(fabs(m))+1; lindex<g.ngps_y()-1; lindex++)
		{

		  index=g.index(rindex,lindex,zindex);
		  indexlpo=g.index(rindex,lindex+1,zindex);
		  indexlmo=g.index(rindex,lindex-1,zindex);
		  tmp+=(complex)(sqrt(((lindex+1.0)*(lindex+1.0)-m*m)/((2.0*lindex+1.0)*(2.0*lindex+3.0))))
		    *conj(start[indexlpo])*start[index]
		    +(complex)(sqrt((lindex*lindex-m*m)/((2.0*lindex+1.0)*(2.0*lindex-1.0))))
		    *conj(start[indexlmo])*start[index];

		};
	      
	      lindex=g.ngps_y()-1;
	      index=g.index(rindex,lindex,zindex);
	      indexlmo=g.index(rindex,lindex-1,zindex);
	      tmp+=(complex)(sqrt((lindex*lindex-m*m)/((2.0*lindex+1.0)*(2.0*lindex-1.0))))
		*conj(start[indexlmo])*start[index];

	    };
	  result[rindex]+=(complex)degeneracies[zindex]*tmp;
	};
    };  
  
  return result;
}

wavefunction wavefunction::calculate_Xi(grid g, const fluid &degeneracies, const fluid &ms)
{
  wavefunction result(g.ngps_x());
  long rindex,lindex,index,indexlpt,indexlmt,zindex;
  double r;
  double rr;
  double m;
  complex tmp(0.0,0.0);

  if (g.ngps_y()<3) fprintf(stdout,"warn: calculate_Xi: g.ngps_y()<3\n");

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      result[rindex]=complex(0.0,0.0);
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  m=ms[zindex];

	  lindex=(long)(fabs(m));
	  index=g.index(rindex,lindex,zindex);
	  indexlpt=g.index(rindex,lindex+2,zindex);
	  tmp=(complex)(3.0/(2.0*(2.0*lindex+3.0))*sqrt((((lindex+1.0)*(lindex+1.0)-m*m)*((lindex+2.0)*(lindex+2.0)-m*m))/((2.0*lindex+1)*(2.0*lindex+5))))*conj(start[indexlpt])*start[index];
	  
	  lindex=(long)(fabs(m))+1;
	  index=g.index(rindex,lindex,zindex);
	  indexlpt=g.index(rindex,lindex+2,zindex);
	  tmp+=(complex)(3.0/(2.0*(2.0*lindex+3.0))*sqrt((((lindex+1.0)*(lindex+1.0)-m*m)*((lindex+2.0)*(lindex+2.0)-m*m))/((2.0*lindex+1)*(2.0*lindex+5))))*conj(start[indexlpt])*start[index];
	  
	  for (lindex=(long)(fabs(m))+2; lindex<g.ngps_y()-2; lindex++)
	    {
	      index=g.index(rindex,lindex,zindex);
	      indexlpt=g.index(rindex,lindex+2,zindex);
	      indexlmt=g.index(rindex,lindex-2,zindex);
	      tmp+=(complex)(3.0/(2.0*(2.0*lindex+3.0))*sqrt((((lindex+1.0)*(lindex+1.0)-m*m)*((lindex+2.0)*(lindex+2.0)-m*m))/((2.0*lindex+1)*(2.0*lindex+5))))*conj(start[indexlpt])*start[index]
		+(complex)(3.0/(2.0*(2.0*lindex-1.0))*sqrt((((lindex-1.0)*(lindex-1.0)-m*m)*(lindex*lindex-m*m))/((2.0*lindex-3.0)*(2.0*lindex+1))))*conj(start[indexlmt])*start[index];
	    };

	  lindex=g.ngps_y()-2;
	  index=g.index(rindex,lindex,zindex);
	  indexlmt=g.index(rindex,lindex-2,zindex);
	  tmp+=(complex)(3.0/(2.0*(2.0*lindex-1.0))*sqrt((((lindex-1.0)*(lindex-1.0)-m*m)*(lindex*lindex-m*m))/((2.0*lindex-3.0)*(2.0*lindex+1))))*conj(start[indexlmt])*start[index];
	      
	  lindex=g.ngps_y()-1;
	  index=g.index(rindex,lindex,zindex);
	  indexlmt=g.index(rindex,lindex-2,zindex);
	  tmp+=(complex)(3.0/(2.0*(2.0*lindex-1.0))*sqrt((((lindex-1.0)*(lindex-1.0)-m*m)*(lindex*lindex-m*m))/((2.0*lindex-3.0)*(2.0*lindex+1))))*conj(start[indexlmt])*start[index];


	  for (lindex=(long)(fabs(m)); lindex<g.ngps_y(); lindex++)
	    {
	      index=g.index(rindex,lindex,zindex);
	      tmp+=(complex)((lindex*(lindex+1.0)-3.0*m*m)/((2.0*lindex-1.0)*(2.0*lindex+3.0)))*conj(start[index])*start[index];
	    };


	  result[rindex]+=(complex)degeneracies[zindex]*tmp;
	};  
    };
  
  return result;
}



fluid wavefunction::calculate_Lambda(grid g, const fluid &degeneracies)
{
  fluid result(g.ngps_x());
  long rindex,lindex,index,zindex;
  double tmp;

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      result[rindex]=0.0;
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  tmp=0.0;
	  for (lindex=0; lindex<g.ngps_y(); lindex++)
	    {
	      index=g.index(rindex,lindex,zindex);
	      tmp+=real(conj(start[index])*start[index]);
	    };
	  result[rindex]+=degeneracies[zindex]*tmp;
	};
    };
      
  return result;
}



wavefunction wavefunction::calculate_hartree_one(grid g, const wavefunction &Thetavector)
{
  wavefunction result(g.ngps_x());
  long rindex,rrindex;
  double r;
  double rr;
  double rm;
  complex upperpart(0.0,0.0);
  complex lowerpart(0.0,0.0);

  rindex=0;
  r=g.r(rindex);
  for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
    {
      rr=g.r(rrindex);
      upperpart+=Thetavector[rrindex]*(complex)r/(rr*rr);
    };
  result[rindex]=upperpart;

  for (rindex=1; rindex<g.ngps_x(); rindex++)
    {      
      r=g.r(rindex);
      rm=g.r(rindex-1);

      lowerpart=lowerpart*(complex)rm*(complex)rm/(r*r)+Thetavector[rindex-1]*(complex)rm/(r*r);
      upperpart=upperpart*(complex)r/rm-Thetavector[rindex-1]*(complex)r/(rm*rm);

      result[rindex]=lowerpart+upperpart;

    };
  
  return result*g.delt_x();
}

wavefunction wavefunction::calculate_hartree_two(grid g, const wavefunction &Xivector)
{
  wavefunction result(g.ngps_x());
  long rindex,rrindex;
  double r;
  double rr;
  double rm;
  complex upperpart(0.0,0.0);
  complex lowerpart(0.0,0.0);

  rindex=0;
  r=g.r(rindex);
  for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
    {
      rr=g.r(rrindex);
      upperpart+=Xivector[rrindex]*(complex)r*r/(rr*rr*rr);
    };
  result[rindex]=upperpart;

  for (rindex=1; rindex<g.ngps_x(); rindex++)
    {      
      r=g.r(rindex);
      rm=g.r(rindex-1);

      lowerpart=lowerpart*(complex)(rm*rm*rm)/(r*r*r)+Xivector[rindex-1]*(complex)(rm*rm)/(r*r*r);
      upperpart=upperpart*(complex)(r*r/(rm*rm))-Xivector[rindex-1]*(complex)(r*r)/(rm*rm*rm);

      result[rindex]=lowerpart+upperpart;

    };
  
  return result*g.delt_x();
}



wavefunction wavefunction::calculate_rruxsigma(grid g, grid q, const fluid &ells, const fluid &degeneracies)
{
  wavefunction result(q.ngps_x());
  wavefunction U_ji(q.ngps_x());

  long rindex,index,i,j,l;
  double clebsch,factor;


  result.nullify();
  

  for (i=0; i<q.ngps_z(); i++)
    {
      for (j=0; j<q.ngps_z(); j++)
	{
	  for (l=(int)(fabs(ells[j]-ells[i])); l<=(int)(ells[j]+ells[i]); l++)
	    {
	      clebsch=CG(ells[j],ells[i],(double)(l),0.0,0.0,0.0);
	      //   std::cout << l << " " << ells[j] << " " << ells[i] << " " << clebsch << std::endl;
	      factor=clebsch*clebsch/(2.0*l+1.0);
	      U_ji=calculate_U_ji_l(g,j,i,l);
	      for (rindex=0; rindex<q.ngps_x(); rindex++)
		{
		  index=q.index(rindex,i,j);
		  result[rindex]-=start[index]*factor*U_ji[rindex]*degeneracies[i]*degeneracies[j]*0.25;
		};
	    };
	};
    };
  
  return result;

}; 


wavefunction wavefunction::calculate_vsxj(grid q, const fluid &Lambdavector, const wavefunction &rruxsigma, int non_inversion)
{
  wavefunction result(q.ngps_z()-non_inversion);

  result.nullify();

  long index,j,rindex;

  for (j=0; j<q.ngps_z()-non_inversion; j++)
    {
      for (rindex=0; rindex<q.ngps_x(); rindex++)
	{
	  index=q.index(rindex,j,j);
	  result[j]+=start[index]/(0.5*Lambdavector[rindex]+1e-100)
	    *rruxsigma[rindex];
	};
    };
  
  return result*q.delt_x();


};

wavefunction wavefunction::calculate_kliav(grid g, grid q, int non_inversion, const fluid &klipot, const fluid &degeneracies)
{
  wavefunction result(q.ngps_z()-non_inversion);
  result.nullify();


  long rindex,index,i;


  for (i=0; i<q.ngps_z()-non_inversion; i++)
    {
      for (rindex=0; rindex<q.ngps_x(); rindex++)
	{
	  index=q.index(rindex,i,i);
	  result[i]+=start[index]*klipot[rindex]*0.5*degeneracies[i];
	};
    };
  
  return result*q.delt_x();


};


wavefunction wavefunction::calculate_uxi(grid g, grid q, int non_inversion, const fluid &ells, const fluid &degeneracies)
{
  wavefunction result(q.ngps_z()-non_inversion);
  wavefunction U_ji(q.ngps_x());

  result.nullify();


  long rindex,index,i,j,l;
  double clebsch,factor;

  

  for (i=0; i<q.ngps_z()-non_inversion; i++)
    {
      for (j=0; j<q.ngps_z(); j++)
	{
	  for (l=(int)(fabs(ells[j]-ells[i])); l<=(int)(ells[j]+ells[i]); l++)
	    {	      
	      clebsch=CG(ells[j],ells[i],(double)(l),0.0,0.0,0.0);
	      factor=clebsch*clebsch/(2.0*l+1.0);
	      U_ji=calculate_U_ji_l(g,j,i,l);
	      for (rindex=0; rindex<q.ngps_x(); rindex++)
		{
		  index=q.index(rindex,i,j);
		  result[i]-=start[index]*factor*U_ji[rindex]*degeneracies[j]*0.5;
		};
	    };
	};
    };
  
  return result*q.delt_x();


};


complex wavefunction::calculate_Mji(grid q, const fluid &Lambdavector, const fluid &ells, const fluid &ms, int j, int i)
{
  complex result(0.0,0.0);

  long rindex,L,index;
  double clebsch_i, clebsch_ii;
  double clebschsum=0.0;

  //   for (L=(int)(fabs(ells[j]-ells[i])); L<=(int)(ells[j]+ells[i]); L++)
  //     {	      
  //       clebsch_i=CG(ells[j],ells[i],(double)(L),0.0,0.0,0.0);
  //       clebsch_ii=CG(ells[j],ells[i],(double)(L),-ms[j],ms[i],ms[i]-ms[j]);
  //       clebschsum+=(2.0*ells[j]+1.0)*(2.0*ells[i]+1.0)/(2.0*L+1.0)
  //       	*clebsch_i*clebsch_i*clebsch_ii*clebsch_ii;
  //       for (rindex=0; rindex<q.ngps_x(); rindex++)
  // 	{
  // 	  index=q.index(rindex,j,i);
  // 	  result+=start[index]*conj(start[index])/(0.5*Lambdavector[rindex]+1e-100)
  // 	    *(2.0*ells[j]+1.0)*(2.0*ells[i]+1.0)/(2.0*L+1.0)
  // 	    *clebsch_i*clebsch_i*clebsch_ii*clebsch_ii;
  // 	};
  //     };

  //   for (L=(int)(fabs(ells[j]-ells[i])); L<=(int)(ells[j]+ells[i]); L++)
  //     {	      
  //       clebsch_i=CG(ells[j],ells[i],(double)(L),0.0,0.0,0.0);
  //       clebschsum+=clebsch_i*clebsch_i;
  //       for (rindex=0; rindex<q.ngps_x(); rindex++)
  // 	{
  // 	  index=q.index(rindex,j,i);
  // 	  result+=start[index]*conj(start[index])/(0.5*Lambdavector[rindex]+1e-100)
  // 	    *clebsch_i*clebsch_i;
  // 	};
  //     };

  //   std::cout << j << ", " << i << " -> clebschsum: " << clebschsum << std::endl;

  for (rindex=0; rindex<q.ngps_x(); rindex++)
    {
      index=q.index(rindex,j,i);
      result+=start[index]*conj(start[index])/(0.5*Lambdavector[rindex]+1e-100);
    };




  return result*q.delt_x();

}




wavefunction wavefunction::calculate_U_ji_l(grid g, long j, long i, long order)
{

  grid q;
  q.set_dim(34); // all the same
  q.set_ngps(g.ngps_x(),g.ngps_z(),g.ngps_z());
  q.set_delt(g.delt_x(),1.0,1.0);
  q.set_offs(0,0,0);


  wavefunction result(g.ngps_x());
  wavefunction stuff(g.ngps_x());
  long rindex,rrindex, k, index;
  double r;
  double rr;
  double rm;
  complex upperpart;
  complex lowerpart;


  stuff=extract_x(q,j,i);


  upperpart=complex(0.0,0.0);
  lowerpart=complex(0.0,0.0);

  switch (order)
    {
    case 0 : // avoiding pow ...
      rindex=0;
      r=g.r(rindex);
      for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
	{
	  rr=g.r(rrindex);
	  upperpart+=stuff[rrindex]*1.0/rr;
	};
      result[rindex]=upperpart;
      
      for (rindex=1; rindex<g.ngps_x(); rindex++)
	{      
	  r=g.r(rindex);
	  rm=g.r(rindex-1);
	  
	  lowerpart=lowerpart*rm/r
	    +stuff[rindex-1]*1.0/r;
	  upperpart=upperpart
	    -stuff[rindex-1]*1.0/rm;
	  
	  result[rindex]=lowerpart+upperpart;
	  
	};
      break;
    case 1 : // avoiding pow ...
      rindex=0;
      r=g.r(rindex);
      for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
	{
	  rr=g.r(rrindex);
	  upperpart+=stuff[rrindex]*r/(rr*rr);
	};
      result[rindex]=upperpart;
      
      for (rindex=1; rindex<g.ngps_x(); rindex++)
	{      
	  r=g.r(rindex);
	  rm=g.r(rindex-1);
	  
	  lowerpart=lowerpart*(rm*rm)/(r*r)
	    +stuff[rindex-1]*rm/(r*r);
	  upperpart=upperpart*r/rm
	    -stuff[rindex-1]*r/(rm*rm);
	  
	  result[rindex]=lowerpart+upperpart;
	  
	};
      break;
    case 2 :  // avoiding pow ...
      rindex=0;
      r=g.r(rindex);
      for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
	{
	  rr=g.r(rrindex);
	  upperpart+=stuff[rrindex]*r*r/(rr*rr*rr);
	};
      result[rindex]=upperpart;
      
      for (rindex=1; rindex<g.ngps_x(); rindex++)
	{      
	  r=g.r(rindex);
	  rm=g.r(rindex-1);
	  
	  lowerpart=lowerpart*(rm*rm*rm)/(r*r*r)
	    +stuff[rindex-1]*rm*rm/(r*r*r);
	  upperpart=upperpart*r*r/(rm*rm)
	    -stuff[rindex-1]*r*r/(rm*rm*rm);
	  
	  result[rindex]=lowerpart+upperpart;
	  
	};
      break;
    case 3 :  // avoiding pow ...
      rindex=0;
      r=g.r(rindex);
      for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
	{
	  rr=g.r(rrindex);
	  upperpart+=stuff[rrindex]*r*r*r/(rr*rr*rr*rr);
	};
      result[rindex]=upperpart;
      
      for (rindex=1; rindex<g.ngps_x(); rindex++)
	{      
	  r=g.r(rindex);
	  rm=g.r(rindex-1);
	  
	  lowerpart=lowerpart*(rm*rm*rm*rm)/(r*r*r*r)
	    +stuff[rindex-1]*rm*rm*rm/(r*r*r*r);
	  upperpart=upperpart*r*r*r/(rm*rm*rm)
	    -stuff[rindex-1]*r*r*r/(rm*rm*rm*rm);
	  
	  result[rindex]=lowerpart+upperpart;
	  
	};
      break;
    case 4 :  // avoiding pow ...
      rindex=0;
      r=g.r(rindex);
      for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
	{
	  rr=g.r(rrindex);
	  upperpart+=stuff[rrindex]*r*r*r*r/(rr*rr*rr*rr*rr);
	};
      result[rindex]=upperpart;
      
      for (rindex=1; rindex<g.ngps_x(); rindex++)
	{      
	  r=g.r(rindex);
	  rm=g.r(rindex-1);
	  
	  lowerpart=lowerpart*(rm*rm*rm*rm*rm)/(r*r*r*r*r)
	    +stuff[rindex-1]*rm*rm*rm*rm/(r*r*r*r*r);
	  upperpart=upperpart*r*r*r*r/(rm*rm*rm*rm)
	    -stuff[rindex-1]*r*r*r*r/(rm*rm*rm*rm*rm);
	  
	  result[rindex]=lowerpart+upperpart;
	  
	};
      break;
    default :
      std::cout << "I AM IN THE U_IJ LOOP WITH POW" << std::endl;
      rindex=0;
      r=g.r(rindex);
      for (rrindex=0; rrindex<g.ngps_x(); rrindex++)
	{
	  rr=g.r(rrindex);
	  upperpart+=stuff[rrindex]*pow(r,order)/pow(rr,order+1);
	};
      result[rindex]=upperpart;
      
      for (rindex=1; rindex<g.ngps_x(); rindex++)
	{      
	  r=g.r(rindex);
	  rm=g.r(rindex-1);
	  
	  lowerpart=lowerpart*pow(rm/r,order+1)
	    +stuff[rindex-1]*pow(rm/r,order)/r;
	  upperpart=upperpart*pow(r/rm,order)
	    -stuff[rindex-1]*pow(r/rm,order)/rm;

	  
	  result[rindex]=lowerpart+upperpart;
	  
	};
	  
    };

  return result*g.delt_x();

}



wavefunction wavefunction::calculate_Lambda_ki(grid g, const fluid &ells)
{

  grid q;
  q.set_dim(34); // all the same
  q.set_ngps(g.ngps_x(),g.ngps_z(),g.ngps_z());
  q.set_delt(g.delt_x(),1.0,1.0);
  q.set_offs(0,0,0);


  wavefunction result(g.ngps_x()*g.ngps_z()*g.ngps_z());
  long jott,kah,rindex,lindex,index,wfindexjott,wfindexkah;

  result.nullify();


  for (jott=0; jott<g.ngps_z(); jott++)
    {
      for (kah=0; kah<g.ngps_z(); kah++)
	{
	  for (rindex=0; rindex<g.ngps_x(); rindex++)
	    {
	      index=q.index(rindex,kah,jott);
	      wfindexjott=g.index(rindex,(long)(ells[jott]),jott);
	      wfindexkah=g.index(rindex,(long)(ells[kah]),kah);
	      result[index]=conj(start[wfindexjott])*start[wfindexkah];
	    };
	};
    };
  
  return result;

}





fluid wavefunction::calculate_x_zero(grid g, const fluid &Lambdavector)
{
  fluid result(g.ngps_x());
  long rindex;
  double r;

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      r=g.r(rindex);
      result[rindex]=-pow(0.075990887*Lambdavector[rindex]/(r*r),0.3333333333);
    };

  return result;
}

fluid wavefunction::calculate_c_zero(grid g, const fluid &Lambdavector)
{
  fluid result(g.ngps_x());
  long rindex;
  double r,rs,sqrtrs;
  double A=0.031091;
  double betaone=7.5957;
  double betatwo=3.5876;
  double betathree=1.6382;
  double betafour=0.49294;
  double alphaone=0.2137;

  double Qone, Qzero, Qonedash;

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      r=g.r(rindex);
      rs=pow(3.0*r*r/(Lambdavector[rindex]+1e-100),0.33333333333);
      sqrtrs=sqrt(rs);
      Qone=2.0*A*(betaone*sqrt(rs)+betatwo*rs+betathree*sqrtrs*sqrtrs*sqrtrs+betafour*rs*rs);
      Qzero=-2.0*A*(1.0+alphaone*rs);
      Qonedash=A*(betaone/sqrtrs+2.0*betatwo+3.0*betathree*sqrtrs+4.0*betafour*rs);
      result[rindex]=-2.0*A*log(1.0+1.0/Qone)*(1.0-0.666666666667*alphaone*rs)
	+0.33333333333*rs*Qzero*Qonedash/(Qone*Qone+Qone);
    };

  return result;
}



wavefunction wavefunction::calculate_x_one(grid g, const fluid &Lambdavector, const wavefunction &Thetavector)
{
  wavefunction result(g.ngps_x());
  long rindex;
  double r;
  double prefact=pow(3.0/(4.0*M_PI*M_PI),0.33333333333);


  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {      
      r=g.r(rindex);
      result[rindex]=-prefact*Thetavector[rindex]/(1e-100+pow(r*Lambdavector[rindex],0.666666667));
    };
  
  return result;
}




// this is 4 pi r^2 n(r)
fluid wavefunction::calculate_radial_density(grid g, const fluid &degeneracies)
{
  fluid result(g.ngps_x());
  long rindex,zindex,lindex,index;
  double r,deg;


  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      r=g.r(rindex);
      result[rindex]=0.0;
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  deg=degeneracies[zindex];
	  for (lindex=0; lindex<g.ngps_y(); lindex++)
	    {
	      index=g.index(rindex,lindex,zindex);
	      result[rindex]+=real(conj(start[index])*start[index])*deg;
	    };
	};
      //      result[rindex]=result[rindex];
    };

  return result;
}


double wavefunction::totalenergy_sic(grid g, const fluid &degeneracies, const fluid &ms)
{
  double result=0.0;

  double tmp_totenerg_x_lda, tmp_hartree_energ;

  // the grid on which Lambda_i(r,t)
  grid q;
  q.set_dim(34); // all the same
  q.set_ngps(g.ngps_x(),1,g.ngps_z());
  q.set_delt(g.delt_x(),1.0,1.0);
  q.set_offs(0,0,0);


  fluid Lambda(g.ngps_x()*g.ngps_z());
  fluid vi(g.ngps_x()*g.ngps_z());
  fluid tmp_hartree_pot(g.ngps_x());
  fluid tmp(g.ngps_x());

  Lambda.nullify();
  vi.nullify();

  long rindex,iindex,index,lindex,wfindex,jindex,iiindex,jjindex;

  // calculate Lambda_i
  for (iindex=0; iindex<g.ngps_z(); iindex++)
    {
      for (rindex=0; rindex<g.ngps_x(); rindex++)
	{
	  index=q.index(rindex,0,iindex);
	  for (lindex=(long)(fabs(ms[iindex])); lindex<g.ngps_y(); lindex++)
	    {
	      wfindex=g.index(rindex,lindex,iindex);
	      Lambda[index]+=real(conj(start[wfindex])*start[wfindex]);
	    };
	};
    };
  
  for (iindex=0; iindex<g.ngps_z(); iindex++)
    {
      tmp=Lambda.extract_x(q,0,iindex);
      tmp_hartree_pot=calculate_hartree_zero(g,tmp);
      tmp_totenerg_x_lda=1.25992105*totalenergy_x_lda(g,tmp);
      tmp_hartree_energ=totalenergy_hartree(g,tmp,tmp_hartree_pot);
      result-=degeneracies[iindex]*(tmp_totenerg_x_lda+tmp_hartree_energ);
    };

  
  return result;


};


double wavefunction::totalenergy_x_lda(grid g, const fluid &Lambdavector)
{
  double result=0.0;
  double r;
  long rindex;

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      r=g.r(rindex);
      result+=pow(Lambdavector[rindex]*Lambdavector[rindex]/r,0.6666666667);
    };

  return -g.delt_x()*0.317674071*result;

}


complex wavefunction::totalenergy_exact_x(grid g,
					  const fluid &ells, const fluid &degeneracies)
{
  complex result=(0.0,0.0);
  long j;

  grid q;
  q.set_dim(34); // all the same
  q.set_ngps(g.ngps_x(),g.ngps_z(),g.ngps_z());
  q.set_delt(g.delt_x(),1.0,1.0);
  q.set_offs(0,0,0);

  wavefunction Lambda(g.ngps_x()*g.ngps_z()*g.ngps_z());  
  wavefunction uxi(g.ngps_z());

  Lambda=calculate_Lambda_ki(g,ells);
  uxi=Lambda.calculate_uxi(g,q,0,ells,degeneracies);

  for (j=0; j<g.ngps_z(); j++)
    result+=uxi[j]*0.5*degeneracies[j];

  return result;
}


double wavefunction::totalenergy_c_lda(grid g, const fluid &Lambdavector)
{
  double result=0.0;
  double r,rs,sqrtrs;
  long rindex;
  double A=0.031091;
  double betaone=7.5957;
  double betatwo=3.5876;
  double betathree=1.6382;
  double betafour=0.49294;
  double alphaone=0.2137;

  double Qone;

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    {
      r=g.r(rindex);
      rs=pow(3.0*r*r/(Lambdavector[rindex]+1e-100),0.33333333333);
      sqrtrs=sqrt(rs);
      Qone=2.0*A*(betaone*sqrt(rs)+betatwo*rs+betathree*sqrtrs*sqrtrs*sqrtrs+betafour*rs*rs);
      result+= -Lambdavector[rindex]*2.0*A*(1.0+alphaone*rs)*log(1.0+1.0/Qone);
    };

  return g.delt_x()*result;

}

double wavefunction::totalenergy_hartree(grid g,
					 const fluid &Lambdavector, const fluid &hartree_zero)
{
  double result=0.0;
  double r;
  long rindex;

  for (rindex=0; rindex<g.ngps_x(); rindex++)
    result+=Lambdavector[rindex]*hartree_zero[rindex];

  return g.delt_x()*0.5*result;
}


double wavefunction::totalenergy_single_part(grid g, const wavefunction &orb_energs, const fluid &degeneracies)
{
  double result=0.0;
  long zindex;

  for (zindex=0; zindex<g.ngps_z(); zindex++)
    {
      result+=degeneracies[zindex]*real(orb_energs[zindex]);
    };

  return result;

};


void wavefunction::calculate_staticpot(grid g, hamop hamil)
{  
  double x,y,z,time;
  int me;
  long index,xindex,yindex,ell,m,zindex;
  complex imagi(0.0,1.0);

  me = 0;
  time = 0.0;

  switch (g.dimens())
    {
    case 34 :
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{ 
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      x=g.r(xindex);
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  start[index]=hamil.scalarpot(x,yindex,0.0,0.0,me)
		    +0.5*yindex*(yindex+1)/(x*x)
		    -imagi*(complex)hamil.imagpot(xindex,yindex,0,time,g);
		};
	    };
	};
      break;
    case 44 :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  x=g.r(xindex);
	  for (ell=0; ell<g.ngps_y(); ell++)
	    { 
	      for (m=-ell; m<=ell; m++)
		{
		  index=g.rlmindex(xindex,ell,m);
		  start[index]=hamil.scalarpot(x,ell,0.0,0.0,me)+
		    0.5*ell*(ell+1)/(x*x)
		    -imagi*(complex)hamil.imagpot(xindex,ell,0,time,g);
		};
	    };
	};
      break;

    default :
      fprintf(stderr, "err: unknown propagation mode (g.dimens())\n");
      exit(-53);
    };


}



void wavefunction::realific()
{
  for(long i=0; i<wf_dim; i++)
    start[i]=complex(real(start[i]),0.0);
}



/*! \fn complex wavefunction::energy(double time, grid g, hamop hamil,int me,
  const wavefunction &staticpot, double charge)
   
  \param time
  \param grid
  \param hamil
  \param me
  \param staticpot
  \param charge
   
  \return Kinetic plus potential energy expectation for the wavefunction.
*/
complex wavefunction::energy(double time, grid g, hamop hamil,int me,
			     const wavefunction &staticpot, double charge)
{
  double Eff_Charge = charge;

  complex result(0.0,0.0);
  double const_diag_part;
  long xindex, yindex, zindex,m;
  double x,y,z,r,field,halfvecpotvecpot;
  complex vpx, vpy, vpz;
  long index, index_xp, index_xm, index_yp, index_ym, index_zp, index_zm;
  double b_upperleft, M_two_upperleft, Delta_two_upperleft;
  wavefunction tmp(g.ngps_x());	
  wavefunction tmp_two(g.ngps_x());	
  wavefunction rhsone(g.size());
  double aaaa, bbbb, cccc;

  switch (g.dimens())
    {
    case 34 :
      aaaa=1.0/(g.delt_x()*g.delt_x());
      cccc=aaaa;
      bbbb=-2.0/(g.delt_x()*g.delt_x());
      Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
      M_two_upperleft=-2.0*(1.0+g.delt_x()*g.delt_x()/12.0*Delta_two_upperleft);

      for (zindex=0; zindex<g.ngps_z(); zindex++) {
	// apply (Delta_2 + M_2 V_eff,l) to Phi
	for (yindex=0; yindex<g.ngps_y(); yindex++) { 
	  xindex=0;
	  index = g.index(xindex,yindex,zindex);
	  index_xp = g.index(xindex+1,yindex,zindex);
	  if (yindex==0) {
	    rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
	      +((complex)M_two_upperleft*(staticpot[index])+Delta_two_upperleft)*start[index];
	  }
	  else {
	    rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
	      +(bbbb-(complex)FOT*(staticpot[index]))*start[index];
	  };
	  
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++) { 
	    index=g.index(xindex,yindex,zindex);
	    index_xp=g.index(xindex+1,yindex,zindex);
	    index_xm=g.index(xindex-1,yindex,zindex);
	    rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
	      +(bbbb-(complex)FOT*(staticpot[index]))*start[index]
	      +(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	  };
	  
	  
	  xindex=g.ngps_x()-1;
	  index=g.index(xindex,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  rhsone[index]=(bbbb-(complex)FOT*(staticpot[index]))*start[index]
	    +(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	};
	
	// solve M_2^-1 times that stuff
	aaaa=-OOS;
	cccc=aaaa;
	bbbb=-FOT;
	
	for (yindex=0; yindex<g.ngps_y(); yindex++) { 
	  xindex=0;
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  
	  if (yindex==0) {
	    b_upperleft=M_two_upperleft;
	  }
	  else {
	    b_upperleft=-FOT;
	  };
	  
	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++) { 
	    index=g.index(xindex,yindex,zindex);
	    tmp_two[xindex]=rhsone[index];
	  };  
	  
	  tmp.solve_toep(aaaa,bbbb,b_upperleft,bbbb,cccc,tmp_two,g.ngps_x());
	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++) { 
	    index=g.index(xindex,yindex,zindex);
	    rhsone[index]=tmp[xindex];
	  };  
	  
	};
	
	
      };
      
      // and finally do the integration 
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  result=result+conj(start[index])*rhsone[index]*(complex)g.delt_x();
		}
	    }
	}
      
      
      break;
    case 44 :
      aaaa=1.0/(g.delt_x()*g.delt_x());
      cccc=aaaa;
      bbbb=-2.0/(g.delt_x()*g.delt_x());
      Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
      M_two_upperleft=-2.0*(1.0+g.delt_x()*g.delt_x()/12.0*Delta_two_upperleft);

      // apply (Delta_2 + M_2 V_eff,l) to Phi
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  for (m=-yindex; m<=yindex; m++)
	    {
	      xindex=0;
	      index=g.rlmindex(xindex,yindex,m);
	      index_xp=g.rlmindex(xindex+1,yindex,m);
	      if (yindex==0)
		{
		  rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		    +((complex)M_two_upperleft*(staticpot[index])+Delta_two_upperleft)*start[index];
		}
	      else
		{
		  rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		    +(bbbb-(complex)FOT*(staticpot[index]))*start[index];
		};
	  
	      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
		{ 
		  index=g.rlmindex(xindex,yindex,m);
		  index_xp=g.rlmindex(xindex+1,yindex,m);
		  index_xm=g.rlmindex(xindex-1,yindex,m);
		  rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		    +(bbbb-(complex)FOT*(staticpot[index]))*start[index]
		    +(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
		};
	  
	  
	      xindex=g.ngps_x()-1;
	      index=g.rlmindex(xindex,yindex,m);
	      index_xm=g.rlmindex(xindex-1,yindex,m);
	      rhsone[index]=(bbbb-(complex)FOT*(staticpot[index]))*start[index]
		+(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	    };
	};
      
      // solve M_2^-1 times that stuff
      aaaa=-OOS;
      cccc=aaaa;
      bbbb=-FOT;
      
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  for (m=-yindex; m<=yindex; m++)
	    {
	      xindex=0;
	      index=g.rlmindex(xindex,yindex,m);
	      index_xp=g.rlmindex(xindex+1,yindex,m);
	      
	      if (yindex==0)
		{
		  b_upperleft=M_two_upperleft;
		}
	      else
		{
		  b_upperleft=-FOT;
		};
	      
	      
	      for (xindex=0; xindex<g.ngps_x(); xindex++)
		{ 
		  index=g.rlmindex(xindex,yindex,m);
		  tmp_two[xindex]=rhsone[index];
		};  
	      
	      tmp.solve_toep(aaaa,bbbb,b_upperleft,bbbb,cccc,tmp_two,g.ngps_x());
	  
	      for (xindex=0; xindex<g.ngps_x(); xindex++)
		{ 
		  index=g.rlmindex(xindex,yindex,m);
		  rhsone[index]=tmp[xindex];
		};  
	      
	    };
	};
      
      
      // and finally do the integration 
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      for (m=-yindex; m<=yindex; m++)
		{
		  index=g.rlmindex(xindex,yindex,m);
		  result=result+conj(start[index])*rhsone[index]*(complex)g.delt_x();
		};
	    };
	};
      
      break;
    };
  
  return result;
}


wavefunction wavefunction::orbital_energies(double time, grid g, 
					    hamop hamil, int me,  const wavefunction &staticpot, double charge)
{
  double Eff_Charge = charge;

  wavefunction result(g.ngps_z());
  double const_diag_part;
  long xindex, yindex, zindex, orb_no;
  double x,y,z,r,field,halfvecpotvecpot;
  complex vpx, vpy, vpz;
  long index, index_xp, index_xm, index_yp, index_ym, index_zp, index_zm;
  double b_upperleft, M_two_upperleft, Delta_two_upperleft;
  wavefunction tmp(g.ngps_x());	
  wavefunction tmp_two(g.ngps_x());	
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double aaaa, bbbb, cccc;

  switch (g.dimens())
    {
    case 34 :
      for (orb_no=0; orb_no<g.ngps_z(); orb_no++)
	{
	  zindex=orb_no;
	  aaaa=1.0/(g.delt_x()*g.delt_x());
	  cccc=aaaa;
	  bbbb=-2.0/(g.delt_x()*g.delt_x());
	  Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
	  M_two_upperleft=-2.0*(1.0+g.delt_x()*g.delt_x()/12.0*Delta_two_upperleft);
	  
	  // apply (Delta_2 + M_2 V_eff,l) to Phi
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    { 
	      xindex=0;
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      if (yindex==0)
		{
		  rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		    +(M_two_upperleft*(staticpot[index])+Delta_two_upperleft)*start[index];
		}
	      else
		{
		  rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		    +(bbbb-FOT*(staticpot[index]))*start[index];
		};
	      
	      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  index_xp=g.index(xindex+1,yindex,zindex);
		  index_xm=g.index(xindex-1,yindex,zindex);
		  rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		    +(bbbb-FOT*(staticpot[index]))*start[index]
		    +(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
		};
	      
	      
	      xindex=g.ngps_x()-1;
	      index=g.index(xindex,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      rhsone[index]=(bbbb-FOT*(staticpot[index]))*start[index]
		+(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	    };
	  
	  // solve M_2^-1 times that stuff
	  aaaa=-OOS;
	  cccc=aaaa;
	  bbbb=-FOT;
	  
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    { 
	      xindex=0;
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      
	      if (yindex==0)
		{
		  b_upperleft=M_two_upperleft;
		}
	      else
		{
		  b_upperleft=-FOT;
		};
	      
	      
	      for (xindex=0; xindex<g.ngps_x(); xindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  tmp_two[xindex]=rhsone[index];
		};  
	      
	      tmp.solve_toep(aaaa,bbbb,b_upperleft,bbbb,cccc,tmp_two,g.ngps_x());
	      
	      for (xindex=0; xindex<g.ngps_x(); xindex++)
		{ 
		  index=g.index(xindex,yindex,zindex);
		  rhsone[index]=tmp[xindex];
		};  
	      
	    };
	  
	  result[zindex]=complex(0.0,0.0);
	  // and finally do the integration 
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  result[zindex]=result[zindex]+conj(start[index])*rhsone[index]*g.delt_x();
		};
	    };
	};
      break;
    };

  return result;
}

complex wavefunction::orbital_energy(double time, grid g, hamop hamil,
				     int me,  const wavefunction &staticpot, double charge, long orb_no)
{
  double Eff_Charge = charge;

  complex result(0.0,0.0);
  double const_diag_part;
  long xindex, yindex, zindex;
  double x,y,z,r,field,halfvecpotvecpot;
  complex vpx, vpy, vpz;
  long index, index_xp, index_xm, index_yp, index_ym, index_zp, index_zm;
  double b_upperleft, M_two_upperleft, Delta_two_upperleft;
  wavefunction tmp(g.ngps_x());	
  wavefunction tmp_two(g.ngps_x());	
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double aaaa, bbbb, cccc;

  switch (g.dimens())
    {
    case 34 :
      zindex=orb_no;
      aaaa=1.0/(g.delt_x()*g.delt_x());
      cccc=aaaa;
      bbbb=-2.0/(g.delt_x()*g.delt_x());
      Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
      M_two_upperleft=-2.0*(1.0+g.delt_x()*g.delt_x()/12.0*Delta_two_upperleft);
	  
      // apply (Delta_2 + M_2 V_eff,l) to Phi
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  xindex=0;
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  if (yindex==0)
	    {
	      rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		+(M_two_upperleft*(staticpot[index])+Delta_two_upperleft)*start[index];
	    }
	  else
	    {
	      rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		+(bbbb-FOT*(staticpot[index]))*start[index];
	    };
	      
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    { 
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		+(bbbb-FOT*(staticpot[index]))*start[index]
		+(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	    };
	      
	      
	  xindex=g.ngps_x()-1;
	  index=g.index(xindex,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  rhsone[index]=(bbbb-FOT*(staticpot[index]))*start[index]
	    +(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	};
	  
      // solve M_2^-1 times that stuff
      aaaa=-OOS;
      cccc=aaaa;
      bbbb=-FOT;
	  
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  xindex=0;
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	      
	  if (yindex==0)
	    {
	      b_upperleft=M_two_upperleft;
	    }
	  else
	    {
	      b_upperleft=-FOT;
	    };
	      
	      
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      index=g.index(xindex,yindex,zindex);
	      tmp_two[xindex]=rhsone[index];
	    };  
	      
	  tmp.solve_toep(aaaa,bbbb,b_upperleft,bbbb,cccc,tmp_two,g.ngps_x());
	      
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      index=g.index(xindex,yindex,zindex);
	      rhsone[index]=tmp[xindex];
	    };  
	      
	};
	  
	  
      // and finally do the integration 
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      result=result+conj(start[index])*rhsone[index]*g.delt_x();
	    };
	};
	  
      break;
    };

  return result;
}






wavefunction wavefunction::orbital_hartrees(double time, grid g, int me, 
					    const fluid &wf_one)
{
  wavefunction result(g.ngps_z());
  long xindex, yindex, zindex, indexlm, indexlp, orb_no;
  double x,y,z,r;
  long index;

  

  switch (g.dimens())
    {
    case 34 :
      for (orb_no=0; orb_no<g.ngps_z(); orb_no++)
	{
	  zindex=orb_no;
	  result[zindex]=complex(0.0,0.0);
	  // the part V_ee^0
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  result[zindex]=result[zindex]+conj(start[index])
		    *wf_one[xindex]*start[index]*g.delt_x();
		};
	    };
	};
      break;
    };

  return result;
}

complex wavefunction::orbital_hartree(double time, grid g, int me, 
				      const fluid &wf_one, long orb_no)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex, indexlm, indexlp;
  double x,y,z,r;
  long index;

  

  switch (g.dimens())
    {
    case 34 :
      zindex=orb_no;
      // the part V_ee^0
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      result=result+conj(start[index])*
		wf_one[xindex]*start[index]*g.delt_x();
	    };
	};
      break;
    };

  return result;
}


fluid wavefunction::orbital_norms(grid g)
{
  fluid result(g.ngps_z());
  long orb_no;
  wavefunction extr_orbital(g.ngps_x()*g.ngps_y());

  for (orb_no=0; orb_no<g.ngps_z(); orb_no++)
    {
      extr_orbital=extract_xy(g,orb_no);
      result[orb_no]=abs((extr_orbital*extr_orbital)*g.delt_x());
    };
  
  return result;

}

// overloaded below --- this does not work in rlm-mode !!!
wavefunction wavefunction::project(grid g, grid gorig, wavefunction &orig)
{
  wavefunction result(gorig.ngps_z());
  wavefunction origextract(gorig.ngps_x()*gorig.ngps_y());
  wavefunction origextractregridded(g.ngps_x()*g.ngps_y());
  long zindex;
  
  if (g.delt_x() != gorig.delt_x())
    {
      std::cerr << "err: unequal delta r in wavefunction::project()." << std::endl;
      exit(-100);
    }

  if (gorig.size()!=orig.wf_size()) 
    {
      std::cerr << "err: gorig.size()!=orig.wf_size() wavefunction::project()."<<std::endl;
      fprintf( stderr, "original_grid_size: %ld, orignial_wf_size: %ld\n", gorig.size(), orig.wf_size() );  
      exit(-101);
    }

  for (zindex=0; zindex<gorig.ngps_z(); zindex++)
    {
      origextract=orig.extract_xy(gorig,zindex);
      origextractregridded.regrid(g,gorig,origextract);
      result[zindex]=origextractregridded*(*this);
      result[zindex]*=g.delt_x();
    }
  return result;
}


// overloaded above --- this is for only one orbital (zindex) and does work in rlm-mode
complex wavefunction::project(grid g, grid gorig, wavefunction &orig, long zindex)
{
  complex result(0.0,0.0);
  long xindex,yindex,xupper,yupper;
  long index,iindex,mindex,m_limit;
  
  if ((g.ngps_z()>1) || (gorig.ngps_z()>1))
    {
      std::cerr << "err: attempt to call complex wavefunction::project() with ngpsz>1" << std::endl;
      exit(-100);
    }
  if (g.delt_x() != gorig.delt_x())
    {
      std::cerr << "err: unequal delta r in complex wavefunction::project()." << std::endl;
      exit(-100);
    }
  if (g.dimens() != gorig.dimens())
    {
      std::cerr << "err: unequal propagation modes for the two grids in complex wavefunction::project()." << std::endl;
      exit(-100);
    }


  if (gorig.ngps_x() > g.ngps_x()){xupper=g.ngps_x();}
  else {xupper=gorig.ngps_x();}
  
  if (gorig.ngps_y() > g.ngps_y()) {yupper=g.ngps_y();}
  else {yupper=gorig.ngps_y();}
  
  for (xindex=0; xindex<xupper; xindex++)
    {
      for (yindex=0; yindex<yupper; yindex++)
	{   
	  if (g.dimens()==34) {m_limit=0;} else {m_limit=yindex;}    
	  for (mindex=-m_limit; mindex<=m_limit; mindex++)
	    {
	      index=g.index(xindex,yindex,mindex,zindex);
	      iindex=gorig.index(xindex,yindex,mindex,zindex);
	      result+=conj(orig[iindex])*start[index];
	    };
	};
    };
  

  return result*g.delt_x();
}



// simple case with just one orbital; overloaded below
void wavefunction::normalize(grid g)
{
  (*this)*=1.0/sqrt(norm(g));;
}



// KS case for orbitals with different quantum numbers m 
// overloaded above
void wavefunction::normalize(grid g, const fluid &ms)
{
  long xindex, yindex, zindex, zzindex;
  wavefunction orbital_one(g.ngps_x()*g.ngps_y());
  wavefunction orbital_two(g.ngps_x()*g.ngps_y());
  complex scalarprod;
  grid gg;
  gg.set_dim(34);
  gg.set_ngps(g.ngps_x(),g.ngps_y(),1);
  gg.set_delt(g.delt_x(),1,1);
  gg.set_offs(0,0,0);
  
  switch (g.dimens())
    {
    case 34 :
      // *** normalize the first orbital
      orbital_one=extract_xy(g,0);
      orbital_one*=1.0/sqrt(orbital_one.norm(gg));
      embed_as_xy(g,0,orbital_one);
      for (zindex=0; zindex<g.ngps_z()-1; zindex++)
	{
	  // *** take orbital no. i
	  orbital_one=extract_xy(g,zindex);
	  for (zzindex=zindex+1; zzindex<g.ngps_z(); zzindex++)
	    {  
	      // *** take orbital no. j>i and normalize it
	      orbital_two=extract_xy(g,zzindex);
	      orbital_two*=1.0/sqrt(orbital_two.norm(gg));
	      // *** if the m quantum numbers are different
	      // they're orthogonal anyway
	      if (ms[zindex]==ms[zzindex])
		{
		  // *** project orbital i out of orbital j,
		  // normalize the result and store it
		  // *** as the  n e w  orbital j
		  scalarprod=orbital_one*orbital_two;
		  orbital_two=orbital_two-scalarprod*g.delt_x()*orbital_one;
		  orbital_two*=1.0/sqrt(orbital_two.norm(gg));
		}
	      embed_as_xy(g,zzindex,orbital_two);
	    }
	}

      break;
    default :
      (*this)*=1.0/sqrt((*this).norm(g));
    };
}

double wavefunction::total_ks_norm(grid g, const fluid &degeneracies)
{
  double result=0.0;
  long xindex, yindex, zindex,index;
  double tmp;

  for (zindex=0; zindex<g.ngps_z(); zindex++)
    { 
      tmp=0.0;
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      tmp+=real(conj(start[index])*start[index]);
	    }
	}
      result+=tmp*degeneracies[zindex];
    }

  result=result*g.delt_x();

  return result;
};

double wavefunction::total_ks_norm_in_sphere(grid g, const fluid &degeneracies, long xlimit)
{
  double result=0.0;
  long xindex, yindex, zindex,index;
  double tmp;

  for (zindex=0; zindex<g.ngps_z(); zindex++)
    { 
      tmp=0.0;
      for (xindex=0; xindex<xlimit; xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      tmp+=real(conj(start[index])*start[index]);
	    }
	}
      result+=tmp*degeneracies[zindex];
    }

  result=result*g.delt_x();

  return result;
}

double wavefunction::total_ks_norm_in_shell(grid g, const fluid &degeneracies, long xlowerlimit, long xupperlimit)
{
  double result=0.0;
  long xindex, yindex, zindex,index;
  double tmp;

  for (zindex=0; zindex<g.ngps_z(); zindex++)
    { 
      tmp=0.0;
      for (xindex=xlowerlimit; xindex<xupperlimit; xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      tmp+=real(conj(start[index])*start[index]);
	    };
	};
      result+=tmp*degeneracies[zindex];
    };

  result=result*g.delt_x();

  return result;
};

/*!

 */
double wavefunction::norm(grid g)
{
  double result=0.0;
  long xindex, yindex, zindex,m;
  long index,indexm,indexmm;

  switch (g.dimens())
    {
    case 34 :

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+real(conj(start[index])*start[index])*g.delt_x();
	    }
	}
      break;

    case 44 :

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      for (m=-yindex; m<=yindex; m++)
		{
		  index=g.rlmindex(xindex,yindex,m);
		  result=result+real(conj(start[index])*start[index])*g.delt_x();
		}
	    }
	}
      break;

    default :
      fprintf(stderr,"err: unknown propagation mode (g.dimens())\n");
      exit(-51);
    }

  //   std::cout << *start << std::endl;
  //   std::cout << "wf(0,0,0)" << start[g.index(0,0,0)] << std::endl;
  //   std::cout << "wf(1,0,0)" << start[g.index(1,0,0)] << std::endl;
  //   std::cout << "wf(2,0,0)" << start[g.index(2,0,0)] << std::endl;
  //   std::cout << "wf(3,0,0)" << start[g.index(3,0,0)] << std::endl;
  //   std::cout << "norm(g) " << result << std::endl;
  //   exit(0);

  return result;
}


void wavefunction::solve(const wavefunction &a, const wavefunction &b, 
			 const wavefunction &c, const wavefunction &psi, long dimens)
{

  wavefunction e(dimens), f(dimens);
  long i;


  e[0]=-b[0]/a[0];
  f[0]=psi[0]/a[0];
  for (i=1; i<dimens; i++)
    {
      e[i]=(-c[i]/e[i-1] - b[i])/a[i];
      f[i]=(psi[i] + c[i]/e[i-1]*f[i-1])/a[i];
    };

  start[dimens-1]=-f[dimens-1]/e[dimens-1];
  for (i=dimens-2; i>=0; i--)
    {
      start[i]=(start[i+1]-f[i])/e[i];
    };
  

}

/*!

 */

void wavefunction::solve_du(const wavefunction &a, const wavefunction &b, 
			    const wavefunction &c, const wavefunction &psi, long dimens)
{

  wavefunction p(dimens-1), q(dimens-1);
  long i;


  p[0]=psi[0]/b[0];
  q[0]=a[0]/b[0];
  for (i=1; i<dimens-1; i++)
    {
      p[i]=(psi[i]-c[i]*p[i-1])/(b[i]-c[i]*q[i-1]);
      q[i]=a[i]/(b[i]-c[i]*q[i-1]);
    };

  start[dimens-1]=(psi[dimens-1]-c[dimens-1]*p[dimens-2])/(b[dimens-1]-c[dimens-1]*q[dimens-2]);

  for (i=dimens-2; i>=0; i--)
    {
      start[i]=p[i]-q[i]*start[i+1];
    };
  
}

void wavefunction::solve_toep(double a, double b, double b_upperleft, double b_lowerright, 
			      double c, const wavefunction &psi, long dimens)
{

  wavefunction f(dimens);
  fluid e(dimens);
  long i;
  double covera, bovera;

  covera=c/a;
  bovera=b/a;


  e[0]=-b_upperleft/a;
  f[0]=psi[0]/a;
  for (i=1; i<dimens-1; i++)
    {
      e[i]=-covera/e[i-1] - bovera;
      f[i]=psi[i]/a + covera/e[i-1]*f[i-1];
    };
  e[dimens-1]=-covera/e[dimens-2] - b_lowerright/a;
  f[dimens-1]=psi[dimens-1]/a + covera/e[dimens-2]*f[dimens-2];


  start[dimens-1]=-f[dimens-1]/e[dimens-1];
  for (i=dimens-2; i>=0; i--)
    {
      start[i]=(start[i+1]-f[i])/e[i];
    };
  

}


void wavefunction::apply(const wavefunction &a, const wavefunction &b, 
			 const wavefunction &c, const wavefunction &psi)
{

  long dimens=psi.wf_size();
  long i;


  start[0]=a[0]*psi[1]+b[0]*psi[0];

  for (i=1; i<dimens-1; i++)
    {
      start[i]=a[i]*psi[i+1]+b[i]*psi[i]+c[i]*psi[i-1];
    };
  
  start[dimens-1]=b[dimens-1]*psi[dimens-1]+c[dimens-1]*psi[dimens-2];


}

wavefunction wavefunction::orb_fieldenergies(double  time, grid g, hamop hamil, int me, const fluid &ms)
{
  wavefunction result(g.ngps_z());
  result.nullify();

  long xindex, lindex, zindex, index, index_lp, index_lm;

  for (zindex=0; zindex<g.ngps_z(); zindex++)
    {
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  lindex=0;
	  index=g.index(xindex,lindex,zindex);
	  index_lp=g.index(xindex,lindex+1,zindex);
	  result[zindex]+=g.delt_x()*g.r(xindex)*hamil.field(time,me)*conj(start[index])*sqrt((lindex+1.0)*(lindex+1.0)-ms[zindex]*ms[zindex])/sqrt((2.0*lindex+1.0)*(2.0*lindex+3.0))*start[index_lp];
	  
	  for (lindex=1; lindex<g.ngps_y()-1; lindex++)
	    {
	      index=g.index(xindex,lindex,zindex);
	      index_lp=g.index(xindex,lindex+1,zindex);
	      index_lm=g.index(xindex,lindex-1,zindex);
	      result[zindex]+=g.delt_x()*g.r(xindex)*hamil.field(time,me)*conj(start[index])*
		(sqrt((lindex+1.0)*(lindex+1.0)-ms[zindex]*ms[zindex])/sqrt((2.0*lindex+1.0)*(2.0*lindex+3.0))*start[index_lp]
		 +sqrt(lindex*lindex-ms[zindex]*ms[zindex])/sqrt((2.0*lindex-1.0)*(2.0*lindex+1.0))*start[index_lm]);
	    };

	  lindex=g.ngps_y()-1;
	  index_lm=g.index(xindex,lindex-1,zindex);
	  result[zindex]+=g.delt_x()*g.r(xindex)*hamil.field(time,me)*conj(start[index])*sqrt(lindex*lindex-ms[zindex]*ms[zindex])/sqrt((2.0*lindex-1.0)*(2.0*lindex+1.0))*start[index_lm];

	};
    };

  return result;

}


// single orbital case (lin. pol. only); overloaded below
complex wavefunction::expect_z(grid g)  // attention!!! works for s states only because m-dependence is suppressed
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index, index_ym, index_yp;
  double z,r;

  switch (g.dimens())
    {
    case 34 : 
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{ 
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      r=g.r(xindex);
	      yindex=0;
	      index=g.index(xindex,yindex,zindex);
	      index_yp=g.index(xindex,yindex+1,zindex);
	      
	      result=result+conj(start[index_yp])*(complex)r*start[index]*(complex)(yindex+1.0)/sqrt((2.0*yindex+1.0)*(2.0*yindex+3.0))*(complex)g.delt_x();

	      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  index_yp=g.index(xindex,yindex+1,zindex);
		  index_ym=g.index(xindex,yindex-1,zindex);
		  
		  result=result+(complex)g.delt_x()*(complex)r*start[index]
		    *(conj(start[index_yp])*(complex)(yindex+1.0)/sqrt((2.0*yindex+1.0)*(2.0*yindex+3.0))
		      +conj(start[index_ym])*(complex)yindex/sqrt((2.0*yindex-1.0)*(2.0*yindex+1.0)));
		  
		};		  
	      
	      yindex=g.ngps_y()-1;
	      index=g.index(xindex,yindex,zindex);
	      index_ym=g.index(xindex,yindex-1,zindex);
	      
	      result=result+(complex)(g.delt_x()*r)*start[index]*conj(start[index_ym])*(complex)yindex/sqrt((2.0*yindex-1.0)*(2.0*yindex+1.0));
	      
	    };
	};

      break;
    };

  return result;

}



// KS case; overloaded above
complex wavefunction::expect_z(grid g, fluid &degeneracies, const fluid &ms)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index, index_ym, index_yp;
  double z,r;

  switch (g.dimens())
    {
    case 34 : 
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{ 
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      if (g.ngps_y() > 1)  {
		r=g.r(xindex);
		yindex=0;
		index=g.index(xindex,yindex,zindex);
		index_yp=g.index(xindex,yindex+1,zindex);
		      
		result=result+conj(start[index_yp])*r*start[index]*sqrt((yindex+1.0)*(yindex+1)-ms[zindex]*ms[zindex])/sqrt((2.0*yindex+1.0)*(2.0*yindex+3.0))*g.delt_x()*degeneracies[zindex];
		      
		for (yindex=1; yindex<g.ngps_y()-1; yindex++)
		  {
		    index=g.index(xindex,yindex,zindex);
		    index_yp=g.index(xindex,yindex+1,zindex);
		    index_ym=g.index(xindex,yindex-1,zindex);
			  
		    result=result+g.delt_x()*degeneracies[zindex]*r*start[index]
		      *((conj(start[index_yp]))*(sqrt((yindex+1.0)*(yindex+1)-ms[zindex]*ms[zindex]))/sqrt((2.0*yindex+1.0)*(2.0*yindex+3.0))
			+(conj(start[index_ym]))*(sqrt(yindex*yindex-ms[zindex]*ms[zindex]))/sqrt((2.0*yindex-1.0)*(2.0*yindex+1.0)));
			  
		  };		  
		      
		yindex=g.ngps_y()-1;
		index=g.index(xindex,yindex,zindex);
		index_ym=g.index(xindex,yindex-1,zindex);
		    
		result=result+(g.delt_x()*degeneracies[zindex]*r)*start[index]*conj(start[index_ym])*(sqrt(yindex*yindex-ms[zindex]*ms[zindex]))/sqrt((2.0*yindex-1.0)*(2.0*yindex+1.0));
	      };    
	    };
	};
	  
      break;
    default :
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  z=g.z(zindex);
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		      
		  result=result+z*
		    real(conj(start[index])*start[index])*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    };
  
  return result;

}



complex wavefunction::expect_cycl_pol_plus(grid g)
{

  complex result(0.0,0.0);
  long xindex, ell, m;
  long index, index_lmmm, index_lpmm;
  double r;

  switch (g.dimens())
    {
    case 44 :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  r=g.r(xindex);
	  ell=0; m=0;
	  index=g.rlmindex(xindex,ell,m);
	  index_lpmm=g.rlmindex(xindex,ell+1,m-1);
	  result=result+r*conj(start[index])*(-start[index_lpmm]
					      *sqrt((ell+1.0)*(ell-m+1.0)*(ell-m+2.0)
						    /((2.0*ell+1)*(2.0*ell+2)*(2.0*ell+3))));
	  for (ell=1; ell<g.ngps_y()-1; ell++)
	    {
	      for (m=-ell; m<=ell; m++)
		{
		  index=g.rlmindex(xindex,ell,m);
		  index_lmmm=g.rlmindex(xindex,ell-1,m-1);
		  index_lpmm=g.rlmindex(xindex,ell+1,m-1);
		  result=result+r*conj(start[index])*(start[index_lmmm]
						      *sqrt((ell+m-1.0)*(ell+m)/(2.0*(2.0*ell-1)*(2.0*ell+1)))
						      -start[index_lpmm]
						      *sqrt((ell+1.0)*(ell-m+1.0)*(ell-m+2.0)
							    /((2.0*ell+1)*(2.0*ell+2)*(2.0*ell+3))));
		};
	    };
	  ell=g.ngps_y()-1;
	  for (m=-ell; m<=ell; m++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_lmmm=g.rlmindex(xindex,ell-1,m-1);
	      index_lpmm=g.rlmindex(xindex,ell+1,m-1);
	      result=result+r*conj(start[index])*(start[index_lmmm]
						  *sqrt((ell+m-1.0)*(ell+m)/(2.0*(2.0*ell-1)*(2.0*ell+1))));

	    };
	};
      result=-SQRTT*g.delt_x()*result;
      break;
    };
  
  return result;

}



complex wavefunction::expect_cycl_pol_minus(grid g)
{

  complex result(0.0,0.0);
  long xindex, ell, m;
  long index, index_lmmp, index_lpmp;
  double r;

  switch (g.dimens())
    {
    case 44 :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  r=g.r(xindex);
	  ell=0; m=0;
	  index=g.rlmindex(xindex,ell,m);
	  index_lpmp=g.rlmindex(xindex,ell+1,m+1);
	  result=result+r*conj(start[index])*(-start[index_lpmp]
					      *sqrt((ell+1.0)*(ell+m+1.0)*(ell+m+2.0)
						    /((2.0*ell+1)*(2.0*ell+2)*(2.0*ell+3))));
	  for (ell=1; ell<g.ngps_y()-1; ell++)
	    {
	      for (m=-ell; m<=ell; m++)
		{
		  index=g.rlmindex(xindex,ell,m);
		  index_lmmp=g.rlmindex(xindex,ell-1,m+1);
		  index_lpmp=g.rlmindex(xindex,ell+1,m+1);
		  result=result+r*conj(start[index])*(start[index_lmmp]
						      *sqrt((ell-m-1.0)*(ell-m)/(2.0*(2.0*ell-1)*(2.0*ell+1)))
						      -start[index_lpmp]
						      *sqrt((ell+1.0)*(ell+m+1.0)*(ell+m+2.0)
							    /((2.0*ell+1)*(2.0*ell+2)*(2.0*ell+3))));
		};
	    };
	  ell=g.ngps_y()-1;
	  for (m=-ell; m<=ell; m++)
	    {
	      index=g.rlmindex(xindex,ell,m);
	      index_lmmp=g.rlmindex(xindex,ell-1,m+1);
	      index_lpmp=g.rlmindex(xindex,ell+1,m+1);
	      result=result+r*conj(start[index])*(start[index_lmmp]
						  *sqrt((ell-m-1.0)*(ell-m)/(2.0*(2.0*ell-1)*(2.0*ell+1))));

	    };
	};
      result=SQRTT*g.delt_x()*result;
      break;
    };
  
  return result;

}





complex wavefunction::expect_rr(grid g, fluid &degeneracies)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;
  double r;


  for (zindex=0; zindex<g.ngps_z(); zindex++)
    { 
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  r=g.r(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      result=result+g.delt_x()*degeneracies[zindex]
		*start[index]*conj(start[index])*r*r;
	    };		  
	};
    };

  return result;

}


void wavefunction::dump_expect_z(grid g, FILE* os, fluid &degeneracies, const fluid &ms)
{
  
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index, index_ym, index_yp;
  double z,r;

  switch (g.dimens())
    {
    case 34 : 
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{ 
	  result=complex(0.0,0.0);
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      r=g.r(xindex);
	      yindex=0;
	      index=g.index(xindex,yindex,zindex);
	      index_yp=g.index(xindex,yindex+1,zindex);

	      result=result+conj(start[index_yp])*r*start[index]*sqrt((yindex+1.0)*(yindex+1)-ms[zindex]*ms[zindex])/sqrt((2.0*yindex+1.0)*(2.0*yindex+3.0))*g.delt_x()*degeneracies[zindex];

	      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  index_yp=g.index(xindex,yindex+1,zindex);
		  index_ym=g.index(xindex,yindex-1,zindex);

		  result=result+g.delt_x()*degeneracies[zindex]*r*start[index]
		    *((conj(start[index_yp]))*(sqrt((yindex+1.0)*(yindex+1)-ms[zindex]*ms[zindex]))/sqrt((2.0*yindex+1.0)*(2.0*yindex+3.0))
		      +(conj(start[index_ym]))*(sqrt(yindex*yindex-ms[zindex]*ms[zindex]))/sqrt((2.0*yindex-1.0)*(2.0*yindex+1.0)));

		}		  

	      yindex=g.ngps_y()-1;
	      index=g.index(xindex,yindex,zindex);
	      index_ym=g.index(xindex,yindex-1,zindex);

	      result=result+(g.delt_x()*
			     degeneracies[zindex]*r)*start[index]*
		conj(start[index_ym])*
		(sqrt(yindex*yindex-ms[zindex]*ms[zindex]))/
		sqrt((2.0*yindex-1.0)*(2.0*yindex+1.0));

	    }
	  fprintf(os,"%e %e ",real(result),imag(result));
	}
      fprintf(os,"\n");
      break;
    }
}
