#include<fluid.h>
#include<hamop.h>
//#include<propmanage.h>
#include<wavefunction.h>
#include<grid.h>

#define pi 3.1415926
#define complex std::complex<double>

class wavefunction;

void fluid::nullify()
{
  for (long i = 0; i < fluid_dim; i++)
    start[i] = 0.0;

}


// *** now some operators follow
fluid &fluid::operator=(const fluid &v)
{
    if (this != &v)
    {
       delete[] start;
       fluid_dim  = v.fluid_dim;
       start = new double[fluid_dim];
       for (long i = 0; i < fluid_dim; i++)
           start[i] = v.start[i];
    }
    return *this;
}

fluid &fluid::operator=(const wavefunction &v)
{
  delete[] start;
  fluid_dim  = v.wf_size();
  start = new double[wf_size()];
  for (long i = 0; i < fluid_dim; i++)
    start[i] = real(v[i]);
  
  return *this;
}

fluid& fluid::operator *= (double z)
{
  for(long i=0; i<fluid_dim; i++)
    start[i]=start[i]*z;
  return *this;
}

fluid operator * (
			 double z, 
			 const fluid &v
			 )
{
  fluid temp=v;
  return temp *= z;
}

fluid operator * (
			 const fluid &v, 
			 double z
			 )
{
  fluid temp=v;
  return temp *= z;
}

double operator * (
			    const fluid &v, 
			    const fluid &w 
			    )
{
  double result=0.0;
  for(long i=0; i<v.wf_size(); i++)
    {
      result+=v[i]*w[i];
    };
  return result;
}

fluid operator + (
			    const fluid &v, 
			    const fluid &w 
			    )
{
  fluid result=v;
  for(long i=0; i<v.wf_size(); i++)
    {
      result[i]=v[i]+w[i];
    };
  return result;
}

fluid operator - (
			    const fluid &v, 
			    const fluid &w 
			    )
{
  fluid result=v;
  for(long i=0; i<v.wf_size(); i++)
    {
      result[i]=v[i]-w[i];
    };
  return result;
}


std::ostream& operator<<(
			 std::ostream& os, 
			 const fluid& v)
{  
  for(long i = 0; i < v.wf_size(); i++)
    {   os << v[i] << std::endl;
    }
  return os;
}


//  istream& operator>>(
//  		    istream& is, 
//  		    fluid& v
//  		    )
//  {  
//    double tmpre;
//    for(long i = 0; i < v.wf_size(); i++)
//      {   
//        is >> tmpre;
//        v[i]=tmpre;
//      }
//    return is;
//  }

//  istream& operator>>(
//  		    ifstream& is, 
//  		    fluid& v
//  		    )
//  {  
//    double tmpre;
//    for(long i = 0; i < v.wf_size(); i++)
//      {   
//        is >> tmpre;
//        v[i]=tmpre;
//      }
//    return is;
//  }

void fluid::dump_to_file_sh(grid g, FILE* os, int stepwidth)
{
  for (long i=0; i<g.ngps_x(); i=i+stepwidth) fprintf(os,"%16.12e\n",start[i]);
}


void fluid::load(FILE* os, int iv)
{
  double value;
  long ctrl_id;

  for (long i=0; i<wf_size(); i++)
  {
    ctrl_id=fscanf(os,"%lf",&value); start[i]=value;
  }
}




void fluid::embed_as_x(grid g, long yindex, long zindex, const fluid &v)
{
  long i,index;
  for (i=0; i<v.wf_size(); i++)
  {
    index=g.index(i,yindex,zindex);
    start[index]=v[i];
  }
}

fluid fluid::extract_x(grid g, long other_one, long other_two)
{

  long i, index;

  fluid result(g.ngps_x());
  for (i=0; i<g.ngps_x(); i++)
  {
    index=g.index(i,other_one,other_two);
    result[i]=start[index];
  }

  return result;
}


int fluid::init(long isize)
{
   fluid_dim = isize;
   start = new double[isize];
   for (long i = 0; i < isize; i++) start[i] = 0.0;
   return 0;
}


double fluid::integrate_along_x(grid g)
{
  long i;
  double result=0.0;

  for (i=0; i<g.ngps_x(); i++)
    {
      result+=start[i]*g.delt_x();
    };

  return result;
}

double fluid::expect_rr(grid g)
{
  long i;
  double result=0.0;
  double r;

  for (i=0; i<g.ngps_x(); i++)
    {
      r=g.r(i);
      result+=start[i]*g.delt_x()*r*r;
    };

  return result;
}

double fluid::expect_oneoverr(grid g)
{
  long i;
  double result=0.0;
  double r;

  for (i=0; i<g.ngps_x(); i++)
    {
      r=g.r(i);
      result+=start[i]*g.delt_x()/r;
    };

  return result;
}


double fluid::calculate_ionkinenergy(const fluid &Q, const double singlecharge, const double singlemass)
{
  double result=0.0;
  int i;

  for (i=0; i<Q.wf_size(); i++)
    {
      result+=Q[i]/singlecharge*singlemass*start[2*i+1]*start[2*i+1];
    };
  
  return 0.5*result;

}

double fluid::calculate_ionpotenergy(const fluid &Q)
{
  double result=0.0;
  int i,j;

  for (i=0; i<Q.wf_size(); i++)
    {
      for (j=i+1; j<Q.wf_size(); j++)
	{
	  result+=Q[i]*Q[j]/start[2*j];
	};
    };
  
  return result;

}
