#include<cmatrix.h>
#include<wavefunction.h>
#include<fluid.h>
#include<grid.h>
#include<hamop.h>

#define complex std::complex<double>

cmatrix &cmatrix::operator=(const cmatrix &v)
{
   if (this != &v)
   {
      delete[] start;
      cm__dim  = v.cm__dim;
      cm__rows  = v.cm__rows;
      cm__cols  = v.cm__cols;
      start = new complex[cm__dim];
      for (long i = 0; i < cm__dim; i++)
	 start[i] = v.start[i];
   }
   return *this;
}


complex cmatrix::elem(long row, long col)
{
   long index;

   index=row*cm__cols+col;

   return start[index];
}

complex cmatrix::elem(long row, long col) const
{
   long index;

   index=row*cm__cols+col;

   return start[index];
}


std::ostream& operator<<(std::ostream& os, const cmatrix& v)
{  
   for(long i = 0; i < v.cm_rows(); i++)
   {
      for(long j = 0; j < v.cm_cols(); j++)
      {
	 if (fabs(real(v.elem(i,j)))<1e-10)
	 {
	    os << "(0,";
	 }
	 else
	 {
	    os << "(" << real(v.elem(i,j)) << ",";
	 };
	 if (fabs(imag(v.elem(i,j)))<1e-10)
	 {
	    os << "0)" << " ";
	 }
	 else
	 {
	    os << imag(v.elem(i,j)) << ") ";
	 };

      };
      os << std::endl;
   };
   os << std::endl;
  
   return os;
}

void cmatrix::set_elem(long row, long col, complex value)
{
   long index;
   index=row*cm__cols+col;
   start[index]=value;
}

void cmatrix::set_elem(long row, long col, double value)
{
   long index;
   index=row*cm__cols+col;
   start[index]=value;
}

void cmatrix::set_unity()
{
   long index,i,j;
   long smaller;
   if (cm__rows>cm__cols) {smaller=cm__cols; } else { smaller=cm__rows;};

   for(i = 0; i < cm__rows; i++)
   {
      for(j = 0; j < cm__cols; j++)
      {  
	 index=i*cm__cols+j;
	 start[index]=0.0;
      };
   };

   for (i=0; i<smaller; i++)
   {
      j=i;
      index=i*cm__cols+j;
      start[index]=1.0;
   }; 

    

}

void cmatrix::set_zero()
{
   long index,i,j;
   long smaller;
   if (cm__rows>cm__cols) {smaller=cm__cols; } else { smaller=cm__rows;};

   for(i = 0; i < cm__rows; i++)
   {
      for(j = 0; j < cm__cols; j++)
      {  
	 index=i*cm__cols+j;
	 start[index]=0.0;
      };
   };

}

cmatrix operator * (complex a, const cmatrix &b)
{
   long i,j;
   cmatrix result(b.cm_rows(),b.cm_cols());
   for (i=0; i<b.cm_rows(); i++)
   {
      for (j=0; j<b.cm_cols(); j++)
      {
	 result.set_elem(i,j,b.elem(i,j)*a);
      };
   };
   return result;
}

cmatrix operator * (const cmatrix &b, complex a)
{
   long i,j;
   cmatrix result(b.cm_rows(),b.cm_cols());
   for (i=0; i<b.cm_rows(); i++)
   {
      for (j=0; j<b.cm_cols(); j++)
      {
	 result.set_elem(i,j,b.elem(i,j)*a);
      };
   };
   return result;
}

cmatrix operator * (double a, const cmatrix &b)
{
   long i,j;
   cmatrix result(b.cm_rows(),b.cm_cols());
   for (i=0; i<b.cm_rows(); i++)
   {
      for (j=0; j<b.cm_cols(); j++)
      {
	 result.set_elem(i,j,b.elem(i,j)*a);
      };
   };
   return result;
}

cmatrix operator * (const cmatrix &b, double a)
{
   long i,j;
   cmatrix result(b.cm_rows(),b.cm_cols());
   for (i=0; i<b.cm_rows(); i++)
   {
      for (j=0; j<b.cm_cols(); j++)
      {
	 result.set_elem(i,j,b.elem(i,j)*a);
      };
   };
   return result;
}


cmatrix operator * (const cmatrix &a, const cmatrix &b)
{
   long i,j,k;
   cmatrix result(a.cm_rows(),b.cm_cols());
   complex entry;

   for (i=0; i<a.cm_rows(); i++)
   {
      for (j=0; j<b.cm_cols(); j++)
      {
	 entry=complex(0.0,0.0);
	 for (k=0; k<a.cm_cols(); k++)
	 {
	    entry+=a.elem(i,k)*b.elem(k,j);
	 };
	 result.set_elem(i,j,entry);
      };
   };
   return result;

}

cmatrix operator + (const cmatrix &a, const cmatrix &b)
{
   long i,j,smaller_rows,smaller_cols;

   if (a.cm_cols()<b.cm_cols()) { smaller_cols=a.cm_cols(); } else { smaller_cols=b.cm_cols(); };
   if (a.cm_rows()<b.cm_rows()) { smaller_rows=a.cm_rows(); } else { smaller_rows=b.cm_rows(); };

   cmatrix result(smaller_rows,smaller_cols);

   for (i=0; i<smaller_rows; i++)
   {
      for (j=0; j<smaller_cols; j++)
      {
	 result.set_elem(i,j,a.elem(i,j)+b.elem(i,j));
      };
   };
   return result;

}

cmatrix operator - (const cmatrix &a, const cmatrix &b)
{
   long i,j,smaller_rows,smaller_cols;

   if (a.cm_cols()<b.cm_cols()) { smaller_cols=a.cm_cols(); } else { smaller_cols=b.cm_cols(); };
   if (a.cm_rows()<b.cm_rows()) { smaller_rows=a.cm_rows(); } else { smaller_rows=b.cm_rows(); };

   cmatrix result(smaller_rows,smaller_cols);

   for (i=0; i<smaller_rows; i++)
   {
      for (j=0; j<smaller_cols; j++)
      {
	 result.set_elem(i,j,a.elem(i,j)-b.elem(i,j));
      };
   };
   return result;

}

cmatrix cmatrix::mult_cm_with_diag(const wavefunction &wf) 
{
   long i,j,k;

   cmatrix result(cm__rows,wf.wf_size());

   for (i=0; i<cm__rows; i++)
   {
      for (j=0; j<wf.wf_size(); j++)
      {
	 result.set_elem(i,j,elem(i,j)*wf[j]);
      };
   };
   return result;

}

cmatrix cmatrix::mult_diag_with_cm(const wavefunction &wf)  
{
   long i,j,k;

   cmatrix result(wf.wf_size(),cm__cols);

   for (i=0; i<wf.wf_size(); i++)
   {
      for (j=0; j<cm__cols; j++)
      {
	 result.set_elem(i,j,elem(i,j)*wf[i]);
      };
   };
   return result;

}

void cmatrix::embedding_diag(const wavefunction &wf)  
{
   long i,j;
   long smaller;
   if (cm__rows>cm__cols) {smaller=cm__cols; } else { smaller=cm__rows;};

   for (i=0; i<smaller; i++)
   {
      for (j=0; j<smaller; j++)
      {
	 set_elem(i,j,0.0);
      };
   };
   for (i=0; i<smaller; i++)
   {
      set_elem(i,i,wf[i]);
   };


}


void cmatrix::ini_tridiag(complex a, complex b, complex c)
{
   long i,j;
   long smaller;
  
   set_zero();
  
   if (cm__rows>cm__cols) {smaller=cm__cols; } else { smaller=cm__rows;};

   set_elem(0,0,b);
   set_elem(0,1,c);

   for (i=1; i<smaller-1; i++)
   {
      set_elem(i,i-1,a);
      set_elem(i,i,b);
      set_elem(i,i+1,c);
   };
  
   set_elem(smaller-1,smaller-1,b);
   set_elem(smaller-1,smaller-2,a);

}

cmatrix operator + (const cmatrix &a, const wavefunction &wf)
{
   long i,j,k;
   long smaller;
   if (a.cm_rows()>a.cm_cols()) {smaller=a.cm_cols(); } else { smaller=a.cm_rows();};

   cmatrix result=a;

   for (i=0; i<smaller; i++)
   {
      result.set_elem(i,i,a.elem(i,i)+wf[i]);
   };
   return result;

}

cmatrix operator + (const wavefunction &wf, const cmatrix &a)
{
   long i,j,k;
   long smaller;
   if (a.cm_rows()>a.cm_cols()) {smaller=a.cm_cols(); } else { smaller=a.cm_rows();};

   cmatrix result=a;

   for (i=0; i<smaller; i++)
   {
      result.set_elem(i,i,a.elem(i,i)+wf[i]);
   };
   return result;

}
cmatrix operator + (const cmatrix &a, complex b)
{
   long i,j,k;
   long smaller;
   if (a.cm_rows()>a.cm_cols()) {smaller=a.cm_cols(); } else { smaller=a.cm_rows();};

   cmatrix result=a;

   for (i=0; i<smaller; i++)
   {
      result.set_elem(i,i,a.elem(i,i)+b);
   };
   return result;

}

cmatrix operator + (complex b, const cmatrix &a)
{
   long i,j,k;
   long smaller;
   if (a.cm_rows()>a.cm_cols()) {smaller=a.cm_cols(); } else { smaller=a.cm_rows();};

   cmatrix result=a;

   for (i=0; i<smaller; i++)
   {
      result.set_elem(i,i,a.elem(i,i)+b);
   };
   return result;

}




cmatrix cmatrix::transp()
{
   cmatrix result(cm__cols,cm__rows);
   long i,j;

   for (i=0; i<cm__cols; i++)
   {
      for (j=0; j<cm__rows; j++)
      { 
	 result.set_elem(i,j,elem(j,i));
      };
   };
	
   return result;  
}

cmatrix cmatrix::adjoint()
{
   cmatrix result(cm__cols,cm__rows);
   long i,j;

   for (i=0; i<cm__cols; i++)
   {
      for (j=0; j<cm__rows; j++)
      { 
	 result.set_elem(i,j,conj(elem(j,i)));
      };
   };
	
   return result;  
}

// Cholesky factorization + solving
void cmatrix::lap_dsytrf_and_dsytrs(fluid &rhs, int info)
{
   char uplo='U';
   long n=cm_rows();
   double a[cm_size()];
   long lda=cm_rows();  
   double b[cm_rows()];
   long nrhs=1;  
   long ldb=cm_rows();
   long ipiv[n];
   long lwork=4*n;
   double work[lwork];

   long i,j;

   for (i=0; i<cm_rows(); i++)
   {
      b[i]=rhs[i];
      //      cout << b[i] << endl;
   };

   for (i=0; i<cm_rows(); i++)
   {
      for (j=0; j<cm_cols(); j++)
      {
	 a[i*cm_cols()+j]=real(elem(j,i));
      };
   };


   dsytrf_( &uplo, &n, a, &lda, ipiv, work, &lwork, &info);
   //  cout << "info's:  " << info;
   dsytrs_( &uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
   //  cout << ", " << info << endl;

   for (i=0; i<cm_rows(); i++)
   {
      rhs[i]=b[i];
   };

}


void cmatrix::lap_dgetrf_and_dgetrs(fluid &rhs, int info)
{
   char trans='N';
   long m=cm_cols();
   long n=cm_rows();
   double a[cm_size()];
   long lda=cm_rows();  
   double b[cm_rows()];
   long nrhs=1;  
   long ldb=cm_rows();
   long ipiv[n];

   long i,j;

   for (i=0; i<cm_rows(); i++)
   {
      b[i]=rhs[i];
      //      cout << b[i] << endl;
   };

   for (i=0; i<cm_rows(); i++)
   {
      for (j=0; j<cm_cols(); j++)
      {
	 a[i*cm_cols()+j]=real(elem(j,i));
      };
   };


   dgetrf_( &m, &n, a, &lda, ipiv, &info);
   //  cout << "info " << info << endl;
   dgetrs_( &trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
   //  cout << "info " << info << endl;

   for (i=0; i<cm_rows(); i++)
   {
      rhs[i]=b[i];
   };

}


void cmatrix::lap_dstegr(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex)
{
   long i,j,index;
   long size=lastindex-firstindex+1;

   double d[size], e[size];

   for (i=0; i<size; i++)          
   {                               
      d[i]=real(elem(firstindex+i,firstindex+i));
   };             
   for (i=0; i<size-1; i++)          
   {                               
      e[i]=real(elem(firstindex+i+1,firstindex+i));
   };             

   char jobz='V';
   char range='A';
   int n=size;
   double vu=0.0;
   double vl=0.0;
   int iu=0;
   int il=0;
   double abstol=0.0;
   int m;
   double w[size];
   double z[size*size];
   int ldz=size;
   int isuppz[2*size];
   double work[18*size];
   int lwork=18*size;
   int iwork[10*size];
   int liwork=10*size;
   int info;

   dstegr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info );

   for (i=0; i<size; i++)          
   {                               
      eval[firstindex+i]=complex(w[i],0.0);
   };


   for (i=0; i<size; i++)
   {
      for (j=0; j<size; j++)
      {
	 evec.set_elem(firstindex+j,firstindex+i,complex(z[j+size*i],0.0));
      };
   };

}


void cmatrix::lap_dsyev(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex)
{
   long i,j,index;
   long size=lastindex-firstindex+1;


   char jobz='V';
   char uplo='U';
   long n=size;
   double a[size*size];

   for (i=0; i<size; i++)
   {
      for (j=0; j<size; j++)
      {
	 a[j+size*i]=real(elem(firstindex+j,firstindex+i));
      };
   };

   long lda=size;
   double w[size];
   long lwork=3*size-1;
   double work[lwork];
   int info;

   dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info );

   for (i=0; i<size; i++)          
   {                               
      eval[firstindex+i]=complex(w[i],0.0);
   };


   for (i=0; i<size; i++)
   {
      for (j=0; j<size; j++)
      {
	 evec.set_elem(firstindex+j,firstindex+i,complex(a[j+size*i],0.0));
      };
   };

}



void cmatrix::lap_zgeev(wavefunction &eval, cmatrix &levec, cmatrix &revec, long firstindex, long lastindex)
{
   long i,j,index;
   long size=lastindex-firstindex+1;
   int lda, ldvl, ldvr,lwork, n, ok;

   double AT[2*size*size],vl[2*size*size],vr[2*size*size];   
   double w[2*size], WORK[4*size], RWORK[2*size];

   char jobvl,jobvr;
   jobvl='V';
   jobvr='V';

   for (i=0; i<size; i++)          
   {                               
      for(j=0; j<size; j++) 
      {
	 AT[2*(j+size*i)]=real(elem(firstindex+j,firstindex+i));
	 AT[2*(j+size*i)+1]=imag(elem(firstindex+j,firstindex+i));
      }             
   }

   n=size;       
   lda=size;
   ldvr=size;    
   ldvl=size; 
   lwork=2*size; 

   zgeev_(&jobvl, &jobvr, &n, AT, &lda, w, vl, &ldvl, vr, &ldvr, WORK, &lwork, RWORK, &ok);

   for (i=0; i<size; i++)
   {
      eval[firstindex+i]=complex(w[2*i],w[2*i+1]);
      for (j=0; j<size; j++)
      {
	 revec.set_elem(firstindex+j,firstindex+i,complex(vr[2*(j+size*i)],vr[2*(j+size*i)+1]));
	 levec.set_elem(firstindex+j,firstindex+i,complex(vl[2*(j+size*i)],vl[2*(j+size*i)+1]));
      };
   };

}

void cmatrix::lap_zheev(wavefunction &eval, cmatrix &evec, long firstindex, long lastindex)
{
   long i,j,index;
   long size=lastindex-firstindex+1;

   char jobz='V';
   char uplo='U';
   long n=size;
   double a[2*size*size];
   long lda=size;
   double w[size];
   long lwork=2*size-1;
   double work[2*lwork];
   double rwork[3*size-2];
   int info;
  
   for (i=0; i<size; i++)          
   {                               
      for(j=0; j<size; j++) 
      {
	 a[2*(j+size*i)]=real(elem(firstindex+j,firstindex+i));
	 a[2*(j+size*i)+1]=imag(elem(firstindex+j,firstindex+i));
      }             
   }

   zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);

   for (i=0; i<size; i++)
   {
      eval[firstindex+i]=complex(w[i],0.0);
      for (j=0; j<size; j++)
      {
	 evec.set_elem(firstindex+j,firstindex+i,complex(a[2*(j+size*i)],a[2*(j+size*i)+1]));
      };
   };

}

wavefunction cmatrix::expand_in_colvecs(const wavefunction  &c)
{
   long i;
   wavefunction result(cm__cols);
   result.nullify();


   for (i=0; i<cm__cols; i++)
   {
      result=result + c[i]*cm_pick_col(i);
   };

   return result;

}


wavefunction cmatrix::cm_pick_col(long j)
{
   long i;
   wavefunction result(cm__rows);

   for (i=0; i<cm__rows; i++)
   {
      result[i]=elem(i,j);
   };
   return result;
}

wavefunction cmatrix::cm_pick_row(long j)
{
   long i;
   wavefunction result(cm__cols);

   for (i=0; i<cm__cols; i++)
   {
      result[i]=elem(j,i);
   };
   return result;
}

wavefunction cmatrix::cm_pick_diag()
{
   long i;
   long smaller;
   if (cm__rows>cm__cols) {smaller=cm__cols; } else { smaller=cm__rows;};  
   wavefunction result(smaller);

   for (i=0; i<smaller; i++)
   {
      result[i]=elem(i,i);
   };
   return result;
}

void cmatrix::build_z_in_sturmian_basis(grid g)
{
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long ll,nn,l,n,index,iindex;
   complex fac;
   complex imagi(0.0,1.0);


   set_zero();  

   // z representation in Sturmian basis
   for (ll=0; ll<=max_l-1; ll++)
   {
      for (nn=1+ll; nn<=max_n; nn++)
      { 
	 index=g.nlindex(nn-1,ll);
	 l=ll+1;
	 for (n=1+l; n<=max_n; n++)
	 { 
	    iindex=g.nlindex(n-1,l);
	    fac=(ll+1.0)/sqrt((2.0*ll+1.0)*(2.0*ll+3.0));
	    if (n==nn)
	    {
	       set_elem(index,iindex,-fac*3.0*sqrt((nn-ll-1.0)*(nn+ll+1.0)));
	    };
	    if (n==nn+1)
	    {
	       set_elem(index,iindex,fac*(2.0*nn-ll)*sqrt(((nn+ll+2.0)*(nn+ll+1.0))/(nn*(nn+1.0))));
	    };
	    if (n==nn-1)
	    {
	       set_elem(index,iindex,fac*(2.0*nn+ll)*sqrt(((nn-ll-2.0)*(nn-ll-1.0))/(nn*(nn-1.0))));
	    };
	    if (n==nn+2)
	    {
	       set_elem(index,iindex,-0.5*fac*sqrt((nn-ll)*(nn+ll+3.0)*(nn+ll+2.0)*(nn+ll+1.0)/(nn*(nn+2.0))));
	    };
	    if (n==nn-2)
	    {
	       set_elem(index,iindex,-0.5*fac*sqrt((nn+ll)*(nn-ll-3.0)*(nn-ll-2.0)*(nn-ll-1.0)/(nn*(nn-2.0))));
	    };
	 };
      };
   };
   for (ll=1; ll<=max_l; ll++)
   {
      for (nn=1+ll; nn<=max_n; nn++)
      { 
	 index=g.nlindex(nn-1,ll);
	 l=ll-1;
	 for (n=1+l; n<=max_n; n++)
	 { 
	    iindex=g.nlindex(n-1,l);
	    fac=1.0*ll/sqrt((2.0*ll-1.0)*(2.0*ll+1.0));
	    if (n==nn)
	    {
	       set_elem(index,iindex,-fac*3.0*sqrt((nn-ll)*(nn+ll)));
	    };
	    if (n==nn+1)
	    {
	       set_elem(index,iindex,fac*(2.0*nn+ll+1.0)*sqrt(((nn-ll)*(nn-ll+1.0))/(nn*(nn+1.0))));
	    };
	    if (n==nn-1)
	    {
	       set_elem(index,iindex,fac*(2.0*nn-ll-1.0)*sqrt(((nn+ll)*(nn+ll-1.0))/(nn*(nn-1.0))));
	    };
	    if (n==nn+2)
	    {
	       set_elem(index,iindex,-0.5*fac*sqrt((nn-ll)*(nn-ll+1.0)*(nn-ll+2.0)*(nn+ll+1.0)/(nn*(nn+2.0))));
	    };
	    if (n==nn-2)
	    {
	       set_elem(index,iindex,-0.5*fac*sqrt((nn+ll)*(nn+ll-1.0)*(nn+ll-2.0)*(nn-ll-1.0)/(nn*(nn-2.0))));
	    };
	 };
      };
   };


};

void cmatrix::build_pz_in_sturmian_basis(grid g)
{
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long ll,nn,l,n,index,iindex;
   complex fac;
   complex imagi(0.0,1.0);


   set_zero();

   // p_z=-i d/dz in Sturmian representation
   for (ll=0; ll<=max_l-1; ll++)
   {
      for (nn=1+ll; nn<=max_n; nn++)
      { 
	 index=g.nlindex(nn-1,ll);
	 l=ll+1;
	 for (n=1+l; n<=max_n; n++)
	 { 
	    iindex=g.nlindex(n-1,l);
	    fac=-imagi*(ll+1.0)/sqrt((2.0*ll+1.0)*(2.0*ll+3.0));
	    if (n==nn+1)
	    {
	       set_elem(index,iindex,fac*sqrt(((nn+ll+2.0)*(nn+ll+1.0))/(nn*(nn+1.0))));
	    };
	    if (n==nn-1)
	    {
	       set_elem(index,iindex,-fac*sqrt(((nn-ll-2.0)*(nn-ll-1.0))/(nn*(nn-1.0))));
	    };
	 };
      };
   };
   for (ll=1; ll<=max_l; ll++)
   {
      for (nn=1+ll; nn<=max_n; nn++)
      { 
	 index=g.nlindex(nn-1,ll);
	 l=ll-1;
	 for (n=1+l; n<=max_n; n++)
	 { 
	    iindex=g.nlindex(n-1,l);
	    fac=-imagi*(ll*1.0)/sqrt((2.0*ll-1.0)*(2.0*ll+1.0));
	    if (n==nn+1)
	    {
	       set_elem(index,iindex,fac*sqrt(((nn-ll)*(nn-ll+1.0))/(nn*(nn+1.0))));
	    };
	    if (n==nn-1)
	    {
	       set_elem(index,iindex,-fac*sqrt(((nn+ll)*(nn+ll-1.0))/(nn*(nn-1.0))));
	    };
	 };
      };
   };

};
      

void cmatrix::build_sturmian_overlap_matrix(grid g)
{
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long l,n,index;

   set_unity();

   for (l=0; l<=max_l; l++)
   {
      for (n=1+l; n<=max_n; n++)
      { 
	 index=g.nlindex(n-1,l);
	 set_elem(index,index,1.0);
	 if (n>1+l)
	 {
	    set_elem(index,index-1,-0.5*sqrt((n-1-l)*(n+l)/((n-1.0)*n)));
	 };	
	 if (n<max_n)
	 {
	    set_elem(index,index+1,-0.5*sqrt((n-l)*(n+l+1)/((n+1.0)*n)));
	 };
      };
   }; 
} 

cmatrix  cmatrix::build_trafo_sturmian_to_orthogonal(grid g,  wavefunction &evals)
{
   long size=g.nlindex(g.ngps_x()-1,g.ngps_y()-1)+1;
   cmatrix cmW(size,size);
   cmatrix res(size,size);
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long l;

   for (l=0; l<=max_l; l++)
   {
      lap_dstegr(evals,cmW,g.nlindex(l+1-1,l),g.nlindex(max_n-1,l));
   };

   res=cmW.mult_cm_with_diag((evals.invert()).sqrtrealwf());

   return res;

}


void cmatrix::build_rdr_matrix_in_sturmian_basis(
   grid g,
   cmatrix &cmS,
   long size
   )
{
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long l,n,index;
   complex imagi(0.0,1.0);

   cmatrix cmB(size,size);

   for (l=0; l<=max_l; l++)
   {
      for (n=1+l; n<=max_n; n++)
      { 
	 index=g.nlindex(n-1,l);
	 cmB.set_elem(index,index,-0.5);
	 if (n>1+l)
	 {
	    cmB.set_elem(index,index-1,0.5*sqrt((n-1-l)*(n+l)*n/((n-1.0))));
	 };	
	 if (n>2+l)
	 {
	    cmB.set_elem(index,index-2,-0.25*sqrt((n+l)*(n-1.0+l)*(n-2.0-l)*(n-1.0-l)/(n*(n-2.0))));
	 };
	 if (n<max_n)
	 {
	    cmB.set_elem(index,index+1,-0.5*sqrt((n-l)*(n+l+1)*n/((n+1.0))));
	 };
	 if (n<max_n-1)
	 {
	    cmB.set_elem(index,index+2,0.25*sqrt((n-l+1.0)*(n-l)*(n+l+1.0)*(n+l+2.0)/(n*(n+2.0))));
	 };
      };
   };

   (*this)=cmB-cmS;


}

cmatrix cmatrix::build_orth_from_symmetric_matrix(cmatrix &cmX, wavefunction &evals, grid g, cmatrix &cmR, long reduced_size)
{

   cmatrix res(reduced_size,reduced_size);
   res=cmR.transp()*(*this)*cmR;  
   res.lap_dsyev(evals,cmX,0,reduced_size-1);
   return res;
}

cmatrix cmatrix::build_orth_from_hermitian_matrix(cmatrix &cmX, wavefunction &evals, grid g, cmatrix &cmR, long reduced_size)
{

   cmatrix res(reduced_size,reduced_size);
   res=cmR.transp()*(*this)*cmR;  
   res.lap_zheev(evals,cmX,0,reduced_size-1);
   return res;
}

cmatrix cmatrix::build_orth_from_general_matrix(cmatrix &cmX, cmatrix &cmXinv, wavefunction &evals, grid g, cmatrix &cmR, long reduced_size)
{

   cmatrix res(reduced_size,reduced_size);
   cmatrix cmLX(reduced_size,reduced_size);
   cmatrix cmRX(reduced_size,reduced_size);
   wavefunction diagD(reduced_size);
   res=cmR.transp()*(*this)*cmR;  
   res.lap_zgeev(evals,cmLX,cmRX,0,reduced_size-1);
   diagD=(cmLX.adjoint()*cmRX).cm_pick_diag();
   cmX=cmRX.mult_cm_with_diag((diagD.invert()).sqrtwf());
   cmXinv=(cmLX.adjoint()).mult_diag_with_cm((diagD.invert()).sqrtwf());
   return res;
}


void cmatrix::inverse(cmatrix &cmY, cmatrix &cmYinv, wavefunction &evals)
{
   (*this)=cmY.mult_cm_with_diag(evals.invert())*cmYinv;
}



void cmatrix::build_atomic_hamiltonian_in_sturmian_basis(
   grid g, 
   double alpha, 
   double theta, 
   double charge, 
   cmatrix &cmS, 
   long reduced_size
   )
{
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long l,n,index;
   long size=g.nlindex(g.ngps_x()-1,g.ngps_y()-1)+1;
   complex imagi(0.0,1.0);
   complex beta=exp(-imagi*theta)/alpha;

   wavefunction V(size);

   // the nuclear potential
   for (l=0; l<=max_l; l++)
   {
      for (n=1+l; n<=max_n; n++)
      { 
	 index=g.nlindex(n-1,l);
	 V[index]=-beta*charge/((double)(n));
      };
   };

   (*this)=-0.5*beta*beta*cmS+beta*beta+V;
}

cmatrix cmatrix::build_kinetic_energy_operator(cmatrix &cmYorth, cmatrix &cmYorthinv, wavefunction &Haevals, grid g, double alpha, double theta, cmatrix &cmR, long reduced_size)
{
   long max_n=g.ngps_x();
   long max_l=g.ngps_y()-1;
   long l,ll,n,nn,index,iindex;
   long size=g.nlindex(g.ngps_x()-1,g.ngps_y()-1)+1;
   complex imagi(0.0,1.0);
   complex beta=exp(-imagi*theta)/alpha;

   wavefunction diagD(reduced_size);
   cmatrix res(reduced_size,reduced_size);
   cmatrix cmH(size,size);
   cmatrix cmLX(reduced_size,reduced_size);
   cmatrix cmRX(reduced_size,reduced_size);

   cmH=-0.5*beta*beta*(*this)+beta*beta;
   res=cmR.transp()*cmH*cmR; 


   res.lap_zgeev(Haevals,cmLX,cmRX,0,reduced_size-1);
   diagD=(cmLX.adjoint()*cmRX).cm_pick_diag();
   cmYorth=cmRX.mult_cm_with_diag((diagD.invert()).sqrtwf());
   cmYorthinv=(cmLX.adjoint()).mult_diag_with_cm((diagD.invert()).sqrtwf());


   return res;

}




void cmatrix::build_kinetic_orthogonal_propagator(grid g, complex timestep, cmatrix &cmY, cmatrix &cmYinv, wavefunction &Haevals)
{
   complex imagi(0.0,1.0);
   long size=g.nlindex(g.ngps_x()-1,g.ngps_y()-1)+1;

   (*this)=cmY.mult_cm_with_diag(Haevals.expwf(-imagi*timestep))*cmYinv;
}

void cmatrix::build_rotation_matrix_from_orthogonal_rdr_in_sturmian_representation(
   grid g,
   long size,
   double theta, 
   cmatrix &cmY, 
   cmatrix &cmYinv, 
   wavefunction &evals,
   cmatrix &cmS,
   cmatrix &cmR
   )
{
   complex imagi(0.0,1.0);
   cmatrix orthres(size,size);
   orthres=cmY.mult_cm_with_diag(evals.expwf(imagi*theta))*cmYinv;
   (*this)=cmR*orthres*cmR.transp()*cmS;
}

void cmatrix::build_rotation_matrix_from_orthogonal_rdr_in_orthogonal_representation(
   grid g,
   long size,
   double theta, 
   cmatrix &cmY, 
   cmatrix &cmYinv, 
   wavefunction &evals
   )
{
   complex imagi(0.0,1.0);
   (*this)=cmY.mult_cm_with_diag(evals.expwf(imagi*theta))*cmYinv;
}

void cmatrix::build_atomic_propagator_from_orthogonal_hamiltonian(
   grid g, 
   complex timestep, 
   cmatrix &cmY, 
   cmatrix &cmYinv, 
   wavefunction &Haevals
   )
{
   complex imagi(0.0,1.0);
   (*this)=cmY.mult_cm_with_diag(Haevals.expwf(-imagi*timestep))*cmYinv;
}


void cmatrix::build_interaction_propagator_in_length_gauge(
   grid g, 
   complex timestep, 
   cmatrix &cmY, 
   cmatrix &cmYinv,
   wavefunction &evals, 
   double fieldampl,
   double alpha,
   double theta
   )
{
   complex imagi(0.0,1.0);
   (*this)=cmY.mult_cm_with_diag(evals.expwf(-imagi*fieldampl*timestep*0.5*alpha*exp(imagi*theta)))*cmYinv;
}


void cmatrix::build_interaction_propagator_in_velocity_gauge(
   grid g, 
   complex timestep, 
   cmatrix &cmY, 
   cmatrix &cmYinv,
   wavefunction &evals, 
   double fieldampl,
   double alpha,
   double theta,
   double time
   )
{
   complex imagi(0.0,1.0);
   (*this)=cmY.mult_cm_with_diag(evals.expwf(-imagi*(-fieldampl*time)*timestep*0.5/alpha*exp(-imagi*theta)))*cmYinv;
}

void cmatrix::build_gauge_transformation_operator(
   grid g, 
   cmatrix &cmY, 
   cmatrix &cmYinv,
   wavefunction &evals, 
   double fieldampl,
   double alpha,
   double theta,
   double time
   )
{
   complex imagi(0.0,1.0);
   (*this)=cmY.mult_cm_with_diag(evals.expwf(imagi*(-fieldampl*time)*0.5*alpha*exp(imagi*theta)))*cmYinv;
}





