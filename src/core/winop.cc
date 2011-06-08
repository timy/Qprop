#include<winop.h>

/*!

 */

int winop_fullchi(wavefunction &fullchi, wavefunction &spec_result_lsub,
		  complex *spec_result_tot, double energ, double energ_width,
		  wavefunction staticpot, fluid V_ee_0,
		  double nuclear_charge, grid g,wavefunction wf,int iv)
{
   long first_ang, last_ang, l_index, i, m_index, m_limit, index, lmindex;
   int modifiedcornerelements, ii, jj, tt;
   complex spec_result_part;

   spec_result_part = complex(0.0, 0.0);

   FILE *file_dbg;

   modifiedcornerelements=1;

   wavefunction wf_of_interest(g.ngps_x());
   wavefunction vec_pot(g.ngps_x());
   wavefunction chi(g.ngps_x());
   wavefunction newchi(g.ngps_x());
   wavefunction m2_wf_of_in(g.ngps_x());
   wavefunction m2_chi(g.ngps_x());
   wavefunction a(g.ngps_x());
   wavefunction b(g.ngps_x());
   wavefunction c(g.ngps_x());
  
   // *****Upper left corrections******
   double const m2_b = -10.0/6.0;
   double const d2_b = -2.0/(g.delt_x()*g.delt_x());

   double const m2_a = -1.0/6.0;
   double const d2_a = 1.0/(g.delt_x()*g.delt_x());

   double const m2_c = m2_a;
   double const d2_c = d2_a;
    
   double const d2_0 = -2.0/(g.delt_x()*g.delt_x())
      *(1.0 - nuclear_charge*g.delt_x()/(12.0 - 10.0*nuclear_charge*g.delt_x()));
   double const m2_0 = -2.0*(1.0 + g.delt_x()*g.delt_x()/12.0*d2_0);

   *spec_result_tot = complex(0.0, 0.0);
   spec_result_lsub.nullify();

   if (iv==1) fprintf(stdout, " energy         l   m     P \n");

   for (l_index=0; l_index<g.ngps_y(); l_index++)
   {
      if (g.dimens()==34) 
      {
	 m_limit = 0;
	 spec_result_lsub[l_index] = complex(0.0,0.0);
      } 
      else 
      {
	 m_limit = l_index;
      }

      for (m_index=-m_limit; m_index<=m_limit; m_index++) 
      {

	 if (g.dimens()==44)
	 {
	    spec_result_lsub[((l_index+1)*(l_index+1)-(l_index+1)+m_index)] 
	       = complex(0.0,0.0);
	 }; 


	 for (i=0; i<g.ngps_x(); i++)
	 { index=g.index(i, l_index, m_index, 0); wf_of_interest[i]=wf[index];}

//     if ((l_index==1) && (m_index==-l_index)) 
//     {
//       file_dbg = fopen("wf_of_interest.dat", "w");
//       wf_of_interest.dump_to_file_sh(g, file_dbg, 1, iv);
//       fclose(file_dbg);
//       exit(-104);
//     }

	 //Multiply M2 and wf;
	 m2_wf_of_in[0] = m2_b*wf_of_interest[0] + m2_a*wf_of_interest[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_wf_of_in[0] = m2_0*wf_of_interest[0] + m2_a*wf_of_interest[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_wf_of_in[ii] = m2_c*wf_of_interest[ii-1] + 
	       m2_b*wf_of_interest[ii] + m2_a*wf_of_interest[ii+1];

	 m2_wf_of_in[g.ngps_x()-1] = m2_c*wf_of_interest[g.ngps_x()-2]
	    + m2_b*wf_of_interest[g.ngps_x()-1];

	 // build up the first matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index = g.index(ii,l_index,m_index,0);
	    vec_pot[ii] = staticpot[index] + V_ee_0[ii] - energ -
	       exp(complex(0.0,1.0)*0.875*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++) a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)   b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)   c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];

	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 // Solve first matrix

	 chi.solve_du(a,b,c, m2_wf_of_in, g.ngps_x());
 
	 // Multiply M2 and chi;

	 m2_chi[0] = m2_b*chi[0] + m2_a*chi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_chi[0] = m2_0*chi[0] + m2_a*chi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_chi[ii] = m2_c*chi[ii-1] + m2_b*chi[ii] + m2_a*chi[ii+1];

	 m2_chi[g.ngps_x()-1] = m2_c*chi[g.ngps_x()-2] + m2_b*chi[g.ngps_x()-1];


	 // build up the second matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii] = staticpot[index] + V_ee_0[ii] - energ
	       + exp(complex(0.0,1.0)*0.875*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++) a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)   b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)   c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 // Solve second matrix

	 newchi.solve_du(a,b,c,m2_chi,g.ngps_x());

	 // and again ...

	 //Multiply M2 and newchi;
	 m2_wf_of_in[0] = m2_b*newchi[0] + m2_a*newchi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_wf_of_in[0] = m2_0*newchi[0] + m2_a*newchi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_wf_of_in[ii] = m2_c*newchi[ii-1] + m2_b*newchi[ii] + m2_a*newchi[ii+1];

	 m2_wf_of_in[g.ngps_x()-1] = m2_c*newchi[g.ngps_x()-2] + 
	    m2_b*newchi[g.ngps_x()-1];

	 // build up the first matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii]= staticpot[index]
	       + V_ee_0[ii] - energ - exp(complex(0.0,1.0)*0.375*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++) a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)   b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)   c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];

	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 //Solve first matrix

	 chi.solve_du(a,b,c, m2_wf_of_in, g.ngps_x());

	 //Multiply M2 and chi;

	 m2_chi[0] = m2_b*chi[0] + m2_a*chi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_chi[0] = m2_0*chi[0] + m2_a*chi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_chi[ii] = m2_c*chi[ii-1] + m2_b*chi[ii] + m2_a*chi[ii+1];

	 m2_chi[g.ngps_x()-1] = m2_c*chi[g.ngps_x()-2] + m2_b*chi[g.ngps_x()-1];

	 // build up the second matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {

	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii]= staticpot[index]
	       + V_ee_0[ii] - energ + exp(complex(0.0,1.0)*0.375*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++) a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)   b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)   c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 //Solve second matrix

	 newchi.solve_du(a,b,c,m2_chi,g.ngps_x());

	 // and again ...

	 //Multiply M2 and newchi;
	 m2_wf_of_in[0] = m2_b*newchi[0] + m2_a*newchi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_wf_of_in[0] = m2_0*newchi[0] + m2_a*newchi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_wf_of_in[ii] = m2_c*newchi[ii-1] + m2_b*newchi[ii] + m2_a*newchi[ii+1];

	 m2_wf_of_in[g.ngps_x()-1] = m2_c*newchi[g.ngps_x()-2] + m2_b*newchi[g.ngps_x()-1];

	 // build up the first matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii] = staticpot[index] + V_ee_0[ii] -
	       energ - exp(complex(0.0,1.0)*0.625*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++)  a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)    b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)    c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];

	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 //Solve first matrix

	 chi.solve_du(a,b,c, m2_wf_of_in, g.ngps_x());

	 //Multiply M2 and chi;

	 m2_chi[0] = m2_b*chi[0] + m2_a*chi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_chi[0] = m2_0*chi[0] + m2_a*chi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_chi[ii] = m2_c*chi[ii-1] + m2_b*chi[ii] + m2_a*chi[ii+1];

	 m2_chi[g.ngps_x()-1] = m2_c*chi[g.ngps_x()-2] + m2_b*chi[g.ngps_x()-1];


	 // build up the second matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii]= staticpot[index] + V_ee_0[ii]
	       - energ + exp(complex(0.0,1.0)*0.625*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++)  a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)    b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)    c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 //Solve second matrix

	 newchi.solve_du(a,b,c,m2_chi,g.ngps_x());

	 // and again ...

	 //Multiply M2 and newchi;
	 m2_wf_of_in[0] = m2_b*newchi[0] + m2_a*newchi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_wf_of_in[0] = m2_0*newchi[0] + m2_a*newchi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_wf_of_in[ii] = m2_c*newchi[ii-1] + m2_b*newchi[ii] + m2_a*newchi[ii+1];

	 m2_wf_of_in[g.ngps_x()-1] = m2_c*newchi[g.ngps_x()-2] + m2_b*newchi[g.ngps_x()-1];

	 // build up the first matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii]= staticpot[index]
	       + V_ee_0[ii] - energ - exp(complex(0.0,1.0)*0.125*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++)  a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)    b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)    c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];

	 if ((l_index==0) && (modifiedcornerelements==1))
	    b[0] = d2_0 + m2_0*vec_pot[0];

	 //Solve first matrix

	 chi.solve_du(a,b,c, m2_wf_of_in, g.ngps_x());

	 //Multiply M2 and chi;

	 m2_chi[0] = m2_b*chi[0] + m2_a*chi[1];
	 if ((l_index==0) && (modifiedcornerelements==1))
	    m2_chi[0] = m2_0*chi[0] + m2_a*chi[1];

	 for (ii=1; ii<g.ngps_x()-1; ii++)
	    m2_chi[ii] = m2_c*chi[ii-1] + m2_b*chi[ii] + m2_a*chi[ii+1];

	 m2_chi[g.ngps_x()-1] = m2_c*chi[g.ngps_x()-2] + m2_b*chi[g.ngps_x()-1];

	 // build up the second matrix

	 for (ii=0; ii<g.ngps_x(); ii++)
	 {
	    index=g.index(ii,l_index,m_index,0);		  
	    vec_pot[ii]= staticpot[index] + V_ee_0[ii]
	       - energ + exp(complex(0.0,1.0)*0.125*M_PI)*energ_width;
	 }

	 for (ii=0; ii<g.ngps_x()-1; ii++)  a[ii] = d2_a + m2_a*vec_pot[ii+1];
	 for (ii=1; ii<g.ngps_x(); ii++)    b[ii] = d2_b + m2_b*vec_pot[ii];
	 for (ii=1; ii<g.ngps_x(); ii++)    c[ii] = d2_c + m2_c*vec_pot[ii-1];

	 b[0] = d2_b + m2_b*vec_pot[0];
	 if ((l_index==0) && (modifiedcornerelements==1)) b[0]=d2_0+m2_0*vec_pot[0];

	 //Solve second matrix

	 newchi.solve_du(a,b,c,m2_chi,g.ngps_x());

	 newchi*=energ_width*energ_width*energ_width*energ_width*
	    energ_width*energ_width*energ_width*energ_width;

	 //
	 // Embedding
	 //
	 fullchi.embed_as_x(g, l_index, m_index, 0, newchi);

	 spec_result_part=newchi*newchi;
	 spec_result_part*=g.delt_x();

	 if (g.dimens()==44) 
	 {
	    lmindex=((l_index+1)*(l_index+1)-(l_index+1)+m_index);
	 }
	 else
	 {
	    lmindex=l_index;
	 };

	 spec_result_lsub[lmindex]+=spec_result_part;	


	 *spec_result_tot+=spec_result_part;

	 if (iv==1) fprintf(stdout, "% le %3ld %3ld  (% le, % le)\n",
			    energ, l_index, m_index,
			    real(spec_result_lsub[lmindex]),
			    imag(spec_result_lsub[lmindex]));
      }
   }

   if (iv==1)
      std::cout << "energ: " << energ << *spec_result_tot << std::endl;

   return(0);

}
