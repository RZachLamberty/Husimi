#include"ed.h"
#include"printing_functions.h"
/////////////////////////////////////////////////////////////////////////////////////

double measure_sz(std::vector<int> const &config)
{ 
  int 		n_up,n_down;
  n_up=(int) count(config.begin(),config.end(),1);
  n_down=config.size()-n_up;
  return double(n_up-n_down)/2.0;
}

/////////////////////////////////////////////////////////////////////////////////////
double compute_spin(double s_s_plus_one)
{return ((sqrt((4.0*s_s_plus_one)+1.0)-1.0)/2.0);}

/////////////////////////////////////////////////////////////////////////////////////

void ed_get_eigs(Ham &h, std::vector<double> &eigs)
{
   
   time_t 				start,end;
   int 					size=pow(2,h.num_sites);
   int 					i,j,k,l,m;
   int 					spin;
   double 				dif;
   double 				tmp;
   complex<double> 			tot,c_first_second;
   std::vector<int> 			config(h.num_sites);
   std::vector< complex<double> > 	hints_list;
   std::vector< std::vector<int> > 	touched_sites_list,vals_on_touched_list;
   std::string 				tmp_string;
   std::string 				eigs_dump_filename;
   std::string 				evecs_dump_filename;
   Matrix 				h_mat(size,size);
   Matrix 				eigenvecs(size,size);
   
   eigs.resize(size);
   
   for (i=0;i<size;i++)
   {
        convert_num_to_vec(i,2,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (k=0;k<touched_sites_list.size();k++)
        {
            j=i;
            for (l=0;l<touched_sites_list[k].size();l++)
            {j=j-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
            
            if (j>=i){h_mat(i,j)=h_mat(i,j)+real(hints_list[k]);h_mat(j,i)=h_mat(i,j);}
        }
   }
   
   time (&start);
   symmetric_diagonalize(h_mat,eigs,eigenvecs);
   time (&end);
   dif=difftime(end,start);

   cout<<"==================================================================="<<endl;
   cout<<"Total time to diagonalize (Exact) was "<<dif<<" seconds"<<endl;
   cout<<"==================================================================="<<endl;
   
   for (j=0;j<size;j++)
   {
      i=0;  	
      while (abs(eigenvecs(i,j))<1.0e-2) {i+=1;}

      //cout<<"Found non zero element i = "<<i<<endl;
        
      convert_num_to_vec(i,2,h.num_sites,config);
      tot=0.0;

      for (int first=0;first<h.num_sites;first++)
      {
       	  for (int second=first+1;second<h.num_sites;second++)
          {  
		touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
		calc_hints_sxsx_sysy(1.0,first,second,config,
			             touched_sites_list,
                                     vals_on_touched_list,
                                     hints_list);
		
		c_first_second=0.0;

		for (k=0;k<touched_sites_list.size();k++)
		{
			m=i;
			for (l=0;l<touched_sites_list[k].size();l++){m=m-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
			c_first_second+=(hints_list[k]*eigenvecs(m,j));
		}
			
		touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
      	  	calc_hints_szsz(1.0,first,second,config,
		      touched_sites_list,vals_on_touched_list,hints_list);
 		
		c_first_second+=(hints_list[0]*eigenvecs(i,j));

		tot+=c_first_second;	  
	}
      }		
     
      tot=tot/eigenvecs(i,j); 
      tmp=(2.0*real(tot))+(0.75*double(h.num_sites));
      //cout<<"S(S+1) = "<<tmp<<endl;
      cout<<"==========================================="<<endl;
      cout<<"Spin   = "<<compute_spin(tmp)<<endl;
      cout<<"S_z    = "<<measure_sz(config)<<endl;
      cout<<"Energy = "<<eigs[j]<<endl;	
      cout<<"==========================================="<<endl;
  }

}

//////////////////////////////////////////////////////////////////////////////
void lanczos_no_sym_get_eigs(Ham &h,
                             int iterations, 
                             std::vector<double> &eigs)
{
   int i,size=pow(2,h.num_sites);
   iterations=min(iterations,size);
   std::vector<int> map,inverse_map;
   std::vector<double> spins;

   for (i=0;i<size;i++){ map.push_back(i);inverse_map.push_back(i);}
   
   lanczos_given_map(h,iterations,map,inverse_map,eigs,spins);
}
//////////////////////////////////////////////////////////////////////////////
void lanczos_requested_sz(Ham &h,
                      int iterations, 
                      std::vector<double> &eigs,
		      double s_z,
		      std::vector<double> &spins,
                      bool measure_s)
{
   int 				i,j;
   int 				n,n_up,n_down;
   int 				size=pow(2,h.num_sites);
   iterations=min(iterations,size);
   std::vector<int> 		map;
   std::vector<int> 		inverse_map(size);
   std::vector<int> 		config(h.num_sites);   
   bool 			change;

   // Initialize
   eigs.clear();spins.clear();
   
   cout<<"==========================================="<<endl;
   cout<<"Starting to make maps"<<endl;
   cout<<"==========================================="<<endl;

   // Make maps
   for (i=0;i<size;i++)
   {
        convert_num_to_vec(i,2,h.num_sites,config);
  	n_up=(int) count(config.begin(),config.end(),1);
	n_down=h.num_sites-n_up;
	if (abs((double(n_up-n_down)/2.0)-s_z)<1.0e-10) 
	{
   		map.push_back(i);
        	inverse_map[i]=map.size()-1;
  	} 
   }

   cout<<"==========================================="<<endl;
   cout<<"Finished making map in lanczos_requested_sz"<<endl;
   cout<<"==========================================="<<endl;
   
   lanczos_given_map(h,iterations,map,inverse_map,
	             eigs,spins,measure_s);

   change=true;
   n=0;
   // Sort (eigs,spins,szs)
   while (change)
   {	
	   n+=1; 
	   change=false;		
	   for (i=0;i<eigs.size()-n;i++)
	   {
		if (eigs[i]>eigs[i+1])
		{
		   change=true;
		   swap(eigs[i],eigs[i+1]);
		   if (measure_s){swap(spins[i],spins[i+1]);} 	
		}
		
		if (eigs[i]==eigs[i+1])
		{
		  if (measure_s)
		  {
			if (spins[i]>spins[i+1])
			{
				change=true;
				swap(spins[i],spins[i+1]);	
			}
		  }
		}
	   }
   }
}


//////////////////////////////////////////////////////////////////////////////
void lanczos_spin_sym(Ham &h,
                      int iterations, 
                      std::vector<double> &eigs,
		      std::vector<double> &szs,
		      std::vector<double> &spins,
                      bool measure_s, bool ipr, 
		      bool only_sz_non_negative)
{
   int 					i,j;
   int 					n,n_up,n_down;
   int 					n_up_max=0;
   int 					size=pow(2,h.num_sites);
   iterations=min(iterations,size);
   double 				sz;
   std::vector< std::vector<int> > 	maps(h.num_sites+1);
   std::vector<int> 			inverse_map(size);
   std::vector<int> 			config(h.num_sites);   
   std::vector<double> 			spins_sz;
   std::vector<double> 			eigs_sz;
   bool 				change;
   bool                                 allowed;
   // Initialize
   eigs.clear();spins.clear();szs.clear();

   // Make maps
   for (i=0;i<size;i++)
   {
        convert_num_to_vec(i,2,h.num_sites,config);
  	n_up=(int) count(config.begin(),config.end(),1);
   	maps[n_up].push_back(i);
	if (n_up>n_up_max) {n_up_max=n_up;}
        inverse_map[i]=maps[n_up].size()-1;
   }

   for (i=0;i<=n_up_max;i++) 
   {
     	     n_down=h.num_sites-i;
     	     sz=double(i-n_down)/2.0;
    
	     allowed=true;
             if (only_sz_non_negative)
	     {
		if (sz>-1.0e-10) {allowed=true;}
		else {allowed=false;}
	     } 
	
             if (allowed)
	     {
		     if (ipr) {cout<<"Lanczos for S_z = "<<sz<<endl;}
		     lanczos_given_map(h,iterations,
				       maps[i],inverse_map,
				       eigs_sz,spins_sz,measure_s,ipr);

		    eigs.insert(eigs.end(),eigs_sz.begin(),eigs_sz.end());
		    if (measure_s) spins.insert(spins.end(),spins_sz.begin(),spins_sz.end());
		    for (j=0;j<eigs_sz.size();j++){szs.push_back(sz);}
		    eigs_sz.clear();spins_sz.clear();	
	     }
   }

   if (ipr)
   {
	   cout<<"==========================================="<<endl;
	   cout<<"Done diagonalizing.... now sorting"<<endl;
	   cout<<"==========================================="<<endl;
   }
   //cout<<"Eigs"<<endl;print_vec(eigs);cout<<"Szs"<<endl;print_vec(szs);
 
   change=true;
   n=0;
   // Sort (eigs,spins,szs)
   /*while (change)
   {	
	   n+=1; 
	   change=false;		
	   for (i=0;i<eigs.size()-n;i++)
	   {
		if (eigs[i]>eigs[i+1])
		{
		   change=true;
		   swap(eigs[i],eigs[i+1]);
		   swap(szs[i],szs[i+1]);
		   if (measure_s){swap(spins[i],spins[i+1]);} 	
		}
		
		if (eigs[i]==eigs[i+1])
		{
		  if (measure_s)
		  {
			if (spins[i]>spins[i+1])
			{
				change=true;
				swap(spins[i],spins[i+1]);	
				swap(szs[i],szs[i+1]);
			}
			if (spins[i]==spins[i+1]) {if (szs[i]>szs[i+1]) {change=true;swap(szs[i],szs[i+1]);} }
		  }
		  else
		  {if (szs[i]>szs[i+1]) {change=true;swap(szs[i],szs[i+1]);}}  	
		}
	   }
   }i*/
   for (j=0;j<eigs.size()-1;j++)
   {	
	   for (i=j+1;i<eigs.size();i++)
	   {
		if (eigs[j]>eigs[i])
		{
		   swap(eigs[j],eigs[i]);
		   swap(szs[j],szs[i]);
		   if (measure_s){swap(spins[j],spins[i]);} 	
		}
		
		if (abs(eigs[j]-eigs[i])<1.0e-8)
		{
		  if (measure_s)
		  {
			if (spins[j]>spins[i])
			{
				swap(spins[j],spins[i]);	
				swap(szs[j],szs[i]);
			}
			if (spins[j]==spins[i]) {if (szs[j]>szs[i]) {swap(szs[j],szs[i]);} }
		  }
		  else
		  {if (szs[j]>szs[i]) {swap(szs[j],szs[i]);}}  	
		}
	   }
   }
}

//////////////////////////////////////////////////////////////////////////////
void lanczos_given_map(Ham &h,
                      int iterations, 
		      std::vector<int> const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &spins,
                      bool measure_s, bool ipr)
{
   time_t 	      				start,end;
   int 		      				i,it,j,k,l,m;
   int 		      				size=map.size();
   double 	      				spin;
   bool 	      				orth_failed;
   iterations=min(iterations,size);
   eigs.resize(iterations);
   spins.clear();
   double 					dif,tmp;
   double					q,alpha,beta,norm;
   complex<double> 				c_first_second,tot;
   std::vector<int> 				config(h.num_sites),new_configs;
   std::vector<double> 				alphas,betas;
   std::vector<double> 				h_dot_v(size),w(size);
   std::vector<double> 				v_p(size),v_o(size);
   std::vector<double> 				v_p_old(size);
   std::vector<double>				ritz_eigenvec(size);
   std::vector< complex<double> > 		hints_list;
   std::vector< std::vector<int> > 		touched_sites_list,vals_on_touched_list;
   std::vector< std::vector<int> > 		vec_new_configs;
   std::vector< std::vector<double> > 		vs;
   std::vector< std::vector< complex<double> > > vec_hints_list;
   Matrix 					t_mat(iterations,iterations);
   Matrix					t_eigenvecs(iterations,iterations);
 
   //ipr=true;
 
   //================= 
   // Initializations
   //================= 
   if (ipr) cout<<"Making Hamiltonian" <<endl;
 
   for (i=0;i<size;i++)
   {
        v_p[i]=uniform_rnd();
        v_o[i]=0.0;w[i]=0.0;
        
        convert_num_to_vec(map[i],2,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
	new_configs.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (k=0;k<touched_sites_list.size();k++)
        {
		j=map[i];
                for (l=0;l<touched_sites_list[k].size();l++)
		{j=j-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
		new_configs.push_back(inverse_map[j]);
        }
       
        vec_new_configs.push_back(new_configs);     // New configs 
        vec_hints_list.push_back(hints_list);       // Hints
   }
   if (ipr) cout<<"TRACE: Finished making Hamiltonian" <<endl;
       
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;
   betas.push_back(beta);

   //vs.push_back(v_p);

   if (ipr)
   {
   	cout<<"iterations = "<<iterations<<endl;
   	cout<<"vs.size()  = "<<vs.size()<<endl;
   }
 
   for (it=0;it<iterations;it++)
   {
       if (ipr)
       {
       		cout<<"================================================================="<<endl;
       		cout<<"Doing Lanczos iteration = "<<it<<endl;
       }

       if (h.num_sites>24) if (vs.size()>32){vs.erase(vs.begin());}
       
       vs.push_back(v_p);     
       time (&start);
       for (i=0;i<size;i++) // Computing H*v_p - This is the bulk of the operation
       {
            for (k=0;k<vec_new_configs[i].size();k++)
            {w[vec_new_configs[i][k]]+=(real(vec_hints_list[i][k])*v_p[i]);}
       }
       
       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
       alphas.push_back(alpha);
       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
       v_o=v_p;
       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
       v_p=w;
       dscal(size,1.0/beta,&*v_p.begin(),1);
       betas.push_back(beta);
       dscal(size,0.0,&*w.begin(),1);

       // Now reorthogonalize vectors
       if (vs.size()<iterations)
       {
	       orth_failed=false;
	       for (i=0;i<vs.size();i++)
	       {
		    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
		    if (abs(q)>1.0e-10)
		    {
			i=vs.size();
			if (ipr)
			{
				cout<<"q (overlap) ="<<q<<endl;
				cout<<"--------------------------------------------------"<<endl;
				cout<<"Orthogonalization failed... choosing random vector"<<endl;
				cout<<"--------------------------------------------------"<<endl;
			}
			orth_failed=true;
		    }
		    else
		    {
			daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
			dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
		    }
	       }
	       
	       norm=0.0;
	       if (orth_failed)
	       {
			while (abs(norm)<1.0e-2)
			{
				for (j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
				norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
				dscal(size,1.0/norm,&*v_p_old.begin(),1);
				//cout<<"norm="<<norm<<endl;
				v_p=v_p_old;
				for (i=0;i<vs.size();i++)
				{
					q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
					//cout<<"q (after orth fail )="<<q<<endl;
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
				}
				norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
			}
			dscal(size,1.0/norm,&*v_p.begin(),1);
	       }
       }

       //vs.push_back(v_p);     
       time (&end);
       dif=difftime(end,start);
       
       if (ipr)
       {
       		cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
       		cout<<"================================================================="<<endl;
       }
   }

   //if (ipr) {cout<<"vs.size() = "<<vs.size()<<endl;}
   //cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
   
   if (ipr) {cout<<"Time to build  T was "<<dif<<" seconds"<<endl;}
 
   //cout<<"Iterations="<<iterations<<endl; 
   
   time (&start);
   symmetric_diagonalize(t_mat,eigs,t_eigenvecs);
   time (&end);
   dif=difftime(end,start);

   if (ipr) {cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}

   if (measure_s)
   {
	   for (j=0;j<iterations;j++)
	   {
	      if (ipr) {cout<<"Making j= "<<j<<" Ritz eigenvector"<<endl;}

	      for (i=0;i<size;i++)
	      {
		ritz_eigenvec[i]=0.0;
		for (k=0;k<iterations;k++){ritz_eigenvec[i]+=vs[k][i]*t_eigenvecs(k,j);}
	      }
	      
	      //cout<<"Finished computing Ritz eigenvector"<<endl;
	       
	      i=0;  	
	      while (abs(ritz_eigenvec[i])<1.0e-3 and i<map.size()) {i+=1;}

	      //cout<<"Found non zero element i = "<<i<<endl;
	      //cout<<"map.size()="<<map.size()<<endl;

	      convert_num_to_vec(map[i],2,h.num_sites,config);
	      tot=0.0;

	      for (int first=0;first<h.num_sites;first++)
	      {
		  for (int second=first+1;second<h.num_sites;second++)
		  {  
			touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
			calc_hints_sxsx_sysy(1.0,first,second,config,
					     touched_sites_list,
					     vals_on_touched_list,
					     hints_list);
			
			c_first_second=0.0;

			for (k=0;k<touched_sites_list.size();k++)
			{
				m=map[i];
				for (l=0;l<touched_sites_list[k].size();l++){m=m-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
				c_first_second+=(hints_list[k]*ritz_eigenvec[inverse_map[m]]);
			}
				
			touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
			calc_hints_szsz(1.0,first,second,config,
			      touched_sites_list,vals_on_touched_list,hints_list);
			
			c_first_second+=(hints_list[0]*ritz_eigenvec[i]);

			tot+=c_first_second;	  
		}
	      }		
	     
	      tot=tot/ritz_eigenvec[i]; 
	      tmp=(2.0*real(tot))+(0.75*double(h.num_sites));
	      //cout<<"S(S+1) = "<<tmp<<endl;
	      spin=compute_spin(tmp);
	      spins.push_back(spin);
       	     
	      if (ipr)
	      {
		      cout<<"===================================="<<endl;
		      cout<<"Spin   = "<<spin<<endl;
		      cout<<"S_z    = "<<measure_sz(config)<<endl;
		      cout<<"Energy = "<<eigs[j]<<endl;	
		      cout<<"===================================="<<endl;
	      }
	}

    }

}

//////////////////////////////////////////////////////////////////////////////
void ed_with_hints_given(std::vector< std::vector<int> > 		const &map,
		      	 std::vector< std::vector<double> > 	        const &hints,
                      	 std::vector<double> 			        &eigs,
			 Matrix 					&eigenvecs,
			 bool 					        ipr)
{
   time_t 				start,end;
   time_t 				start_0,end_0;
   int 					i,k;
   int 					size=map.size();
   Matrix                               ham(size,size);
   double 				dif;
 
   cout<<" Number of states in Matrix diag is "<<size<<endl;
   for (i=0;i<size*size;i++) ham[i]=0.0;
   time(&start_0);
   eigs.clear();
   eigenvecs.clear();
   eigs.resize(size);
   eigenvecs.resize(size,size);    
   time (&start_0);
   for (i=0;i<size;i++) 
   {
   	for (k=0;k<map[i].size();k++){ham(i,map[i][k])+=(hints[i][k]);}
   }
   
   symmetric_diagonalize(ham,eigs,eigenvecs);
   time(&end_0);
   dif=difftime(end_0,start_0);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void lanczos_with_hints_given(int 			      		iterations,
			      int 			      		how_many_eigenvecs, 
		      	      std::vector< std::vector<int> > 		const &map,
		      	      std::vector< std::vector<double> > 	const &hints,
                      	      std::vector<double> 			&eigs,
			      Matrix 					&eigenvecs,
			      bool 					ipr)
{
   time_t 				start,end;
   time_t 				start_0,end_0;
   int 					i,it,j,k,l,m;
   int 					size=map.size();
   
   time(&start_0);
   iterations=min(iterations,size);
   eigs.resize(iterations);
   
   double 				dif,tmp;
   double				q,alpha,beta,norm;
   std::vector<double> 			alphas,betas;
   std::vector<double> 			h_dot_v(size),w(size);
   std::vector<double> 			v_p(size),v_o(size);
   std::vector< std::vector<double> >   vs;
   Matrix 				t_mat(iterations,iterations);
   Matrix				t_eigenvecs(iterations,iterations);

   //vs.clear();
   ipr=true;
   // Initializations
   for (i=0;i<size;i++) {v_p[i]=uniform_rnd();v_o[i]=0.0;w[i]=0.0;}
       
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;betas.push_back(beta);

   for (it=0;it<iterations;it++)
   {
       //if (ipr) cout<<"Doing Lanczos iteration = "<<it<<endl;
       if (vs.size()<size) vs.push_back(v_p);     
       
       time (&start);
       for (i=0;i<size;i++) 
       // Computing H*v_p - This is the bulk of the operation
       {
	      for (k=0;k<map[i].size();k++)
	      {w[map[i][k]]+=(hints[i][k]*v_p[i]);}
       }
       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
       alphas.push_back(alpha);
       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
       v_o=v_p;
       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
       v_p=w;
       dscal(size,1.0/beta,&*v_p.begin(),1);
       betas.push_back(beta);
       dscal(size,0.0,&*w.begin(),1);
       time (&end);
       dif=difftime(end,start);
       
       if (ipr) cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
   }

   if (ipr) cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (j=0;j<iterations*iterations;j++) t_mat[j]=0.0;

   for (j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
  
   if (ipr) cout<<"Time to build T was "<<dif<<" seconds"<<endl;
 
   time (&start);
   symmetric_diagonalize(t_mat,eigs,t_eigenvecs);
   //if (ipr) cout<<"Tmat"<<endl;
   //if (ipr) print_real_mat(t_mat);
   time (&end);
   dif=difftime(end,start);

   if (ipr) cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;
	  
   how_many_eigenvecs=min(size,how_many_eigenvecs);
   eigenvecs.resize(size,how_many_eigenvecs);

   //if (ipr) cout<<"how_many_eigenvecs="<<how_many_eigenvecs<<endl;//if (ipr) cout<<"size="<<size<<endl;

   for (j=0;j<how_many_eigenvecs;j++)
   {
	      if (ipr) cout<<"Making "<<j<<"th Ritz eigenvector"<<endl;
	      for (i=0;i<size;i++)
	      {
		eigenvecs(i,j)=0.0;
		for (k=0;k<iterations;k++){eigenvecs(i,j)+=vs[k][i]*t_eigenvecs(k,j);}
	      }
   }      
   //f (ipr) print_real_mat(eigenvecs);
   if (ipr) cout<<"Done computing Ritz eigenvectors"<<endl;
   time(&end_0);
   dif=difftime(end_0,start_0);
   if (ipr) cout<<"Total time for all Lanczos iterations was "<<dif<<" seconds"<<endl;
}
*/

void lanczos_with_hints_given(int 			      		iterations,
			      int 			      		how_many_eigenvecs, 
		      	      std::vector< std::vector<int> > 		const &map,
		      	      std::vector< std::vector<double> > 	const &hints,
                      	      std::vector<double> 			&eigs,
			      Matrix 					&eigenvecs,
			      bool 					ipr)
{
   time_t 				start,end;
   time_t 				start_0,end_0;
   int 					i,it,j,k,l,m;
   int 					size=map.size();
   bool 				orth_failed;
   
   time(&start_0);
   iterations=min(iterations,size);
   eigs.resize(iterations);
   
   double 				dif,tmp;
   double				q,alpha,beta,norm;
   std::vector<double> 			alphas,betas;
   std::vector<double> 			h_dot_v(size),w(size);
   std::vector<double> 			v_p(size),v_o(size),v_p_old(size),tmpv;
   std::vector< std::vector<double> >   vs;
   Matrix 				t_mat(iterations,iterations);
   Matrix				t_eigenvecs(iterations,iterations);

   /*   PRINT STATEMENT
   ipr=true;
   */
   // Initializations
   for (i=0;i<size;i++)
   {
        v_p[i]=uniform_rnd();
        v_o[i]=0.0;w[i]=0.0;
   }
       
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;betas.push_back(beta);

   for (it=0;it<iterations;it++)
   {
       //if (ipr) cout<<"Doing Lanczos iteration = "<<it<<endl;
       if (vs.size()<size) vs.push_back(v_p);     
       
       time (&start);
       for (i=0;i<size;i++) 
       // Computing H*v_p - This is the bulk of the operation
       {
	      for (k=0;k<map[i].size();k++)
	      {w[map[i][k]]+=(hints[i][k]*v_p[i]);}
       }
       
       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
       alphas.push_back(alpha);
       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
       v_o=v_p;
       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
       v_p=w;
       dscal(size,1.0/beta,&*v_p.begin(),1);
       betas.push_back(beta);
       dscal(size,0.0,&*w.begin(),1);

       // Now reorthogonalize vectors
       if (vs.size()<iterations)
       {
	       orth_failed=false;
	       v_p_old=v_p;
	       for (i=0;i<vs.size();i++)
	       {
		    q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
		    if (abs(q)>1.0e-10)
		    {
			i=vs.size();
			if (ipr)
			{
				cout<<"q (overlap) ="<<q<<endl;
				cout<<"--------------------------------------------------"<<endl;
				cout<<"Orthogonalization failed... choosing random vector"<<endl;
				cout<<"--------------------------------------------------"<<endl;
			}
			orth_failed=true;
		    }
		    else
		    {
			daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
			//dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
		    }
	       }
	       norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
	       dscal(size,1.0/norm,&*v_p.begin(),1);
		
	       //cout<<"vs.size() = "<<vs.size()<<endl;
	       norm=0.0;
	       if (orth_failed)
	       {
			while (abs(norm)<1.0e-2)
			{
				for (j=0;j<size;j++)  {v_p_old[j]=(uniform_rnd());}
				norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
				dscal(size,1.0/norm,&*v_p_old.begin(),1);
				//cout<<"norm="<<norm<<endl;
				v_p=v_p_old;
				for (i=0;i<vs.size();i++)
				{
					q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
					//cout<<"q (after orth fail )="<<q<<endl;
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
				}
				norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				//cout<<"norm="<<norm<<endl;
				dscal(size,1.0/norm,&*v_p.begin(),1);
				for (i=0;i<vs.size();i++)
				{
					q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
					cout<<"q (after orth fail )="<<q<<endl;
				}			
			}
			//dscal(size,1.0/norm,&*v_p.begin(),1);
	       }
       }

       time (&end);
       dif=difftime(end,start);
       
       if (ipr) cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
   }

   if (ipr) cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (j=0;j<iterations*iterations;j++) t_mat[j]=0.0;

   for (j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
  
   if (ipr) cout<<"Time to build T was "<<dif<<" seconds"<<endl;
 
   time (&start);
   symmetric_diagonalize(t_mat,eigs,t_eigenvecs);
   //cout<<"Tmat"<<endl;
   //print_real_mat(t_mat);
   time (&end);
   dif=difftime(end,start);

   if (ipr) cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;
	  
   how_many_eigenvecs=min(size,how_many_eigenvecs);
 
   eigenvecs.resize(size,how_many_eigenvecs);

   //if (ipr) cout<<"how_many_eigenvecs="<<how_many_eigenvecs<<endl;//if (ipr) cout<<"size="<<size<<endl;

   for (j=0;j<how_many_eigenvecs;j++)
   {
	      if (ipr) cout<<"Making "<<j<<"th Ritz eigenvector"<<endl;
	      for (i=0;i<size;i++)
	      {
		eigenvecs(i,j)=0.0;
		for (k=0;k<iterations;k++){eigenvecs(i,j)+=vs[k][i]*t_eigenvecs(k,j);}
	      }
   }      
   //f (ipr) print_real_mat(eigenvecs);
   if (ipr) cout<<"Done computing Ritz eigenvectors"<<endl;
   time(&end_0);
   dif=difftime(end_0,start_0);
   if (ipr) cout<<"Total time for all Lanczos iterations was "<<dif<<" seconds"<<endl;
}


