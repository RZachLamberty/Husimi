#include"lanczos_more.h"
#include"ed.h"
#include"printing_functions.h"
#include"hamiltonian_spin_functions.h"
#include"global.h"
#include"number_functions.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
void lanczos_given_map_evecs(Ham &h,
                      int iterations, 
		      std::vector<int>  const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &spins,
                      std::vector< std::vector<double> > &lowest_evecs,
                      int states_requested,
		      bool ipr)
{
   time_t 				start,end;
   int 					i,it,j,k,l,m;
   int 					size=map.size();
   iterations=min(iterations,size);
   
   double 				spin;
   double 				dif;
   double               		tmp;
   double 				q,alpha,beta,norm;
   complex<double> 			c_first_second,tot;
   bool 				orth_failed;
   std::vector<int> 			config(h.num_sites),new_configs;
   std::vector<double>  		alphas,betas;
   std::vector<double>  		h_dot_v(size),w(size);
   std::vector<double>			v_p(size),v_o(size);
   std::vector<double> 		        v_p_old(size);
   std::vector<double>			ritz_eigenvec(size);
   std::vector< complex<double> > 	hints_list;
   std::vector< std::vector<int> > 	touched_sites_list,vals_on_touched_list,vec_new_configs;
   std::vector< std::vector<double> >   vs;
   std::vector< std::vector< complex<double> > > vec_hints_list;
   Matrix 				t_mat(iterations,iterations);
   Matrix 				t_eigenvecs(iterations,iterations);
   bool 				measure_s=true;
  
   eigs.resize(iterations);
   spins.clear();
   lowest_evecs.clear();
   states_requested=min(states_requested,size);
    
   // Initializations
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
       
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;
   betas.push_back(beta);
       
   for (it=0;it<iterations;it++)
   {
       if (ipr) cout<<"Doing Lanczos iteration = "<<it<<endl;
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
	       for (i=0;i<it;i++)
	       {
		    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
		    if (abs(q)>1.0e-10)
		    {
			i=it;
			if(ipr) {cout<<"q (overlap) ="<<q<<endl;
				cout<<"--------------------------------------------------"<<endl;
				cout<<"Orthogonalization failed... choosing random vector"<<endl;
				cout<<"--------------------------------------------------"<<endl;
				orth_failed=true;
				}
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
				for (j=0;j<size;j++)  {v_p_old[j]=1.0+uniform_rnd();}
				norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
				dscal(size,1.0/norm,&*v_p_old.begin(),1);
				v_p=v_p_old;
				for (i=0;i<it;i++)
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

       time (&end);
       dif=difftime(end,start);
       
       //cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
   }

   if (ipr) cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
  
   if (ipr) {cout<<"============================================"<<endl; 
	     cout<<"Time to build  T was "<<dif<<" seconds"<<endl;
	     cout<<"============================================"<<endl; 
	     cout<<"Iterations="<<iterations<<endl; 
   	    }

   time (&start);
   symmetric_diagonalize(t_mat,eigs,t_eigenvecs);
   time (&end);
   dif=difftime(end,start);
  
   if (ipr)
   { 
	   cout<<endl;
	   cout<<"============================================"<<endl; 
	   cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;
	   cout<<"============================================"<<endl; 
	   cout<<endl;
   }

   if (measure_s)
   {
	   for (j=0;j<states_requested;j++)
	   {
	      if (ipr) cout<<"Making "<<j<<" th Ritz eigenvector"<<endl;

	      for (i=0;i<size;i++)
	      {
		ritz_eigenvec[i]=0.0;
		for (k=0;k<iterations;k++){ritz_eigenvec[i]+=vs[k][i]*t_eigenvecs(k,j);}
	      }
	      
	      if (ipr) cout<<"Finished computing Ritz eigenvector"<<endl;
	      if (ipr) cout<<"Recording lowest eigenvec"<<endl;
	      lowest_evecs.push_back(ritz_eigenvec);

	      //print_vec(ritz_eigenvec,true);
              
	      //cout<<"h.num_sites="<<h.num_sites<<endl; 

	      i=0;  	
	      //while (abs(ritz_eigenvec[i])<1.0e-2 and ritz_eigenvec.size()>1) {i+=1;}
	      while (abs(ritz_eigenvec[i])<1.0e-3) {i+=1;}

	      if (ipr) cout<<"Found non zero element i = "<<i<<endl;
	      //cout<<"map[i] = "<<map[i]<<endl; 
             
              //cout<<"h.num_sites="<<h.num_sites<<endl; 
              //cout<<"Converting num to vector"<<endl;
	      //convert_num_to_vec(1,2,h.num_sites,config);
	      //print_vec(config);

	      //convert_num_to_vec(map[i],2,h.num_sites,config);
	      tot=0.0;
	     
              //cout<<"Converted num to vector"<<endl;
              //if (ritz_eigenvec.size()>1)
              //{
	      	      convert_num_to_vec(map[i],2,h.num_sites,config);
		
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
				//cout<<"FINISHED hint calculaion "<<endl;
				c_first_second+=(hints_list[0]*ritz_eigenvec[i]);
				tot+=c_first_second;	  
			}
		      }		
		     
		      //cout<<"tot="<<tot<<endl;
		      tot=tot/ritz_eigenvec[i]; 
		      tmp=(2.0*real(tot))+(0.75*double(h.num_sites));
		      //cout<<"S(S+1) = "<<tmp<<endl;
		      spin=compute_spin(tmp);
		      spins.push_back(spin);
		      cout<<endl;
		      if (ipr)
		      {cout<<"======================================"<<endl;
		      cout<<"Spin of state ="<<spin<<endl;
		      cout<<"S_z  of state ="<<measure_sz(config)<<endl;
		      cout<<"E    of state ="<<eigs[j]<<endl;	
		      cout<<"======================================"<<endl;
		      cout<<endl;
		      }
	     //}
	     /*else
	     {
		cout<<"Got here"<<endl;
	        spins.clear();
		eigs.clear();
		spins.push_back(0);
		eigs.push_back(1.0e10);
		cout<<"Finished recording"<<endl;
	     }*/
    	}
    }

    //cout<<"Map at end"<<endl;
    //for (int n=0;n<map.size();n++) cout<<map[n]<<endl;
}


