#include"nrg_spin.h"
#include"xxz_model.h"
#include"ed.h"
#include<math.h>

using namespace std;


double splusfn(double spin,int ind)
{
	double sz=double(ind)-spin;
	//return 1;
	return sqrt(abs((spin-sz)*(spin+sz+1)));
}

double splus_sminus_fn(double spin,int ind)
{
	double sz=double(ind)-spin;
        //return 1;
	return sqrt(abs((spin+sz)*(spin-sz+1)*(spin-sz+1)*(spin+sz)));
}

double sminus_splus_fn(double spin,int ind)
{
	double sz=double(ind)-spin;
        //return 1;
	return sqrt(abs((spin-sz)*(spin+sz+1)*(spin-sz)*(spin+sz+1)));
}


double sminusfn(double spin,int ind)
{
	double sz=double(ind)-spin;
        //return 1;
	return sqrt(abs((spin+sz)*(spin-sz+1)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void nrg_spin(Spin_Half_XXZ 					xxz,
              std::vector<double>                               target_szs,
              std::vector<int>                                  target_numbers,
	      std::vector< std::vector<double> > 		&eigs,
              std::vector<double>				&spins,
    	      std::vector< std::vector<int> >                   &inverse_subspace_map,
    	      std::vector< std::vector<int> >                   &inverse_sz_map,
	      std::vector< std::vector< Matrix > >  	        &transforms,
              std::vector< std::vector< Matrix > >              &s_plus_block,
              std::vector< std::vector< Matrix > >              &s_minus_block,
              std::vector< std::vector< Matrix > > 	        &s_z_block,
              std::vector< std::vector< Matrix > > 	        &s_z2_block,
	      std::vector< std::vector< Matrix > > 	        &ham_block,
	      int 						max_states, 
	      bool 					        ipr)
{
    int                                i,j,size;
    int 			       retained;
    std::vector<int>                   nstates;
    std::vector<double> 	       szs,energies,tmpszs;
    Matrix                             single_ham(1,1),eye(1,1),null(0,0),mat(1,1);	
    std::vector< std::vector<double> > szs_on_block;
    double                             sz;
    bool                               measure=true;
    double 			       spin=xxz.spin;
    double 			       gen=xxz.gen;
    double 			       j_x=xxz.j_x;
    double 			       j_z=xxz.j_z;
    double 			       d=xxz.d;
    Hier 			       hier;

    hier.init(gen);
    cout<<"Spin on each site  ="<<spin<<endl;    
    cout<<"Max Gen of cactus  ="<<gen<<endl;    
    max_states=min(60,max_states);
    cout<<"max states (inside RG routine)"<<max_states<<endl;
    // Clear everything that enters this routine
    eigs.clear(); spins.clear();transforms.clear();ham_block.clear();
    s_plus_block.clear();s_minus_block.clear();s_z_block.clear();s_z2_block.clear();
    inverse_subspace_map.clear(); inverse_sz_map.clear();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    eye(0,0)=1.0;single_ham(0,0)=0.0;

    cout<<"Preparing..............."<<endl;
    for (i=0;i<gen+1;i++)
    {
	s_plus_block.push_back(std::vector<Matrix>());
	s_minus_block.push_back(std::vector<Matrix>());
	s_z_block.push_back(std::vector<Matrix>());
	s_z2_block.push_back(std::vector<Matrix>());
        ham_block.push_back(std::vector<Matrix>());
        transforms.push_back(std::vector<Matrix>());
        szs_on_block.push_back(std::vector<double>());
        inverse_subspace_map.push_back(std::vector<int>());
        inverse_sz_map.push_back(std::vector<int>());
    }		
    cout<<"Finished Preparing 1....."<<endl;


    // Initial spins of spin length S - i.e. boundary spins which comprise generation 0 - generalize to spin>=1/2 which has in general 2*S+1 states 

    for (i=0;i<int(2.0*spin+1.0e-6);i++) 
    {
	//sz=double(i)-spin;
	eye(0,0)=splusfn(spin,i);
	s_plus_block[0].push_back(eye);
    }
    s_plus_block[0].push_back(null);  // Highest s_z state can not be raised

    s_minus_block[0].push_back(null); // Lowest s_z state can not be lowered
    for (i=1;i<int(2*spin+1.0e-6)+1;i++) 
    {
	//sz=double(i)-spin;
	eye(0,0)=sminusfn(spin,i);
	s_minus_block[0].push_back(eye);
    }

    for (i=0;i<int(2*spin+1.0e-6)+1;i++)    // All S_z states from -S to +S S=spin of site
    { 
	mat(0,0)=double(i)-spin;
        s_z_block[0].push_back(mat);
	mat(0,0)=pow(double(i)-spin,2.0);
        s_z2_block[0].push_back(mat);
        szs_on_block[0].push_back(double(i)-spin); 
        inverse_sz_map[0].push_back(0);
        inverse_subspace_map[0].push_back(i);
    }
    
    for (i=0;i<int(2*spin+1.0e-6)+1;i++) 
    {
	single_ham(0,0)=d*pow((double(i)-spin),2.0);
	ham_block[0].push_back(single_ham);  // Hamiltonian in various S_z sectors is just 0 for free spins, changes when you have anisotropy i.e. d is non zero
    }
    cout<<"Finished Preparing 2 ....."<<endl;
    hier.s_plus_hier[0].push_back(s_plus_block[0]);
    hier.s_minus_hier[0].push_back(s_minus_block[0]);
    hier.s_z_hier[0].push_back(s_z_block[0]);
    hier.s_z2_hier[0].push_back(s_z2_block[0]);

    for (int g=0;g<gen;g++)
    {
	    nstates.clear();
            retained=0;
            for (i=0;i<s_z_block[g].size();i++) retained=retained+s_z_block[g][i].NRows();
	    for (i=0;i<4;i++) nstates.push_back(retained);
            nstates.push_back(int(2*spin+1.0e-6)+1);
	    /*  PRINT STATEMENT
	    cout<<endl;
            cout<<"Cactus at generation "<<g<<endl; 		
	    cout<<"Szs at present generation "<<endl; 		
            print_vec(szs_on_block[g]);
	    cout<<"Hams at present generation "<<endl; 		
            for (i=0;i<ham_block[g].size();i++) print_real_mat(ham_block[g][i]);
        */
	    nrg_dmrg_spin_cactus(xxz,g,target_szs,target_numbers,nstates,max_states,
					 szs_on_block[g],inverse_subspace_map[g],inverse_sz_map[g],s_plus_block[g],s_minus_block[g],s_z_block[g],s_z2_block[g],ham_block[g],
				         szs_on_block[g+1],inverse_subspace_map[g+1],inverse_sz_map[g+1],s_plus_block[g+1],s_minus_block[g+1],s_z_block[g+1],s_z2_block[g+1],ham_block[g+1],
                                         transforms[g+1],hier,4*max_states,measure,ipr);
	    //Szs are compact - now expand them!!!!
            tmpszs.clear();
	    for (i=0;i<s_z_block[g+1].size();i++)
	    { 
		for (j=0;j<s_z_block[g+1][i].NRows();j++)
		{
			tmpszs.push_back(szs_on_block[g+1][i]);
		}
            }
            /*  PRINT STATEMENT
            cout<<"Szs at generation + 1  "<<endl; 		
            print_vec(szs_on_block[g+1]);
	        cout<<endl;
	        */
            szs_on_block[g+1].clear();
            szs_on_block[g+1]=tmpszs;
            /*  PRINT STATEMENT
            cout<<"Szs at generation + 1 (AFTER EXPANSION) "<<endl; 		
            print_vec(szs_on_block[g+1]);
            */
            //cout<<"Hams at generation +1  "<<endl; 		
            //for (i=0;i<ham_block[g+1].size();i++) print_real_mat(ham_block[g+1][i]);
            //cout<<"Transforms at generation +1  "<<endl; 		
            //for (i=0;i<transforms[g+1].size();i++)print_real_mat(transforms[g+1][i]);
    }		
    /*if (ipr) {cout<<endl;}
    cout<<"TRACE: End result of i-DMRG (pre sweep DMRG) calculation"<<endl;
    for (int s=0;s<szs.size();s++)
    {
	    if (energies[s].size()>0) {cout<<"Energies for Sz = "<<szs[s]<<endl;
	    			       print_vec_acc(energies[s],true);	
	   			       cout<<endl;}
    }
    
    for (int s=0;s<szs.size();s++)
    {
	if (energies[s].size()>0 and (szs[s]>0 or abs(szs[s])<1.0e-10)) 
   	{
	  for (int i=0;i<energies[s].size();i++)
   	  {cout<<boost::format("%+.10f") %energies[s][i]<<"  "<<boost::format("%+i") %closest_int(szs[s])<<endl;}
	}
    }
    cout<<endl;cout<<endl;
    eigs=energies;*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void nrg_dmrg_spin_cactus(Spin_Half_XXZ				xxz,
			  int                                   gen,
			  std::vector<double> 			input_target_szs,
			  std::vector<int> 			input_target_states,
			  std::vector<int> 			block_num_states,
			  int 				        max_trunc_states, 
			  std::vector<double>                   &szs_on_block,
    	 		  std::vector<int>                      &inverse_subspace_map,
    	 		  std::vector<int>                      &inverse_sz_map,
			  std::vector< Matrix >                 &s_plus_block,
			  std::vector< Matrix >                 &s_minus_block,
			  std::vector< Matrix >                 &s_z_block,
			  std::vector< Matrix >                 &s_z2_block,
			  std::vector< Matrix >                 &local_hams,
			  std::vector<double>                   &szs,
    	 		  std::vector<int>                      &inverse_subspace_map_out,
    	 		  std::vector<int>                      &inverse_sz_map_out,
			  std::vector<Matrix> 		        &s_p_psite, 
			  std::vector<Matrix> 		        &s_m_psite, 
			  std::vector<Matrix> 		        &s_z_psite, 
			  std::vector<Matrix> 		        &s_z2_psite, 
			  std::vector<Matrix>                   &ham_out,
			  std::vector<Matrix> 		        &opt_evecs, 
			  Hier                                  &hier,
			  int 				        lanc_dav_it,
			  bool 					measure,
			  bool					ipr)
{
	double j_x=xxz.j_x;
        double j_z=xxz.j_z;
	double d=xxz.d;
        double spin=xxz.spin;
	ipr=true;
        if (ipr) {cout<<boost::format("J_z,J_x,D = %+.10f %+.10f %.10f") %j_z %j_x %d<<endl;} 
    ipr=false;
	// Clear unnecessary info
	ham_out.clear();
        inverse_subspace_map_out.clear();inverse_sz_map_out.clear();
        s_p_psite.clear();s_m_psite.clear();s_z_psite.clear();s_z2_psite.clear();

	// Variables
	time_t								start,start_0,end,end_0;
	bool 								flg,spin_all_found,spin_sys_found;
	bool 								env_changed,target_states_not_found;
        bool								bfound;
	int 								num,final_num;
	int 								subspace_all,subspace_sys;
	int 								env_state,prev_env_state;
	int 								ind_1,ind_2,i1,c1,c2;
	int 								n_all=1,n_sys=1,n_env1,n_env2;
	int 								num_spaces_all=0,num_spaces_sys=0;
	int 								num_blocks,total_states;
        int								state_location,subspace_location,subspace1_location,subspace2_location,ind_subspace_location,ind_subspace1_location,ind_subspace2_location;	
        int								state1_location,state1_new_location,state2_location,state2_new_location,state_1_location,state_2_location;	
	int 								plus_one_subsp,minus_one_subsp,present_subspace;
	int								shift,m_1,m_2,m_1_inv,m_2_inv,highest,selected;
	int 								sz_ind;
        int 								rind_2;
        double 								hint,cutoff,dif;
	double								s_z_total,s_z_system,s_z_env1,s_z_env2,s_z_rdm;
        std::vector<double>						tmp_en;
	std::vector<int> 						vec_1,vec_2;
        std::vector< std::vector<int> >					list_of_converted_vecs;	
	std::vector<int>                                                inverse_sz_map_all,inverse_subspace_map_all;
	std::vector< std::vector<int> >                                 map_for_states_sys;
	std::vector<int>                                                inverse_subspace_map_sys,inverse_sz_map_sys;
        std::vector<int>                                                states_sys;
        int                                                             block_plus_subspace,block_minus_subspace; 
        for (int i=0;i<block_num_states.size();i++){n_all=n_all*block_num_states[i];}
	for (int i=2;i<block_num_states.size();i++){n_sys=n_sys*block_num_states[i];states_sys.push_back(block_num_states[i]);}
	if (ipr) cout<<"TRACE: Total number of states env1+env2+sys1+sys2+site ="<<n_all<<endl;
	if (ipr) cout<<"TRACE: Total number of states sys1+sys2+site ="<<n_sys<<endl;

	map_for_states_sys.clear();
        inverse_subspace_map_sys.clear();inverse_subspace_map_sys.resize(n_sys);
        inverse_sz_map_sys.clear();inverse_sz_map_sys.resize(n_sys);
	inverse_sz_map_all.clear();inverse_sz_map_all.resize(n_all);
	inverse_subspace_map_all.clear();inverse_subspace_map_all.resize(n_all);

	std::vector< std::vector<double> >				eigs_sys,eigs_all;
	std::vector< std::vector<double> >				rdm_eigs;
	std::vector< std::vector<int> > 				map_for_states_all;
	std::vector< std::vector< std::vector<int> > > 			map_for_hints_all,map_for_hints_sys;
        std::vector<int>                                                target_states;
	std::vector<double>						target_state;
	std::vector<int>   				                shifts;
	std::vector< std::vector< std::vector<double> > >		hints_all,hints_sys;
	Matrix 								eigenvecs;
	std::vector<Matrix>						opt_evecs_2_all,opt_evecs_2_sys;
	std::vector<Matrix>						rdm_evecs,opt_evecs_rdm;
	std::vector< std::vector<double> > 				empty_dbl;
	std::vector< std::vector<int> > 				empty_int;
        std::vector<int>  				                empty_int_vec;
	std::vector<double>                                             szs_all,newsz_all,newsz_sys;
	std::vector<double>						tmp_eigs,tmp_dmat_eigs;
	std::vector<Matrix>						rdms,avg_rdms;
        std::vector<Matrix>                                             htransforms;
	Matrix								target_state_rdm;	
        double								h=0.01; 
        
	//We consider each subspace of S_z = fixed and diagonalize in that subspace
	//====================================================================================================================
	//                                            Diagonal terms and creation of maps
	//====================================================================================================================
	if (ipr) cout<<"Computing diagonal terms "<<endl;
	shift=0;
	for (int j=0;j<s_z_block.size();j++)
	{
		shifts.push_back(shift);
		//cout<<"j (shift)="<<j<<endl;
		shift=shift+s_z_block[j].NRows();
	}
	/*  PRINT STATEMENT
	print_vec(shifts);
	*/

	for (ind_1=0;ind_1<n_all;ind_1++)
	{
		vec_1=convert_ind_to_vec(ind_1,block_num_states);
		list_of_converted_vecs.push_back(vec_1);

		env_state=vec_1[0]+vec_1[1];
 	        //cout<<"env_state="<<env_state<<endl;
	
		s_z_total=0.0;
		for (int i=0;i<vec_1.size()-1;i++) s_z_total=s_z_total+szs_on_block[vec_1[i]];

                s_z_total+=(double(vec_1[vec_1.size()-1])-spin);
		s_z_system=s_z_total-szs_on_block[vec_1[0]]-szs_on_block[vec_1[1]]; // Subtract off spins of 2 environments
		//if (ipr) cout<<"S_z_total on system + env ="<<s_z_total<<endl;if (ipr) cout<<"S_z_total on system alone ="<<s_z_system<<endl;

		spin_all_found=false; spin_sys_found=false;

		for (int j=0;j<map_for_states_all.size();j++)
		{
			if (abs(newsz_all[j]-s_z_total)<1.0e-10) 
			{
				subspace_all=j;
				spin_all_found=true;
				j=map_for_states_all.size();
				// If found, add a new entry to existing subspace
			}
		}

		if (env_state==0)
		{
			//cout<<"Got here (env=0)"<<endl;
			for (int j=0;j<map_for_states_sys.size();j++)
			{
				if (abs(newsz_sys[j]-s_z_system)<1.0e-10) 
				{
					subspace_sys=j;
					spin_sys_found=true;
					j=map_for_states_sys.size();
					// If found, add a new entry to existing subspace
				}
			}
			//if (ipr) cout<<"Spin for system found = "<<spin_sys_found<<endl;
			if (not spin_sys_found) // Add a new entry for maps and hints
			{
				map_for_states_sys.push_back(empty_int_vec);
				map_for_hints_sys.push_back(empty_int);
				hints_sys.push_back(empty_dbl);
				newsz_sys.push_back(s_z_system);
				szs.push_back(s_z_system);
				subspace_sys=map_for_states_sys.size()-1;
				num_spaces_sys=num_spaces_sys+1;
				opt_evecs_2_sys.push_back(eigenvecs);
				opt_evecs.push_back(eigenvecs);
				rdm_evecs.push_back(eigenvecs);
				eigs_sys.push_back(std::vector<double> ());
				s_p_psite.push_back(eigenvecs);
				s_m_psite.push_back(eigenvecs);
				s_z_psite.push_back(eigenvecs);
				s_z2_psite.push_back(eigenvecs);
			}
			//if (ipr) cout<<"Made additions for new spin found"<<endl;
			map_for_states_sys[subspace_sys].push_back(ind_1);// Dont need to change ind_1 since env_state=0!
			inverse_sz_map_sys[ind_1]=map_for_states_sys[subspace_sys].size()-1;
			inverse_subspace_map_sys[ind_1]=subspace_sys;
			map_for_hints_sys[subspace_sys].push_back(std::vector<int> ());
			hints_sys[subspace_sys].push_back(std::vector<double> ());
			map_for_hints_sys[subspace_sys][map_for_states_sys[subspace_sys].size()-1].push_back(map_for_states_sys[subspace_sys].size()-1);
			hints_sys[subspace_sys][map_for_hints_sys[subspace_sys].size()-1].push_back(0.0);
			//if (ipr) cout<<"Maps updated after new spin found"<<endl;
			
			for (int i=2;i<vec_1.size()-1;i++)                         // the first 2 indices correspond to environments
			{
				state_location=inverse_sz_map[vec_1[i]];           //cout<<"State Location (diag) = "<<state_location<<endl;
				subspace_location=inverse_subspace_map[vec_1[i]];  //cout<<"Subspace Location (diag) = "<<subspace_location<<endl;
				if (local_hams[subspace_location].NRows()>0)
				{hints_sys[subspace_sys][map_for_hints_sys[subspace_sys].size()-1][0]+=local_hams[subspace_location](state_location,state_location);
				 hints_sys[subspace_sys][map_for_hints_sys[subspace_sys].size()-1][0]+=(j_z*(double(vec_1[vec_1.size()-1])-spin)*s_z_block[subspace_location](state_location,state_location));
				}
			}
			hints_sys[subspace_sys][map_for_hints_sys[subspace_sys].size()-1][0]+=(d*pow(double(vec_1[vec_1.size()-1])-spin,2.0)); // Magnetic field coupling to S_z^2 of central site
                        state1_location=inverse_sz_map[vec_1[2]];state2_location=inverse_sz_map[vec_1[3]];
                        subspace1_location=inverse_subspace_map[vec_1[2]];subspace2_location=inverse_subspace_map[vec_1[3]];
			hints_sys[subspace_sys][map_for_hints_sys[subspace_sys].size()-1][0]+=(j_z*s_z_block[subspace1_location](state1_location,state1_location)*s_z_block[subspace2_location](state2_location,state2_location));
		}
		if (not spin_all_found) // Add a new entry for maps and hints
		{
			map_for_states_all.push_back(empty_int_vec);
			map_for_hints_all.push_back(empty_int);
			hints_all.push_back(empty_dbl);
			newsz_all.push_back(s_z_total);
			szs_all.push_back(s_z_total);
			subspace_all=map_for_states_all.size()-1;
			num_spaces_all=num_spaces_all+1;
			opt_evecs_2_all.push_back(eigenvecs);
			eigs_all.push_back(std::vector<double> ());
			target_states_not_found=true;
			for (int i=0;i<input_target_states.size();i++)
			{
				if (abs(input_target_szs[i]-s_z_total)<1.0e-10) {target_states.push_back(input_target_states[i]);target_states_not_found=false;}
				
			}
			if (target_states_not_found) {target_states.push_back(0);}
		}
		//cout<<"Subspace ="<<subspace<<endl;
		map_for_states_all[subspace_all].push_back(ind_1);
		inverse_sz_map_all[ind_1]=map_for_states_all[subspace_all].size()-1;
		inverse_subspace_map_all[ind_1]=subspace_all;
		
		map_for_hints_all[subspace_all].push_back(std::vector<int> ());
		hints_all[subspace_all].push_back(std::vector<double> ());

		map_for_hints_all[subspace_all][map_for_states_all[subspace_all].size()-1].push_back(map_for_states_all[subspace_all].size()-1);
		hints_all[subspace_all][map_for_hints_all[subspace_all].size()-1].push_back(0.0);
		
		for (int i=0;i<vec_1.size()-1;i++)
		{
			state_location=inverse_sz_map[vec_1[i]];
			//cout<<"State Location (diag) = "<<state_location<<endl;
                        subspace_location=inverse_subspace_map[vec_1[i]];
			//cout<<"Subspace Location (diag) = "<<subspace_location<<endl;
			if (local_hams[subspace_location].NRows()>0)
			{hints_all[subspace_all][map_for_hints_all[subspace_all].size()-1][0]+=local_hams[subspace_location](state_location,state_location);
			hints_all[subspace_all][map_for_hints_all[subspace_all].size()-1][0]+=(j_z*(double(vec_1[vec_1.size()-1])-spin)*s_z_block[subspace_location](state_location,state_location));
			}
		}
		//if (ipr) cout<<"Finished diagonal calculations for state with ind_1 = "<<ind_1<<endl;
		hints_all[subspace_all][map_for_hints_all[subspace_all].size()-1][0]+=(d*pow(double(vec_1[vec_1.size()-1])-spin,2.0)); // Magnetic field coupling to S_z^2 of central site
                state1_location=inverse_sz_map[vec_1[2]];state2_location=inverse_sz_map[vec_1[3]];
                subspace1_location=inverse_subspace_map[vec_1[2]];subspace2_location=inverse_subspace_map[vec_1[3]];
		hints_all[subspace_all][map_for_hints_all[subspace_all].size()-1][0]+=(j_z*s_z_block[subspace1_location](state1_location,state1_location)*s_z_block[subspace2_location](state2_location,state2_location));
                state1_location=inverse_sz_map[vec_1[0]];state2_location=inverse_sz_map[vec_1[1]];
                subspace1_location=inverse_subspace_map[vec_1[0]];subspace2_location=inverse_subspace_map[vec_1[1]];
		hints_all[subspace_all][map_for_hints_all[subspace_all].size()-1][0]+=(j_z*s_z_block[subspace1_location](state1_location,state1_location)*s_z_block[subspace2_location](state2_location,state2_location));
	}

	//if (ipr) cout<<"Inv sz map of sys"<<endl; if (ipr) print_vec(inverse_sz_map_sys);
	//if (ipr) cout<<"Szs found"<<endl; if (ipr) print_vec(szs_all);
        //for (i=0;i<hints_sys.size();i++) print_mat_double(hints_sys[s]);
	//===================================================================================================================================================
	//                                                  Off diagonal terms for sys - P1 (number of loops can be reduced!)
	//===================================================================================================================================================
	for (int subspace=0;subspace<map_for_states_sys.size();subspace++)
	{
		if (ipr) cout<<"Off diagonal terms of system alone (Part 1) being evaluated for subspace = "<<subspace<<" which has "<<map_for_states_sys[subspace].size()<<" states"<<endl;
		for (int m=0;m<map_for_states_sys[subspace].size();m++)
		{
			ind_1=map_for_states_sys[subspace][m];
			vec_1=list_of_converted_vecs[ind_1];
			//cout<<"ind_1="<<ind_1<<endl;//print_vec(vec_1);
			for (i1=2;i1<vec_1.size()-1;i1++)                                                           // Skip the environments
			{
				subspace_location=inverse_subspace_map[vec_1[i1]];
				state_1_location=inverse_sz_map[vec_1[i1]];
				
				vec_2=vec_1;
                                sz_ind=vec_1[vec_1.size()-1];
				if (vec_1[vec_1.size()-1]!=int(2*spin+1.0e-6))
				{
					vec_2[vec_1.size()-1]=vec_1[vec_1.size()-1]+1;
					ind_subspace_location=subspace_location-1;
					if (ind_subspace_location>=0)
					{
						for (int j=0;j<s_z_block[ind_subspace_location].NRows();j++)
						{
							vec_2[i1]=j+shifts[ind_subspace_location];  //cout<<"vec_2[i1] = "<<vec_2[i1]<<endl;//cout<<"shift[i1][ind_subspace_location = "<<shifts[i1][ind_subspace_location]<<endl;
							state_2_location=inverse_sz_map[vec_2[i1]]; // must be j; //cout<<"j="<<j<<" state_2_location="<<state_2_location<<endl;
							ind_2=convert_vec_to_ind(vec_2,block_num_states);
							hint=0.0;
							if (s_minus_block[subspace_location].NRows()>=1) hint+=(j_x*(0.5*splusfn(spin,sz_ind)*s_minus_block[subspace_location](state_2_location,state_1_location)));
							if ((ind_1!=ind_2) and (hint!=0.0)) 
							{
								map_for_hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(inverse_sz_map_sys[ind_2]);
								hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(hint);
							}
						}
					}
				}
	
				vec_2=vec_1;
				if (vec_1[vec_1.size()-1]!=0)
				{
					vec_2[vec_1.size()-1]=vec_1[vec_1.size()-1]-1;
					ind_subspace_location=subspace_location+1;
					if (ind_subspace_location<s_z_block.size())
					{
						for (int j=0;j<s_z_block[ind_subspace_location].NRows();j++)
						{
							vec_2[i1]=j+shifts[ind_subspace_location];  //cout<<"vec_2[i1] = "<<vec_2[i1]<<endl;
							state_2_location=inverse_sz_map[vec_2[i1]]; // must be j; //cout<<"j="<<j<<" state_2_location="<<state_2_location<<endl;
							ind_2=convert_vec_to_ind(vec_2,block_num_states);
							hint=0.0;
							if (s_plus_block[subspace_location].NRows()>=1) hint+=(j_x*(0.5*sminusfn(spin,sz_ind)*s_plus_block[subspace_location](state_2_location,state_1_location)));
							if ((ind_1!=ind_2) and (hint!=0.0)) 
							{
								map_for_hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(inverse_sz_map_sys[ind_2]);
								hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(hint);
							}
						 }
					 }
				 }
					
				 vec_2=vec_1;	
				 ind_subspace_location=subspace_location;
				 for (int j=0;j<s_z_block[ind_subspace_location].NRows();j++)
				 {
						vec_2[i1]=j+shifts[ind_subspace_location];
						state_2_location=inverse_sz_map[vec_2[i1]]; // must be j;
						ind_2=convert_vec_to_ind(vec_2,block_num_states);
						hint=0.0;
						if (s_z_block[subspace_location].NCols()>=1)
						{
							hint+=(j_z*(double(vec_1[vec_1.size()-1])-spin)*s_z_block[subspace_location](state_1_location,state_2_location));
							hint+=local_hams[subspace_location](state_1_location,state_2_location);
						}
						if ((ind_1!=ind_2) and (hint!=0.0)) 
						{
							map_for_hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(inverse_sz_map_sys[ind_2]);
							hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(hint);
						}
				  }
			}		
		}
	}

	//===================================================================================================================================================
	//                                                             Off diagonal terms for sys - P2
	//===================================================================================================================================================
	
        for (int subspace=0;subspace<map_for_states_sys.size();subspace++)
	{ //sub
		if (ipr) cout<<"Off diagonal terms of system alone (Part 2) being evaluated for subspace = "<<subspace<<" which has "<<map_for_states_sys[subspace].size()<<" states"<<endl;
		for (int m=0;m<map_for_states_sys[subspace].size();m++)
		{ //m
			ind_1=map_for_states_sys[subspace][m];
			vec_1=list_of_converted_vecs[ind_1];
			subspace1_location=inverse_subspace_map[vec_1[2]];
			state1_location=inverse_sz_map[vec_1[2]];
			subspace2_location=inverse_subspace_map[vec_1[3]];
			state2_location=inverse_sz_map[vec_1[3]];
			vec_2=vec_1;
			ind_subspace1_location=subspace1_location-1;
			ind_subspace2_location=subspace2_location+1;
			if (ind_subspace1_location>=0 and ind_subspace2_location<s_z_block.size())
			{
				for (int j=0;j<s_z_block[ind_subspace1_location].NRows();j++)
				{
					vec_2[2]=j+shifts[ind_subspace1_location];
					state1_new_location=inverse_sz_map[vec_2[2]]; // must be j;
					for (int k=0;k<s_z_block[ind_subspace2_location].NRows();k++)
					{
						vec_2[3]=k+shifts[ind_subspace2_location];
						state2_new_location=inverse_sz_map[vec_2[3]]; // must be k;
						ind_2=convert_vec_to_ind(vec_2,block_num_states);
						hint=0.0;
						if (s_minus_block[subspace1_location].NRows()>=1 and s_plus_block[subspace2_location].NRows()>=1) 
						{hint+=(j_x*(0.5*s_plus_block[subspace2_location](state2_new_location,state2_location)*s_minus_block[subspace1_location](state1_new_location,state1_location)));}
						if ((ind_1!=ind_2) and (hint!=0.0)) 
						{
							map_for_hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(inverse_sz_map_sys[ind_2]);
							hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(hint);
						}
					 }
				 }
			}
			vec_2=vec_1;
			ind_subspace1_location=subspace1_location+1;
			ind_subspace2_location=subspace2_location-1;
			if (ind_subspace2_location>=0 and ind_subspace1_location<s_z_block.size())
			{
				for (int j=0;j<s_z_block[ind_subspace1_location].NRows();j++)
				{
					vec_2[2]=j+shifts[ind_subspace1_location];
					state1_new_location=inverse_sz_map[vec_2[2]]; // must be j;
					for (int k=0;k<s_z_block[ind_subspace2_location].NRows();k++)
					{
						vec_2[3]=k+shifts[ind_subspace2_location];
						state2_new_location=inverse_sz_map[vec_2[3]]; // must be k;
						ind_2=convert_vec_to_ind(vec_2,block_num_states);
						hint=0.0;
						if (s_plus_block[subspace1_location].NRows()>=1 and s_minus_block[subspace2_location].NRows()>=1) 
						{hint+=(j_x*(0.5*s_minus_block[subspace2_location](state2_new_location,state2_location)*s_plus_block[subspace1_location](state1_new_location,state1_location)));}
						if ((ind_1!=ind_2) and (hint!=0.0)) 
						{
							map_for_hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(inverse_sz_map_sys[ind_2]);
							hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(hint);
						}
					}
				}
			}
			vec_2=vec_1;
			ind_subspace1_location=subspace1_location;
			ind_subspace2_location=subspace2_location;
			for (int j=0;j<s_z_block[ind_subspace1_location].NRows();j++)
			{
				vec_2[2]=j+shifts[ind_subspace1_location];
				state1_new_location=inverse_sz_map[vec_2[2]]; // must be j;
				for (int k=0;k<s_z_block[ind_subspace2_location].NRows();k++)
				{
					vec_2[3]=k+shifts[ind_subspace2_location];
					state2_new_location=inverse_sz_map[vec_2[3]]; // must be k;
					ind_2=convert_vec_to_ind(vec_2,block_num_states);
					hint=0.0;
					hint+=(j_z*(s_z_block[subspace2_location](state2_new_location,state2_location)*s_z_block[subspace1_location](state1_new_location,state1_location)));
					if ((ind_1!=ind_2) and (hint!=0.0)) 
					{
						map_for_hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(inverse_sz_map_sys[ind_2]);
						hints_sys[subspace][inverse_sz_map_sys[ind_1]].push_back(hint);
					}
				}
			}
					
		}
	} // sub

	//===================================================================================================================================================
	//                                                             Off diagonal terms for sys + environment
	//===================================================================================================================================================
	for (int subspace=0;subspace<map_for_states_all.size();subspace++)
	{
		time(&start_0);
		if (target_states[subspace]>0) // if target
		{
		if (ipr) {cout<<"Off diagonal terms of system+ environment (Part 1) being evaluated for subspace = "<<subspace<<" which has "<<map_for_states_all[subspace].size()<<" states"<<endl;}
		
		for (int m=0;m<map_for_states_all[subspace].size();m++)  // m loop
		{
			ind_1=map_for_states_all[subspace][m];
			vec_1=list_of_converted_vecs[ind_1];
			//cout<<"ind_1="<<ind_1<<endl;print_vec(vec_1);
			for (i1=0;i1<vec_1.size()-1;i1++)               // i1 loop
			{
				subspace_location=inverse_subspace_map[vec_1[i1]];
				state_1_location=inverse_sz_map[vec_1[i1]];
				
                                sz_ind=vec_1[vec_1.size()-1];
				vec_2=vec_1;
				if (vec_1[vec_1.size()-1]!=int(2*spin+1.0e-6))
				{
					vec_2[vec_1.size()-1]=vec_1[vec_1.size()-1]+1;
					ind_subspace_location=subspace_location-1;
					if (ind_subspace_location>=0)
					{
						for (int j=0;j<s_z_block[ind_subspace_location].NRows();j++)
						{
							vec_2[i1]=j+shifts[ind_subspace_location];
							state_2_location=inverse_sz_map[vec_2[i1]]; // must be j;
							ind_2=convert_vec_to_ind(vec_2,block_num_states);
							hint=0.0;
							if (s_minus_block[subspace_location].NRows()>=1) hint+=(j_x*(0.5*splusfn(spin,sz_ind)*s_minus_block[subspace_location](state_2_location,state_1_location)));
							if ((ind_1!=ind_2) and (hint!=0.0)) 
							{
								map_for_hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(inverse_sz_map_all[ind_2]);
								hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(hint);
							}
						}
					}
				}
	
				vec_2=vec_1;
				if (vec_1[vec_1.size()-1]!=0)
				{
					vec_2[vec_1.size()-1]=vec_1[vec_1.size()-1]-1;
					ind_subspace_location=subspace_location+1;

					if (ind_subspace_location<s_z_block.size())
					{
						for (int j=0;j<s_z_block[ind_subspace_location].NRows();j++)
						{
							vec_2[i1]=j+shifts[ind_subspace_location];
							state_2_location=inverse_sz_map[vec_2[i1]]; // must be j;
							ind_2=convert_vec_to_ind(vec_2,block_num_states);
							hint=0.0;
							if (s_plus_block[subspace_location].NRows()>=1) hint+=(j_x*(0.5*sminusfn(spin,sz_ind)*s_plus_block[subspace_location](state_2_location,state_1_location)));
							if ((ind_1!=ind_2) and (hint!=0.0)) 
							{
								map_for_hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(inverse_sz_map_all[ind_2]);
								hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(hint);
							}
						 }
					 }
				 }
					
				 vec_2=vec_1;	
				 ind_subspace_location=subspace_location;
				 for (int j=0;j<s_z_block[ind_subspace_location].NRows();j++)
				 {
						vec_2[i1]=j+shifts[ind_subspace_location];
						state_2_location=inverse_sz_map[vec_2[i1]]; // must be j;
						ind_2=convert_vec_to_ind(vec_2,block_num_states);
						hint=0.0;
						if (s_z_block[subspace_location].NCols()>=1)
						{
							hint+=(j_z*(double(vec_1[vec_1.size()-1])-spin)*s_z_block[subspace_location](state_1_location,state_2_location));
							hint+=local_hams[subspace_location](state_1_location,state_2_location);
						}
						if ((ind_1!=ind_2) and (hint!=0.0)) 
						{
							map_for_hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(inverse_sz_map_all[ind_2]);
							hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(hint);
						}
				  }
			         			
			}// i1 loop
		     }// m loop
		} //if target
        
		for (c1=0;c1<3;c1+=2)  // c1 loop
		{	
			c2=c1+1;
			//===================================================================================================================================================
			//                                                             Off diagonal terms for sys + env - P2
			//===================================================================================================================================================
			if (target_states[subspace]>0) // if target 
			{
				if (ipr and c1==0) {cout<<"Off diagonal terms of system+ environment (Part 2a) being evaluated for subspace = "<<subspace<<" which has "<<map_for_states_all[subspace].size()<<" states"<<endl;}
				if (ipr and c1==2) {cout<<"Off diagonal terms of system+ environment (Part 2b) being evaluated for subspace = "<<subspace<<" which has "<<map_for_states_all[subspace].size()<<" states"<<endl;}
				for (int m=0;m<map_for_states_all[subspace].size();m++) // mloop
				{ //m
					ind_1=map_for_states_all[subspace][m];
					vec_1=list_of_converted_vecs[ind_1];
					subspace1_location=inverse_subspace_map[vec_1[c1]];
					state1_location=inverse_sz_map[vec_1[c1]];
					subspace2_location=inverse_subspace_map[vec_1[c2]];
					state2_location=inverse_sz_map[vec_1[c2]];
					vec_2=vec_1;
					ind_subspace1_location=subspace1_location-1;
					ind_subspace2_location=subspace2_location+1;
					if (ind_subspace1_location>=0 and ind_subspace2_location<s_z_block.size())
					{
						for (int j=0;j<s_z_block[ind_subspace1_location].NRows();j++)
						{
							vec_2[c1]=j+shifts[ind_subspace1_location];
							state1_new_location=inverse_sz_map[vec_2[c1]]; // must be j;
							for (int k=0;k<s_z_block[ind_subspace2_location].NRows();k++)
							{
								vec_2[c2]=k+shifts[ind_subspace2_location];
								state2_new_location=inverse_sz_map[vec_2[c2]]; // must be k;
								ind_2=convert_vec_to_ind(vec_2,block_num_states);
								hint=0.0;
								if (s_minus_block[subspace1_location].NRows()>=1 and s_plus_block[subspace2_location].NRows()>=1) 
								{hint+=(j_x*(0.5*s_plus_block[subspace2_location](state2_new_location,state2_location)*s_minus_block[subspace1_location](state1_new_location,state1_location)));}
								if ((ind_1!=ind_2) and (hint!=0.0)) 
								{
									map_for_hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(inverse_sz_map_all[ind_2]);
									hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(hint);
								}
							 }
						 }
					}
					//cout<<"Done S_^+ S-"<<endl;
					vec_2=vec_1;
					ind_subspace1_location=subspace1_location+1;
					ind_subspace2_location=subspace2_location-1;
					if (ind_subspace2_location>=0 and ind_subspace1_location<s_z_block.size())
					{
						for (int j=0;j<s_z_block[ind_subspace1_location].NRows();j++)
						{
							vec_2[c1]=j+shifts[ind_subspace1_location];
							state1_new_location=inverse_sz_map[vec_2[c1]]; // must be j;
							for (int k=0;k<s_z_block[ind_subspace2_location].NRows();k++)
							{
								vec_2[c2]=k+shifts[ind_subspace2_location];
								state2_new_location=inverse_sz_map[vec_2[c2]]; // must be k;
								ind_2=convert_vec_to_ind(vec_2,block_num_states);
								hint=0.0;
								if (s_plus_block[subspace1_location].NRows()>=1 and s_minus_block[subspace2_location].NRows()>=1) 
								{hint+=(j_x*(0.5*s_minus_block[subspace2_location](state2_new_location,state2_location)*s_plus_block[subspace1_location](state1_new_location,state1_location)));}
								if ((ind_1!=ind_2) and (hint!=0.0)) 
								{
									map_for_hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(inverse_sz_map_all[ind_2]);
									hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(hint);
								}
							}
						}
					}
					//cout<<"Done S_- S+"<<endl;
					vec_2=vec_1;
					ind_subspace1_location=subspace1_location;
					ind_subspace2_location=subspace2_location;
					for (int j=0;j<s_z_block[ind_subspace1_location].NRows();j++)
					{
						vec_2[c1]=j+shifts[ind_subspace1_location];
						state1_new_location=inverse_sz_map[vec_2[c1]]; // must be j;
						for (int k=0;k<s_z_block[ind_subspace2_location].NRows();k++)
						{
							vec_2[c2]=k+shifts[ind_subspace2_location];
							state2_new_location=inverse_sz_map[vec_2[c2]]; // must be k;
							ind_2=convert_vec_to_ind(vec_2,block_num_states);
							hint=0.0;
							hint+=(j_z*(s_z_block[subspace2_location](state2_new_location,state2_location)*s_z_block[subspace1_location](state1_new_location,state1_location)));
							if ((ind_1!=ind_2) and (hint!=0.0)) 
							{
								map_for_hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(inverse_sz_map_all[ind_2]);
								hints_all[subspace][inverse_sz_map_all[ind_1]].push_back(hint);
							}
						}
					}
					//cout<<"Done S_z S_z"<<endl;
				}// m loop
			}// if target
		}// c1 loop
		if (map_for_states_all[subspace].size()>0 and target_states[subspace]>0)
		{
			time(&end_0);
			dif=difftime(end_0,start_0);
        		if (ipr) cout<<"Time to compute effective Hamiltonian in relevant spin sector was "<<dif<<" seconds"<<endl;
			if (ipr) cout<<"Diagonalizing effective Hamiltonian of system + environment (Lanczos) in relevant spin sector"<<endl;
			if (map_for_states_all[subspace].size()>1000) 
			{
				/*  PRINT STATEMENT
				cout<<"Calling Lanczos..............."<<endl;
				cout<<endl;
				*/
				lanczos_with_hints_given((150*target_states[subspace])+1,min(target_states[subspace],int(map_for_states_all[subspace].size())),map_for_hints_all[subspace],hints_all[subspace],eigs_all[subspace],opt_evecs_2_all[subspace],false);
			}
			else
			{
			    /*  PRINT STATEMENT
			    cout<<"Calling Exact Matrix Diag.........."<<endl;
				cout<<endl;
				*/
				ed_with_hints_given(map_for_hints_all[subspace],hints_all[subspace],eigs_all[subspace],opt_evecs_2_all[subspace],false);
			}
			if (ipr) cout<<"Energy eigenvalues from Lanczos in subspace with S_z_total = "<<szs_all[subspace]<<endl;
			if (ipr) print_vec_acc(eigs_all[subspace],true,30);
			if (ipr) cout<<endl;
			if (ipr) cout<<endl;
			//if (ipr) cout<<"Inserting Eigs from Lanczos in subsp"<<endl;
			tmp_eigs.insert(tmp_eigs.end(),eigs_all[subspace].begin(),eigs_all[subspace].end());
			
			if (measure)
			{
				cout<<"======================================================================================================================"<<endl;	
				// Now measure in ground state of every spin sector
				cout<<"Now measuring observables in the ground state with S_z "<<szs_all[subspace]<<" of the cactus upto generation "<<gen+1<<endl;
				cout<<"(Note in HJC's notation generation 1 is the bowtie)"<<endl;
				cout<<"(Note in HJC's notation for printing Last(highest) site is boundary, 0 corresponds to center)"<<endl;
				cout<<"(Note in HJC's notation for storage  Last(highest) site is center,   0 corresponds to boundary)"<<endl;
				cout<<"======================================================================================================================"<<endl;	
			        cout<<endl;
				std::vector<double> s_z_avg(gen+2),s_z2_avg(gen+2),siz_scz_avg(gen+2),sip_scm_avg(gen+2),sim_scp_avg(gen+2),si_sc_avg(gen+2);
				double sum=0.0,sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0;
				double szvalue,smvalue,spvalue;
				int j_prime;
                                int block_subspace;
                                int shift_state;
				std:vector<int> tmp_vec;
				target_state.clear();
				for (int j=0;j<opt_evecs_2_all[subspace].NRows();j++) {target_state.push_back(opt_evecs_2_all[subspace](j,0));} // Lowest energy eigenvector
				// Central site measurement
				for (int j=0;j<target_state.size();j++)
				{
					j_prime=j;
					tmp_vec=list_of_converted_vecs[map_for_states_all[subspace][j]];
					szvalue=double(tmp_vec[tmp_vec.size()-1])-spin;	
					sum=sum+(target_state[j])*(target_state[j_prime])*(szvalue);
					sum1=sum1+(target_state[j])*(target_state[j_prime])*(szvalue*szvalue);
					sum2=sum2+(target_state[j])*(target_state[j_prime])*splus_sminus_fn(spin,tmp_vec[tmp_vec.size()-1]);
					sum3=sum3+(target_state[j])*(target_state[j_prime])*sminus_splus_fn(spin,tmp_vec[tmp_vec.size()-1]);
				}
				s_z_avg[gen+1]=sum;
                                s_z2_avg[gen+1]=sum1;
                                siz_scz_avg[gen+1]=s_z2_avg[gen+1];
				sip_scm_avg[gen+1]=sum2;
				sim_scp_avg[gen+1]=sum3;
                                // Measurements of other sites
				for (int st=0;st<gen+1;st++)
				{
					sum=0.0;sum1=0.0;sum2=0.0;sum3=0.0;sum4=0.0;
					for (int j=0;j<target_state.size();j++)
					{
						vec_1=list_of_converted_vecs[map_for_states_all[subspace][j]];
					        block_subspace=inverse_subspace_map[vec_1[0]]; // Arbitrarily choose 1st block since all 4 blocks are equivalent
						shift_state=0;
						for (int k=0;k<block_subspace;k++) {shift_state+=hier.s_z_hier[gen][st][k].NCols();}
						szvalue=double(vec_1[vec_1.size()-1])-spin;	
						vec_2=vec_1;
						for (int k=0;k<hier.s_z_hier[gen][st][block_subspace].NRows();k++)
						{	
							vec_2[0]=k+shift_state;
							j_prime=convert_vec_to_ind(vec_2,block_num_states); //Compact index
                                                        j_prime=inverse_sz_map_all[j_prime]; // Where does the compact index lie in the S_z ground state?
							sum=sum+(target_state[j])*(target_state[j_prime])*(hier.s_z_hier[gen][st][block_subspace](vec_1[0]-shift_state,k));
							sum1=sum1+(target_state[j])*(target_state[j_prime])*(hier.s_z2_hier[gen][st][block_subspace](vec_1[0]-shift_state,k));
							sum2=sum2+(target_state[j])*(target_state[j_prime])*(hier.s_z_hier[gen][st][block_subspace](vec_1[0]-shift_state,k))*szvalue;
						}
					}
					s_z_avg[st]=sum;
                                	s_z2_avg[st]=sum1;
					siz_scz_avg[st]=sum2;
				}
		
				// Si+ Sc-
				for (int st=0;st<gen+1;st++)
				{
					sum=0.0;sum1=0.0;sum2=0.0;
					for (int j=0;j<target_state.size();j++)
					{
						vec_1=list_of_converted_vecs[map_for_states_all[subspace][j]];
					        block_subspace=inverse_subspace_map[vec_1[0]]; // Arbitrarily choose 1st block since all 4 blocks are equivalent
						block_plus_subspace=block_subspace+1;
						if (block_plus_subspace<hier.s_z_hier[gen][st].size())
						{
							shift_state=0;
							for (int k=0;k<block_plus_subspace;k++) {shift_state+=hier.s_z_hier[gen][st][k].NCols();}
							szvalue=double(vec_1[vec_1.size()-1])-spin;	
							vec_2=vec_1;
							smvalue=sminusfn(spin,vec_1[vec_1.size()-1]);
							vec_2[vec_2.size()-1]=vec_1[vec_1.size()-1]-1;
							if (vec_2[vec_2.size()-1]>=0)
							{
								for (int k=0;k<hier.s_plus_hier[gen][st][block_subspace].NRows();k++)
								{	
									vec_2[0]=k+shift_state;
									j_prime=convert_vec_to_ind(vec_2,block_num_states); //Compact index
									j_prime=inverse_sz_map_all[j_prime]; // Where does the compact index lie in the S_z ground state?
									sum=sum+((target_state[j])*(target_state[j_prime])*(hier.s_plus_hier[gen][st][block_subspace](k,inverse_sz_map[vec_1[0]]))*smvalue);
								}
							}
						}
					}
					sip_scm_avg[st]=sum;
				}
				// Si- Sc+
				for (int st=0;st<gen+1;st++)
				{
					sum=0.0;sum1=0.0;sum2=0.0;
					for (int j=0;j<target_state.size();j++)
					{
						vec_1=list_of_converted_vecs[map_for_states_all[subspace][j]];
					        block_subspace=inverse_subspace_map[vec_1[0]]; // Arbitrarily choose 1st block since all 4 blocks are equivalent
						block_minus_subspace=block_subspace-1;
						if (block_minus_subspace>=0)
						{
							shift_state=0;
							for (int k=0;k<block_minus_subspace;k++) {shift_state+=hier.s_z_hier[gen][st][k].NCols();}
							szvalue=double(vec_1[vec_1.size()-1])-spin;	
							vec_2=vec_1;
							spvalue=splusfn(spin,vec_1[vec_1.size()-1]);
							vec_2[vec_2.size()-1]=vec_1[vec_1.size()-1]+1;
							if (vec_2[vec_2.size()-1]<int(2.0*spin+1.0+1.0e-6))
							{
								for (int k=0;k<hier.s_minus_hier[gen][st][block_subspace].NRows();k++)
								{	
									vec_2[0]=k+shift_state;
									j_prime=convert_vec_to_ind(vec_2,block_num_states); //Compact index
									j_prime=inverse_sz_map_all[j_prime]; // Where does the compact index lie in the S_z ground state?
									sum=sum+((target_state[j])*(target_state[j_prime])*(hier.s_minus_hier[gen][st][block_subspace](k,inverse_sz_map[vec_1[0]]))*spvalue);
								}
							}
						}
					}
					sim_scp_avg[st]=sum;
				}

				for (int st=0;st<gen+2;st++) {si_sc_avg[st]=0.5*(sim_scp_avg[st]+sip_scm_avg[st])+siz_scz_avg[st];}

				cout<<"Site               <Sz>             <S_z^2>           <Siz Scz>            <Si+ Sc->            <Si-Sc+>            <Si dot Sj>"<<endl;
				for (int st=gen+1;st>=0;st--)
				{
    					cout<<boost::format("%i         %+.10f       %+.10f       %+.10f       %+.10f       %+.10f      %+.10f") %((gen+1)-st) %s_z_avg[st] % s_z2_avg[st] % siz_scz_avg[st] % sip_scm_avg[st] % sim_scp_avg[st]  % si_sc_avg[st];
					cout<<endl;
				}
				cout<<endl;
				sum=0.0;
				cout<<"Now doing some sanity checks...... i.e. do these measurements make sense?"<<endl;
				cout<<"Checking S_z of all sites combined... This must match up with S_z of the state"<<endl;
				int number_of_sites=pow(2,(gen+2));
				for (int st=0;st<gen+1;st++) {sum=sum+s_z_avg[st]*(number_of_sites); number_of_sites=number_of_sites/2;}
				sum=sum+s_z_avg[gen+1]; // Central site
				cout<<"Sum(S_z) for all sites in Ground state of subspace "<<subspace<<" which has S_z_total "<<szs_all[subspace]<<" is = "<<sum<<endl;
				cout<<"======================================================================================================================"<<endl;	
				cout<<endl;
				cout<<endl;
			} // if measure
		}// if target
		/*  PRINT STATEMENT
		cout<<"Diagonalization of full system and measurements done.... hence deleting maps and hints"<<endl;
		*/
		map_for_hints_all[subspace].clear();
		hints_all[subspace].clear();
	} // subspace

         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Compute RDMS of states you want to average --> organize this by spin
	total_states=0;
	if (ipr) cout<<"Target states"<<endl;
	if (ipr) print_vec(target_states);
	avg_rdms.clear();

	for (int s=0;s<target_states.size();s++)
	{
		for (int ctr=0;ctr<target_states[s];ctr++)
		{	
			if (ipr) cout<<"TRACE: State averaging being performed in subspace with S_z_total = "<<szs_all[s]<<".....considering state "<<ctr<<endl;
			target_state.clear();
			for (int j=0;j<opt_evecs_2_all[s].NRows();j++) {target_state.push_back(opt_evecs_2_all[s](j,ctr));}
			if (ipr) cout<<endl;
		        //cout<<"State whose RDM being calculated"<<endl;
			//print_vec(target_state);
			if (ipr) cout<<"TRACE: Calculating state density matrix for each spin sector"<<endl;
			n_env1=block_num_states[0];
			n_env2=block_num_states[1];	
			if (ipr)
			{
				cout<<"=================================================="<<endl;
				cout<<"             Total N_env     =     "<<n_env1*n_env2<<endl;
				cout<<"=================================================="<<endl;
			}
			rdms.clear();
			rdm_eigs.clear();
			for (int i=0;i<szs.size();i++)       // S_zs_of_system
			{
				rdms.push_back(Matrix(0,0));
				rdms[i].resize(map_for_states_sys[i].size(),map_for_states_sys[i].size());
				rdm_eigs.push_back(std::vector<double>());
				rdm_eigs[i].resize(map_for_states_sys[i].size());
			}

			time(&start);	
			for (int env1=0;env1<n_env1;env1++)
			{
				s_z_env1=szs_on_block[env1];
				for (int env2=0;env2<n_env2;env2++)
				{
					s_z_env2=szs_on_block[env2];
                                        s_z_rdm=szs_all[s]-s_z_env1-s_z_env2;
					
					shift=0;
					flg=false;
					for (int i=0;i<szs.size();i++)
					{
						if (abs(szs[i]-s_z_rdm)<1.0e-10) {selected=i;i=szs.size();flg=true;}
						else {shift=shift+map_for_states_sys[i].size();}
					}

					if (flg)
					{
						for(int state_1_location=0;state_1_location<map_for_states_sys[selected].size();state_1_location++)
						{
							m_1=(n_sys*n_env2*env1)+(n_sys*env2)+map_for_states_sys[selected][state_1_location];
							m_1=inverse_sz_map_all[m_1];	
							
							if (abs(target_state[m_1])>1.0e-12)
							{
							     for(int state_2_location=state_1_location;state_2_location<map_for_states_sys[selected].size();state_2_location++)
							     {
								m_2=(n_sys*env1*n_env2)+(n_sys*env2)+map_for_states_sys[selected][state_2_location];
								m_2=inverse_sz_map_all[m_2];	
								
								rdms[selected](state_1_location,state_2_location)+=(target_state[m_1]*target_state[m_2]);
								rdms[selected](state_2_location,state_1_location)=rdms[selected](state_1_location,state_2_location);
							     }
							}
						}
					}
			        }
			}
			time(&end);		

			dif=difftime(end,start);

			if (ipr) 
			{cout<<"==================================================================="<<endl;
			cout<<"Time to compute Density matrix (S_z arranged) was "<<dif<<" seconds"<<endl;
			cout<<"==================================================================="<<endl;
			}
			
			//cout<<"TRACE: FINISHED calculating reduced density matrices for each spin. Now summing in existing sums"<<endl;
			//cout<<endl;
 			// Lines added below by HJC on March 6 2012... might be removed later
                        //////////////////////////////////////////////////////////////////////////////////////////
			double ee=0.0;
			for (int m=0;m<rdms.size();m++)
		        {
			    if (ipr) cout<<endl;	
				if (ipr) cout<<"(No avg) RDM with S_z_total = "<<szs_all[s]<<".....considering state "<<ctr<<endl;
				/*  PRINT STATEMENT
				cout<< "(No avg) RDM diagonalize for S_z = "<<szs[m]<<endl;
				*/
				symmetric_diagonalize(rdms[m],rdm_eigs[m],rdm_evecs[m]);
				/*  PRINT STATEMENT
				cout<<" (No Avg) Reduced density Matrix eigenvalues in subsp = "<<m<<endl;
				print_vec(rdm_eigs[m]);
				*/
				for (int a=0;a<rdm_eigs[m].size();a++)
				{
					if (abs(rdm_eigs[m][a])>1.0e-7) {ee=ee-(rdm_eigs[m][a]*log(rdm_eigs[m][a]));}
				}
			        if (ipr) cout<<endl;	
			}
		        /*  PRINT STATEMENT
		        cout <<"Entanglement entropy="<<ee<<endl;
		        */
                        //////////////////////////////////////////////////////////////////////////////////////////
			
			if (avg_rdms.size()==0){avg_rdms=rdms;total_states=total_states+1;}
			else
			{
				total_states=total_states+1;
				for (int m=0;m<rdms.size();m++)
				{
					for (int j=0;j<rdms[m].NRows()*rdms[m].NRows();j++)  avg_rdms[m][j]+=rdms[m][j];
				}
			}
		}
	}
        
	// Keep averaging these RDMS as you get them
	if (ipr) cout<<"TRACE: Averaging Density matrices in each spin sector"<<endl;	
	for (int m=0;m<avg_rdms.size();m++)
	{
				for (int j=0;j<rdms[m].NRows()*rdms[m].NRows();j++) avg_rdms[m][j]=avg_rdms[m][j]/double(total_states);
	}

		
	if (ipr)
	{
		cout<<"TRACE: Final State averaged density matrix in each spin sector"<<endl;
		for (int subsp=0;subsp<avg_rdms.size();subsp++) print_real_mat(avg_rdms[subsp]);
	}
        
	// Diagonalize these RDMS store density matrix eigenvalues and make a list of tmp_dmat_eigs which will be sorted later

	time(&start);
	if (ipr) cout<<"TRACE: Diagonalizing state averaged density matrix in each S_z sector         "<<endl;
	for (int subsp=0;subsp<avg_rdms.size();subsp++)
	{
		if(ipr) cout<<"Subspace whose RDM being diagonalized "<<subsp<<"which has S_z "<<szs[subsp]<<endl;
		if(ipr) cout<<"Size of RDM being diagonalized is "<<avg_rdms[subsp].NRows()<<endl;
		/*if (avg_rdms[subsp].NRows()<256) {symmetric_diagonalize(avg_rdms[subsp],rdm_eigs[subsp],rdm_evecs[subsp]);}
		else {matrix_lanczos(max_trunc_states/2,max_trunc_states/2,avg_rdms[subsp],rdm_eigs[subsp],rdm_evecs[subsp],false);}*/
		symmetric_diagonalize(avg_rdms[subsp],rdm_eigs[subsp],rdm_evecs[subsp]);
		//cout<<"Reduced density Matrix eigenvalues in subsp = "<<subsp<<endl;
		//print_vec(rdm_eigs[subsp]);
		//if(ipr) cout<<"Inserting density matrix eigenvalues"<<endl;
		tmp_dmat_eigs.insert(tmp_dmat_eigs.end(),rdm_eigs[subsp].begin(),rdm_eigs[subsp].end());
	}
	if(ipr) cout<<"TRACE: FINISHED diagonalizing state averaged density matrix in all S_z sectors"<<endl;
	time(&end);
	dif=difftime(end,start);	

	if (ipr) cout<<"Total time to diagonalize all the spin arranged RDMs was "<<dif<<" seconds"<<endl;

	//cout<<"Before sort"<<endl;
	//print_vec(tmp_dmat_eigs);
	
	// Sort DM weights in descending order
	sort(tmp_dmat_eigs.begin(),tmp_dmat_eigs.end());
	reverse(tmp_dmat_eigs.begin(),tmp_dmat_eigs.end());	
	//cout<<"After descending sort"<<endl;
	//print_vec(tmp_dmat_eigs);
	
	if(ipr) cout<<"tmp_dmat_eigs.size()="<<tmp_dmat_eigs.size()<<endl;

	//num=min(n_sys,max_trunc_states/2);
	num=min(n_sys,max_trunc_states);
	cutoff=tmp_dmat_eigs[num-1];

	if (ipr) cout<<"Number of states expected to be retained = "<<num<<endl;
	if (ipr) cout<<"Expected density matrix weight cutoff = "<<cutoff<<endl;
	
	final_num=num;	
	for (int i=num;i<tmp_dmat_eigs.size();i++)
	{
		if (final_num<2*max_trunc_states and abs(cutoff)>1.0e-11) 
		{
			cutoff=tmp_dmat_eigs[i];
			final_num=final_num+1;
		}
		else
		{
			i=tmp_dmat_eigs.size();
		}
	}

	if (ipr) cout<<"Corrected number of states to be retained = "<<final_num<<endl;
	if (ipr) cout<<"Corrected Density matrix weight cutoff = "<<cutoff<<endl;

	//num_spaces_sys=map_for_states_sys.size();
	
	for (int subsp=0;subsp<num_spaces_sys;subsp++)
	{
		highest=0;
		if (ipr) {cout<<"rdm_evecs[subsp].NCols() = "<<rdm_evecs[subsp].NCols()<<endl;}
		//print_vec(eigs[subsp]);
		for (int i=rdm_evecs[subsp].NCols()-1;i>=0;i--)
		{
			if (rdm_eigs[subsp][i]>=cutoff) 
			{ 
				highest=highest+1;
				//energies[subsp].push_back(eigs[subsp][i]);
			}
			else
			{
				i=-1;	
			}
		}
		
		//cout<<"highest="<<highest<<endl;
		if (highest>0) {opt_evecs[subsp].resize(rdm_evecs[subsp].NRows(),highest);}
		else {opt_evecs[subsp].resize(0,0);}

		if (ipr) cout<<"TRACE: Number of vectors in subspace"<<subsp<<" = "<<highest<<endl;
		if (highest>0)
		{
			for (int i=0;i<rdm_evecs[subsp].NRows();i++)
			{
				for (int j=0;j<highest;j++)
				{
					opt_evecs[subsp](i,j)=rdm_evecs[subsp](i,rdm_evecs[subsp].NCols()-1-j);
				}
			}
		}
		//cout<<"TRACE: Finished Printing Opt_evecs"<<endl;
	}
	
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (ipr) cout<<"TRACE: Calculating S+,S-,S_z,S_z^2 matrix elements"<<endl; 
	// Calc S+, S- , S_z matrices 
	if (ipr) cout<<"Number of spaces for system = "<<num_spaces_sys<<endl;

	for (int subsp=0;subsp<num_spaces_sys;subsp++)
	{
		plus_one_subsp=subsp+1;
                minus_one_subsp=subsp-1;

		if (plus_one_subsp<opt_evecs.size()) 
		{
		        if (opt_evecs[plus_one_subsp].NCols()>0 and opt_evecs[subsp].NCols()>0)
			{
				s_p_psite[subsp].resize(opt_evecs[plus_one_subsp].NCols(),opt_evecs[subsp].NCols());	
				for (int m=0;m<map_for_states_sys[subsp].size();m++)
				{
					ind_1=map_for_states_sys[subsp][m];
					state_1_location=inverse_sz_map_sys[ind_1];
					// state_1_location must be m
					//cout<<"m = "<<m<<" state_1_location = "<<state_1_location<<endl;
					vec_1=convert_ind_to_vec(ind_1,block_num_states);
                                        sz_ind=vec_1[vec_1.size()-1];
					if (vec_1[vec_1.size()-1]!=int(2*spin+1.0e-6)) 
					{
						vec_2=vec_1;
						vec_2[vec_2.size()-1]=vec_1[vec_1.size()-1]+1;
						ind_2=convert_vec_to_ind(vec_2,block_num_states);
						state_2_location=inverse_sz_map_sys[ind_2];
				
						for(int i=0;i<opt_evecs[plus_one_subsp].NCols();i++)
						{	
							for(int j=0;j<opt_evecs[subsp].NCols();j++) s_p_psite[subsp](i,j)+=opt_evecs[subsp](state_1_location,j)*splusfn(spin,sz_ind)*opt_evecs[plus_one_subsp](state_2_location,i);
						}
					}
				 }
			}
			else
			{
				s_p_psite[subsp].clear();	
			}
		}
		//cout<<"Done with S+"<<endl;
		//{calc_s_plus_matrix_spin(opt_evecs[subsp],opt_evecs[subsp+1],s_p_psite[subsp]);}
		if (minus_one_subsp>=0) 
		{
		        if (opt_evecs[minus_one_subsp].NCols()>0 and opt_evecs[subsp].NCols()>0)
		        {
				s_m_psite[subsp].resize(opt_evecs[minus_one_subsp].NCols(),opt_evecs[subsp].NCols());	
				for (int m=0;m<map_for_states_sys[subsp].size();m++)
				{
					ind_1=map_for_states_sys[subsp][m];
					state_1_location=inverse_sz_map_sys[ind_1];
					// state_1_location must be m
					vec_1=convert_ind_to_vec(ind_1,block_num_states);
					if (vec_1[vec_1.size()-1]!=0) 
					{
						vec_2=vec_1;
                                                sz_ind=vec_1[vec_1.size()-1];
						vec_2[vec_2.size()-1]=vec_1[vec_1.size()-1]-1;
						ind_2=convert_vec_to_ind(vec_2,block_num_states);
						state_2_location=inverse_sz_map_sys[ind_2];
				
						for(int i=0;i<opt_evecs[minus_one_subsp].NCols();i++)
						{	
							for(int j=0;j<opt_evecs[subsp].NCols();j++) s_m_psite[subsp](i,j)+=opt_evecs[subsp](state_1_location,j)*sminusfn(spin,sz_ind)*opt_evecs[minus_one_subsp](state_2_location,i);
						}
					}
						
				}
			}
			else
			{
				s_m_psite[subsp].clear();	
			}
		}
		//cout<<"Done with S-"<<endl;
		//if (subsp-1>=0) {calc_s_minus_matrix_spin(opt_evecs[subsp],opt_evecs[subsp-1],s_m_psite[subsp]);}
		if (opt_evecs[subsp].NCols()>0)
		{	
			s_z_psite[subsp].resize(opt_evecs[subsp].NCols(),opt_evecs[subsp].NCols());
			s_z2_psite[subsp].resize(opt_evecs[subsp].NCols(),opt_evecs[subsp].NCols());
			/*cout<<"Done with resize S_z"<<endl;cout<<"Subspace="<<subsp<<endl;*/
			for (int m=0;m<map_for_states_sys[subsp].size();m++)
			{
					ind_1=map_for_states_sys[subsp][m];
					state_1_location=inverse_sz_map_sys[ind_1];
					vec_1=convert_ind_to_vec(ind_1,block_num_states);
                                        sz_ind=vec_1[vec_1.size()-1];
					for(int i=0;i<opt_evecs[subsp].NCols();i++)
					{	
						for(int j=0;j<opt_evecs[subsp].NCols();j++) s_z_psite[subsp](i,j)+=opt_evecs[subsp](state_1_location,j)*(double(sz_ind)-spin)*opt_evecs[subsp](state_1_location,i);
						for(int j=0;j<opt_evecs[subsp].NCols();j++) s_z2_psite[subsp](i,j)+=opt_evecs[subsp](state_1_location,j)*pow((double(sz_ind)-spin),2.0)*opt_evecs[subsp](state_1_location,i);
					}
			}
		}
		else
		{
			s_z_psite[subsp].resize(0,0);
			s_z2_psite[subsp].resize(0,0);
		}
		
	}
	
	if (ipr) {cout<<"TRACE: FINISHED calculating S+,S-,S_z,S_z^2 matrix elements of central site"<<endl; cout<<endl;}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (ipr) {cout<<"TRACE: Now calculating S+,S-,S_z,S_z^2 matrix elements of all sites at a lower generation (including central site)"<<endl; cout<<endl;}
	 
        if (ipr) {cout<<"TRACE: Presently calculating S+ of all sites at a lower generation (including central site)"<<endl; cout<<endl;}
	compute_splus_higher_gen(gen,szs,map_for_states_sys,states_sys,inverse_sz_map_sys,inverse_subspace_map_sys,inverse_sz_map,inverse_subspace_map,opt_evecs,hier.s_plus_hier[gen],hier.s_z_hier[gen],s_p_psite,hier.s_plus_hier[gen+1]);
	cout<<"Finished computing heirarchical S_+ matrix elements.... Now printing them"<<endl;
	cout<<endl;
	for (int site=0;site<hier.s_plus_hier[gen+1].size();site++)
	{
			for (int subspace=0;subspace<hier.s_plus_hier[gen+1][site].size();subspace++)
			{
				cout<<"Site       ="<<site<<endl;
				cout<<"Subspace   ="<<subspace<<endl;
				cout<<endl;
				print_real_mat(hier.s_plus_hier[gen+1][site][subspace]);
				cout<<endl;
			}
	}
        
        if (ipr) {cout<<"TRACE: Presently calculating S- of all sites at a lower generation (including central site)"<<endl; cout<<endl;}
	compute_sminus_higher_gen(gen,szs,map_for_states_sys,states_sys,inverse_sz_map_sys,inverse_subspace_map_sys,inverse_sz_map,inverse_subspace_map,opt_evecs,hier.s_minus_hier[gen],hier.s_z_hier[gen],s_m_psite,hier.s_minus_hier[gen+1]);
	cout<<"Finished computing heirarchical S_- matrix elements.... Now printing them"<<endl;
	cout<<endl;
	for (int site=0;site<hier.s_minus_hier[gen+1].size();site++)
	{
			for (int subspace=0;subspace<hier.s_minus_hier[gen+1][site].size();subspace++)
			{
				cout<<"Site       ="<<site<<endl;
				cout<<"Subspace   ="<<subspace<<endl;
				cout<<endl;
				print_real_mat(hier.s_minus_hier[gen+1][site][subspace]);
				cout<<endl;
			}
	}
        
	if (ipr) {cout<<"TRACE: Presently calculating Sz of all sites at a lower generation (including central site)"<<endl; cout<<endl;}
	compute_sz_higher_gen(gen,szs,map_for_states_sys,states_sys,inverse_sz_map_sys,inverse_subspace_map_sys,inverse_sz_map,inverse_subspace_map,opt_evecs,hier.s_z_hier[gen],s_z_psite,hier.s_z_hier[gen+1]);
	cout<<"Finished computing heirarchical S_z matrix elements.... Now printing them"<<endl;
	cout<<endl;
	for (int site=0;site<hier.s_z_hier[gen+1].size();site++)
	{
			for (int subspace=0;subspace<hier.s_z_hier[gen+1][site].size();subspace++)
			{
				cout<<"Site       ="<<site<<endl;
				cout<<"Subspace   ="<<subspace<<endl;
				cout<<endl;
				print_real_mat(hier.s_z_hier[gen+1][site][subspace]);
				cout<<endl;
			}
	}
	
	if (ipr) {cout<<"TRACE: Presently calculating Sz^2 of all sites at a lower generation (including central site)"<<endl; cout<<endl;}
	compute_sz2_higher_gen(gen,szs,map_for_states_sys,states_sys,inverse_sz_map_sys,inverse_subspace_map_sys,inverse_sz_map,inverse_subspace_map,opt_evecs,hier.s_z2_hier[gen],s_z2_psite,hier.s_z2_hier[gen+1]);
	cout<<"Finished computing heirarchical S_z^2 matrix elements.... Now printing them"<<endl;
	cout<<endl;
	for (int site=0;site<hier.s_z2_hier[gen+1].size();site++)
	{
			for (int subspace=0;subspace<hier.s_z2_hier[gen+1][site].size();subspace++)
			{
				cout<<"Site       ="<<site<<endl;
				cout<<"Subspace   ="<<subspace<<endl;
				cout<<endl;
				print_real_mat(hier.s_z2_hier[gen+1][site][subspace]);
				cout<<endl;
			}
	}
	
	

	if (ipr) {cout<<"TRACE: Calculating H times transforms"<<endl;}

	htransforms.clear();
	ham_out.clear();

	for(int s=0;s<num_spaces_sys;s++)
	{
		//cout<<"Hamiltonian without transformation in system subspace "<<s<<" which has Sz"<<endl;
                //print_mat_double(hints_sys[s]);
		htransforms.push_back(Matrix(0,0));
		if (opt_evecs[s].NCols()>0) htransforms[s].resize(opt_evecs[s].NRows(),opt_evecs[s].NCols());

		for (int j=0;j<opt_evecs[s].NCols();j++)
		{
			for (int i=0;i<opt_evecs[s].NRows();i++)
			{
				for (int k=0;k<hints_sys[s][i].size();k++)
				{
					htransforms[s](i,j)+=(hints_sys[s][i][k]*opt_evecs[s](map_for_hints_sys[s][i][k],j));
				}
			}
		}
		ham_out.push_back(Matrix(0,0));
		if (opt_evecs[s].NCols()>0) ham_out[s].resize(opt_evecs[s].NCols(),opt_evecs[s].NCols());

		for (int col=0;col<opt_evecs[s].NCols();col++)
		{
			for (int row=0;row<opt_evecs[s].NCols();row++)
			{
				for (int k=0;k<opt_evecs[s].NRows();k++)
				{
					ham_out[s](row,col)+=(opt_evecs[s](k,row)*htransforms[s](k,col));
				}
			}
		}
		//cout<<"Hamiltonian for S_z "<<szs[s]<<endl;
		//print_real_mat(ham_out[s]);
		//cout<<"Eigenvalues of this Hamiltonian for S_z "<<szs[s]<<endl;
		//tmp_en.resize(ham_out[s].NRows());
		//symmetric_diagonalize(ham_out[s],tmp_en,eigenvecs);
		//print_vec(tmp_en);
	}

	inverse_subspace_map_out.clear();
	inverse_sz_map_out.clear();

	for (int subsp=0;subsp<num_spaces_sys;subsp++)
	{
		for (int i=0;i<opt_evecs[subsp].NCols();i++)
		{
			inverse_subspace_map_out.push_back(subsp);
			inverse_sz_map_out.push_back(i);
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_splus_higher_gen(int gen,
			      std::vector<double>                       szs_two_blk,
			      std::vector< std::vector<int> >  		&map_two_blk,
			      std::vector<int>                		states_two_blk,
                              std::vector<int>                		inverse_sz_map_two_blk,
			      std::vector<int>                		inverse_subspace_map_two_blk,
			      std::vector<int>                		inverse_sz_map_one_blk,
			      std::vector<int>                		inverse_subspace_map_one_blk,
                              std::vector<Matrix>             		transforms,
			      std::vector< std::vector<Matrix> >       &s_plus_all_lower,
			      std::vector< std::vector<Matrix> >       &s_z_all_lower,
			      std::vector<Matrix>                      &s_plus_central,
			      std::vector< std::vector<Matrix> >       &s_plus_all_higher)
{
	int 			block_plus_subsp,block_subspace,plus_one_subsp;
        int			ind_1,ind_2,s_1,s_2;
	std::vector<int>        vec_1,vec_2;

	s_plus_all_higher.clear();
	for (int site=0;site<gen+2;site++)
	{
		s_plus_all_higher.push_back(std::vector<Matrix> ());
                for (int i=0;i<s_plus_central.size();i++) s_plus_all_higher[site].push_back(Matrix(s_plus_central[i].NRows(),s_plus_central[i].NCols()));
		if (site==gen+1) {s_plus_all_higher[site]=s_plus_central;}
		else
		{
			for (int subspace=0;subspace<transforms.size();subspace++)
			{
			      plus_one_subsp=subspace+1;
			      if (plus_one_subsp<transforms.size())
			      {
				for (int state_1=0;state_1<transforms[subspace].NCols();state_1++)
				{
				   for (int state_2=0;state_2<transforms[plus_one_subsp].NCols();state_2++)
				   {
					for (int mapind=0;mapind<map_two_blk[subspace].size();mapind++)
					{
						ind_1=map_two_blk[subspace][mapind];
						vec_1=convert_ind_to_vec(ind_1,states_two_blk);	
						vec_2=vec_1;
						block_subspace=inverse_subspace_map_one_blk[vec_1[0]]; // Choose one of the 2 blocks that make up the system. They are equivalent anyway
						block_plus_subsp=block_subspace+1;
						if (block_plus_subsp<s_plus_all_lower[site].size())
						{
							s_1=inverse_sz_map_one_blk[vec_1[0]];
							for (s_2=0;s_2<s_plus_all_lower[site][block_subspace].NRows();s_2++)
							{
								vec_2[0]=s_2;
								for (int i=0;i<block_plus_subsp;i++) {vec_2[0]+=s_z_all_lower[site][i].NCols();}
								ind_2=convert_vec_to_ind(vec_2,states_two_blk);
								s_plus_all_higher[site][subspace](state_2,state_1)+=transforms[subspace](inverse_sz_map_two_blk[ind_1],state_1)*transforms[plus_one_subsp](inverse_sz_map_two_blk[ind_2],state_2)*s_plus_all_lower[site][block_subspace](s_2,s_1);
							}
						}
					 }
				   }
				}
		              }
	         	 }
	        }	
          }
}

void compute_sminus_higher_gen(int gen,
			      std::vector<double>                       szs_two_blk,
			      std::vector< std::vector<int> >  		&map_two_blk,
			      std::vector<int>                		states_two_blk,
                              std::vector<int>                		inverse_sz_map_two_blk,
			      std::vector<int>                		inverse_subspace_map_two_blk,
			      std::vector<int>                		inverse_sz_map_one_blk,
			      std::vector<int>                		inverse_subspace_map_one_blk,
                              std::vector<Matrix>             		transforms,
			      std::vector< std::vector<Matrix> >       &s_minus_all_lower,
			      std::vector< std::vector<Matrix> >       &s_z_all_lower,
			      std::vector<Matrix>                      &s_minus_central,
			      std::vector< std::vector<Matrix> >       &s_minus_all_higher)
{
	int 			block_minus_subsp,block_subspace,minus_one_subsp;
        int			ind_1,ind_2,s_1,s_2;
	std::vector<int>        vec_1,vec_2;

	s_minus_all_higher.clear();
	for (int site=0;site<gen+2;site++)
	{
		s_minus_all_higher.push_back(std::vector<Matrix> ());
                for (int i=0;i<s_minus_central.size();i++) s_minus_all_higher[site].push_back(Matrix(s_minus_central[i].NRows(),s_minus_central[i].NCols()));
		if (site==gen+1) {s_minus_all_higher[site]=s_minus_central;}
		else
		{
			for (int subspace=0;subspace<transforms.size();subspace++)
			{
			      minus_one_subsp=subspace-1;
			      if (minus_one_subsp>=0)
			      {
				for (int state_1=0;state_1<transforms[subspace].NCols();state_1++)
				{
				   for (int state_2=0;state_2<transforms[minus_one_subsp].NCols();state_2++)
				   {
					for (int mapind=0;mapind<map_two_blk[subspace].size();mapind++)
					{
						ind_1=map_two_blk[subspace][mapind];
						vec_1=convert_ind_to_vec(ind_1,states_two_blk);	
						vec_2=vec_1;
						block_subspace=inverse_subspace_map_one_blk[vec_1[0]]; // Choose one of the 2 blocks that make up the system. They are equivalent anyway
						block_minus_subsp=block_subspace-1;
						if (block_minus_subsp>=0)
						{
							s_1=inverse_sz_map_one_blk[vec_1[0]];
							for (s_2=0;s_2<s_minus_all_lower[site][block_subspace].NRows();s_2++)
							{
								vec_2[0]=s_2;
								for (int i=0;i<block_minus_subsp;i++) {vec_2[0]+=s_z_all_lower[site][i].NCols();}
								ind_2=convert_vec_to_ind(vec_2,states_two_blk);
								s_minus_all_higher[site][subspace](state_2,state_1)+=transforms[subspace](inverse_sz_map_two_blk[ind_1],state_1)*transforms[minus_one_subsp](inverse_sz_map_two_blk[ind_2],state_2)*s_minus_all_lower[site][block_subspace](s_2,s_1);
							}
						}
					 }
				   }
				}
		              }
	         	 }
	        }	
          }
}


void compute_sz_higher_gen(int gen,
			      std::vector<double>                       szs_two_blk,
			      std::vector< std::vector<int> >  		&map_two_blk,
			      std::vector<int>                		states_two_blk,
                              std::vector<int>                		inverse_sz_map_two_blk,
			      std::vector<int>                		inverse_subspace_map_two_blk,
			      std::vector<int>                		inverse_sz_map_one_blk,
			      std::vector<int>                		inverse_subspace_map_one_blk,
                              std::vector<Matrix>             		transforms,
			      std::vector< std::vector<Matrix> >       &s_z_all_lower,
			      std::vector<Matrix>                      &s_z_central,
			      std::vector< std::vector<Matrix> >       &s_z_all_higher)
{
	int 			block_subspace;
        int			ind_1,ind_2,s_1,s_2;
	std::vector<int>        vec_1,vec_2;

	s_z_all_higher.clear();
	for (int site=0;site<gen+2;site++)
	{
		s_z_all_higher.push_back(std::vector<Matrix> ());
                for (int i=0;i<s_z_central.size();i++) s_z_all_higher[site].push_back(Matrix(s_z_central[i].NRows(),s_z_central[i].NCols()));
		if (site==gen+1) {s_z_all_higher[site]=s_z_central;}
		else
		{
			for (int subspace=0;subspace<transforms.size();subspace++)
			{
				for (int state_1=0;state_1<transforms[subspace].NCols();state_1++)
				{
				   for (int state_2=0;state_2<transforms[subspace].NCols();state_2++)
				   {
					for (int mapind=0;mapind<map_two_blk[subspace].size();mapind++)
					{
						ind_1=map_two_blk[subspace][mapind];
						vec_1=convert_ind_to_vec(ind_1,states_two_blk);	
						vec_2=vec_1;
						block_subspace=inverse_subspace_map_one_blk[vec_1[0]]; // Choose one of the 2 blocks that make up the system. They are equivalent anyway
						s_1=inverse_sz_map_one_blk[vec_1[0]];
						for (s_2=0;s_2<s_z_all_lower[site][block_subspace].NRows();s_2++)
						{
							vec_2[0]=s_2;
							for (int i=0;i<block_subspace;i++) {vec_2[0]+=s_z_all_lower[site][i].NCols();}
							ind_2=convert_vec_to_ind(vec_2,states_two_blk);
							s_z_all_higher[site][subspace](state_2,state_1)+=transforms[subspace](inverse_sz_map_two_blk[ind_1],state_1)*transforms[subspace](inverse_sz_map_two_blk[ind_2],state_2)*s_z_all_lower[site][block_subspace](s_2,s_1);
						}
					 }
				   }
				}
	         	 }
	        }	
          }
}

void compute_sz2_higher_gen(int gen,
			      std::vector<double>                       szs_two_blk,
			      std::vector< std::vector<int> >  		&map_two_blk,
			      std::vector<int>                		states_two_blk,
                              std::vector<int>                		inverse_sz_map_two_blk,
			      std::vector<int>                		inverse_subspace_map_two_blk,
			      std::vector<int>                		inverse_sz_map_one_blk,
			      std::vector<int>                		inverse_subspace_map_one_blk,
                              std::vector<Matrix>             		transforms,
			      std::vector< std::vector<Matrix> >       &s_z2_all_lower,
			      std::vector<Matrix>                      &s_z2_central,
			      std::vector< std::vector<Matrix> >       &s_z2_all_higher)
{
	int 			block_subspace;
        int			ind_1,ind_2,s_1,s_2;
	std::vector<int>        vec_1,vec_2;

	s_z2_all_higher.clear();
	for (int site=0;site<gen+2;site++)
	{
		s_z2_all_higher.push_back(std::vector<Matrix> ());
                for (int i=0;i<s_z2_central.size();i++) s_z2_all_higher[site].push_back(Matrix(s_z2_central[i].NRows(),s_z2_central[i].NCols()));
		if (site==gen+1) {s_z2_all_higher[site]=s_z2_central;}
		else
		{
			for (int subspace=0;subspace<transforms.size();subspace++)
			{
				for (int state_1=0;state_1<transforms[subspace].NCols();state_1++)
				{
				   for (int state_2=0;state_2<transforms[subspace].NCols();state_2++)
				   {
					for (int mapind=0;mapind<map_two_blk[subspace].size();mapind++)
					{
						ind_1=map_two_blk[subspace][mapind];
						vec_1=convert_ind_to_vec(ind_1,states_two_blk);	
						vec_2=vec_1;
						block_subspace=inverse_subspace_map_one_blk[vec_1[0]]; // Choose one of the 2 blocks that make up the system. They are equivalent anyway
						s_1=inverse_sz_map_one_blk[vec_1[0]];
						for (s_2=0;s_2<s_z2_all_lower[site][block_subspace].NRows();s_2++)
						{
							vec_2[0]=s_2;
							for (int i=0;i<block_subspace;i++) {vec_2[0]+=s_z2_all_lower[site][i].NCols();}
							ind_2=convert_vec_to_ind(vec_2,states_two_blk);
							s_z2_all_higher[site][subspace](state_2,state_1)+=transforms[subspace](inverse_sz_map_two_blk[ind_1],state_1)*transforms[subspace](inverse_sz_map_two_blk[ind_2],state_2)*s_z2_all_lower[site][block_subspace](s_2,s_1);
						}
					 }
				   }
				}
	         	 }
	        }	
          }
}


