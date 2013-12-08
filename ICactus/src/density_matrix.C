#include"density_matrix.h"
#include"number_functions.h"
#include"printing_functions.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////
void density_matrix_calc_with_map(std::vector<double> const &evec, 
			 std::vector<int> const &map,
			 std::vector<int> const &inverse_map,
			 std::vector<int> const &all_states,
			 std::vector<int> const &which_indices_to_keep,
			 Matrix &den_mat)
{

	time_t 			start,end;
	double 			dif;
	int 			n=1,n_dm=1,n_env=1;
	int 			c_ind_1,c_ind_2,m_1,m_2;
	int 			size;
	std::vector<int> 	states,env_states,indices_copy,site_indices,env_indices;
	std::vector<int> 	vec_1,vec_2,vec_env,vec_1_t(all_states.size()),vec_2_t(all_states.size());

	size=pow(2,all_states.size());
	indices_copy=which_indices_to_keep;
	sort(indices_copy.begin(),indices_copy.end());

	site_indices=indices_copy;
	
	for (int i=0;i<all_states.size();i++) n=n*all_states[i];
	
	for (int i=0;i<all_states.size();i++) 
	{
		if (indices_copy[0]==i)
		{
			n_dm=n_dm*all_states[i];
		 	states.push_back(all_states[i]);
		 	indices_copy.erase(indices_copy.begin());
		}
		else
		{
			n_env=n_env*all_states[i];
		 	env_states.push_back(all_states[i]);
			env_indices.push_back(i);
		}
	}
	
	cout<<"========================================="<<endl;
	cout<<"             Total N_dm      =     "<<n_dm<<endl;
	cout<<"             Total N_env     =     "<<n_env<<endl;
	cout<<"             Total N_states  =     "<<n<<endl;
	cout<<"========================================="<<endl;

	den_mat.clear();
	den_mat.resize(n_dm,n_dm);
	
	cout<<"==================================="<<endl;
	cout<<"	          Site indices	          "<<endl;
	print_vec(site_indices);
	cout<<"		  Env indices		  "<<endl;
	print_vec(env_indices);
	cout<<"==================================="<<endl;
				
	
        time(&start);		
	for (int env=0;env<n_env;env++)
	{
		vec_env=convert_ind_to_vec(env,env_states);
		for (int i=0;i<env_states.size();i++) vec_1_t[env_indices[i]]=vec_env[i];

		vec_2_t=vec_1_t;

		for(int ind1=0;ind1<n_dm;ind1++)
		{
			vec_1=convert_ind_to_vec(ind1,states);
			for (int i=0;i<states.size();i++) vec_1_t[site_indices[i]]=vec_1[i];
				
                        c_ind_1=convert_vec_to_ind(vec_1_t,all_states);
		        m_1=inverse_map[c_ind_1];

                        if (map[m_1]==c_ind_1)
			{
				if (abs(evec[m_1])>1.0e-10)
				{
					for(int ind2=ind1;ind2<n_dm;ind2++)
					{
						vec_2=convert_ind_to_vec(ind2,states);
						for (int i=0;i<states.size();i++) vec_2_t[site_indices[i]]=vec_2[i];
					
						c_ind_2=convert_vec_to_ind(vec_2_t,all_states);
						m_2=inverse_map[c_ind_2];
							
						if (map[m_2]==c_ind_2) den_mat(ind1,ind2)+=(evec[m_1]*evec[m_2]);

						den_mat(ind2,ind1)=den_mat(ind1,ind2);
					}
				}
			}
		}
	}
        time(&end);		

	dif=difftime(end,start);

        cout<<"==================================================================="<<endl;
        cout<<"Total to compute DM (diagonal) was "<<dif<<" seconds"<<endl;
        cout<<"==================================================================="<<endl;

}

//////////////////////////////////////////////////////////////////////////////
void off_diag_density_matrix_calc_with_map(std::vector<double> const &evec_1, 
			std::vector<int> const &map_1,
			std::vector<double> const &evec_2, 
			std::vector<int> const &map_2,
			std::vector<int> const &inverse_map,
			std::vector<int> const &all_states,
			std::vector<int> const &which_indices_to_keep,
			Matrix &den_mat)
{

	time_t 			start,end;
	double 			dif;
	int 			i,ind1,ind2;
	int 			n=1,n_dm=1,n_env=1;
	int 			c_ind_1,c_ind_2,m_1,m_2;
	int 			size;
	std::vector<int> 	states,env_states,indices_copy,site_indices,env_indices;
	std::vector<int> 	vec_1,vec_2,vec_env,vec_1_t(all_states.size()),vec_2_t(all_states.size());

	size=pow(2,all_states.size());

	indices_copy=which_indices_to_keep;
	sort(indices_copy.begin(),indices_copy.end());

	site_indices=indices_copy;
	
	for (i=0;i<all_states.size();i++) {n=n*all_states[i];}
	
	for (i=0;i<all_states.size();i++) 
	{
		if (indices_copy[0]==i)
		{
			n_dm=n_dm*all_states[i];
		 	states.push_back(all_states[i]);
		 	indices_copy.erase(indices_copy.begin());
		}
		else
		{
			n_env=n_env*all_states[i];
		 	env_states.push_back(all_states[i]);
			env_indices.push_back(i);
		}
	}
	
	cout<<"========================================="<<endl;
	cout<<"             Total N_dm      =     "<<n_dm<<endl;
	cout<<"             Total N_env     =     "<<n_env<<endl;
	cout<<"             Total N_states  =     "<<n<<endl;
	cout<<"========================================="<<endl;

	den_mat.clear();
	den_mat.resize(n_dm,n_dm);
	
	cout<<"==================================="<<endl;
	cout<<"	          Site indices	          "<<endl;
	print_vec(site_indices);
	cout<<"		  Env indices		  "<<endl;
	print_vec(env_indices);
	cout<<"==================================="<<endl;
				

	time(&start);
	for (int env=0;env<n_env;env++)
	{
		vec_env=convert_ind_to_vec(env,env_states);
		for (i=0;i<env_states.size();i++) vec_1_t[env_indices[i]]=vec_env[i];

		vec_2_t=vec_1_t;

		for(ind1=0;ind1<n_dm;ind1++)
		{
			vec_1=convert_ind_to_vec(ind1,states);
			for (i=0;i<states.size();i++) vec_1_t[site_indices[i]]=vec_1[i];
			
			c_ind_1=convert_vec_to_ind(vec_1_t,all_states);
			m_1=inverse_map[c_ind_1];

			if (map_1[m_1]==c_ind_1)
			{
				if(abs(evec_1[m_1])>1.0e-10)
				{
					for(ind2=0;ind2<n_dm;ind2++)
					{
						vec_2=convert_ind_to_vec(ind2,states);
						for (i=0;i<states.size();i++) vec_2_t[site_indices[i]]=vec_2[i];
					
						c_ind_2=convert_vec_to_ind(vec_2_t,all_states);

						m_2=inverse_map[c_ind_2];
							
						if (map_2[m_2]==c_ind_2) den_mat(ind1,ind2)+=(evec_1[m_1]*evec_2[m_2]);
					}
				}
			}
		}
	}
   	time (&end);
   	
	dif=difftime(end,start);
        
	cout<<"==================================================================="<<endl;
        cout<<"Total to compute DM (off diagonal) was "<<dif<<" seconds"<<endl;
        cout<<"==================================================================="<<endl;

}

/////////////////////////////////////////////////////////////////////////////////////////
void density_matrix_calc(std::vector<double> const &evec, 
			 std::vector<int> const &all_states,
			 std::vector<int> const &which_indices_to_keep,
			 Matrix &den_mat)
{

        time_t			start,end;
	int 			i,ind1,ind2;
	int 			n=1,n_dm=1,n_env=1;
	int 			c_ind_1,c_ind_2,m_1,m_2;
	int 			size;
        double 			dif;
	std::vector<int> 	states,env_states,indices_copy,site_indices,env_indices;
	std::vector<int> 	vec_1,vec_2,vec_env,vec_1_t(all_states.size()),vec_2_t(all_states.size());

	size=pow(2,all_states.size());

	indices_copy=which_indices_to_keep;
	sort(indices_copy.begin(),indices_copy.end());
	site_indices=indices_copy;
	
	for (i=0;i<all_states.size();i++) {n=n*all_states[i];}
	
	for (i=0;i<all_states.size();i++) 
	{
		if (indices_copy[0]==i)
		{
			n_dm=n_dm*all_states[i];
		 	states.push_back(all_states[i]);
		 	indices_copy.erase(indices_copy.begin());
		}
		else
		{
			n_env=n_env*all_states[i];
		 	env_states.push_back(all_states[i]);
			env_indices.push_back(i);
		}
	}
	
	cout<<"========================================="<<endl;
	cout<<"             Total N_dm      =     "<<n_dm<<endl;
	cout<<"             Total N_env     =     "<<n_env<<endl;
	cout<<"             Total N_states  =     "<<n<<endl;
	cout<<"========================================="<<endl;

	den_mat.clear();
	den_mat.resize(n_dm,n_dm);
	
	cout<<"==================================="<<endl;
	cout<<"	          Site indices	          "<<endl;
	print_vec(site_indices);
	cout<<"		  Env indices		  "<<endl;
	print_vec(env_indices);
	cout<<"==================================="<<endl;

        time(&start);

        std::vector< std::vector<int> > inds;
	for(ind1=0;ind1<n_dm;ind1++)
	{
	 	vec_1=convert_ind_to_vec(ind1,states);
		inds.push_back(vec_1);
        }
		
	for (int env=0;env<n_env;env++)
	{
		vec_env=convert_ind_to_vec(env,env_states);
		for (i=0;i<env_states.size();i++) vec_1_t[env_indices[i]]=vec_env[i];

		vec_2_t=vec_1_t;

		for(ind1=0;ind1<n_dm;ind1++)
		{
			//vec_1=convert_ind_to_vec(ind1,states);
			vec_1=inds[ind1];
			
			for (i=0;i<states.size();i++) vec_1_t[site_indices[i]]=vec_1[i];
				
                        c_ind_1=convert_vec_to_ind(vec_1_t,all_states);
                    
			if (abs(evec[c_ind_1])>1.0e-10)
			{
				for(ind2=ind1;ind2<n_dm;ind2++)
				{
					//vec_2=convert_ind_to_vec(ind2,states);
					vec_2=inds[ind2];
					
					for (i=0;i<states.size();i++) vec_2_t[site_indices[i]]=vec_2[i];
				
					c_ind_2=convert_vec_to_ind(vec_2_t,all_states);

					den_mat(ind1,ind2)+=(evec[c_ind_1]*evec[c_ind_2]);

					den_mat(ind2,ind1)=den_mat(ind1,ind2);
				}
			}
		}
	}
	time(&end);
	dif=difftime(end,start);
        
	cout<<"==================================================================="<<endl;
        cout<<"Total time to compute DM was "<<dif<<" seconds"<<endl;
        cout<<"==================================================================="<<endl;

}


