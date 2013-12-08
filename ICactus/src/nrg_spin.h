#ifndef NRG_SPIN_HEADER
#define NRG_SPIN_HEADER

#include"global.h"
#include"matrix_functions.h"
#include"printing_functions.h"
#include"xxz_model.h"

class Hier
{
	public:
	std::vector< std::vector< std::vector<Matrix> > >       s_plus_hier;
	std::vector< std::vector< std::vector<Matrix> > >       s_minus_hier;
	std::vector< std::vector< std::vector<Matrix> > >       s_z_hier;
	std::vector< std::vector< std::vector<Matrix> > >       s_z2_hier;

        void init(int gen)
	{
		for (int i=0;i<gen+2;i++)
		{
			this->s_plus_hier.push_back(std::vector< std::vector<Matrix> > ());
			this->s_minus_hier.push_back(std::vector< std::vector<Matrix> > ());
			this->s_z_hier.push_back(std::vector< std::vector<Matrix> > ());
			this->s_z2_hier.push_back(std::vector< std::vector<Matrix> > ());
		}
	}
};

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
			      std::vector< std::vector<Matrix> >       &s_plus_all_higher);

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
			      std::vector< std::vector<Matrix> >       &s_minus_all_higher);

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
			      std::vector< std::vector<Matrix> >       &s_z_all_higher);

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
			      std::vector< std::vector<Matrix> >       &s_z2_all_higher);

////////////////////////////////////////////////////////////////////////////////////////////////
void nrg_spin(Spin_Half_XXZ                                     xxz,
         std::vector<double>                                    target_szs,
         std::vector<int>                                       target_numbers,
	 std::vector< std::vector<double> > 		   	&eigs,
         std::vector<double>				   	&spins,
    	 std::vector< std::vector<int> >                        &inverse_subspace_map,
    	 std::vector< std::vector<int> >                        &inverse_sz_map,
	 std::vector< std::vector< Matrix > > 	                &transforms,
         std::vector< std::vector< Matrix > >                   &s_plus_block,
         std::vector< std::vector< Matrix > >                   &s_minus_block,
         std::vector< std::vector< Matrix > > 	                &s_z_block,
         std::vector< std::vector< Matrix > > 	                &s_z2_block,
	 std::vector< std::vector< Matrix > > 	                &ham_block,
	 int 							max_states, 
	 bool 							ipr=true);

////////////////////////////////////////////////////////////////////////////////////////////////
void nrg_dmrg_spin_cactus(Spin_Half_XXZ                         xxz,
			  int 					gen,
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
			  bool                                  measure=true,
			  bool					ipr=true);

#endif
