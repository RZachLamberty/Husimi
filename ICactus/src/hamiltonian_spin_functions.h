#ifndef HAMILTONIAN_SPIN_FUNCTIONS_HEADER
#define HAMILTONIAN_SPIN_FUNCTIONS_HEADER

#include"global.h"
using namespace std;

void calc_hints_sxsx_sysy(double coupling, 
                        int first, 
                        int second, 
                        std::vector<int> const &config,
                        std::vector< std::vector<int> > &touched_sites_list,
                        std::vector< vector<int> > &vals_on_touched_list,
                        std::vector< complex<double> > &hints_list);

void calc_hints_szsz(double coupling, 
                    int first, 
                    int second, 
                    std::vector<int> const &config,
                    std::vector< std::vector<int> > &touched_sites_list,
                    std::vector< std::vector<int> > &vals_on_touched_list,
                    std::vector< complex<double> > &hints_list);

void calc_hints_szsz_all(double coupling, 
                        std::vector< std::vector<int> > const &pairs, 
                        std::vector<int> const &config,
                        std::vector< std::vector<int> > &touched_sites_list,
                        std::vector< std::vector<int> > &vals_on_touched_list,
                        std::vector< complex<double> > &hints_list);

void calc_hints_stag_sz(double hstag, 
                         std::vector<int> const &eta, 
                         std::vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list);

void calc_hints_sx(double coupling, 
                   int site, 
                   std::vector<int> const &config,
                   std::vector< std::vector<int> > &touched_sites_list,
                   std::vector< vector<int> > &vals_on_touched_list,
                   std::vector< complex<double> > &hints_list);

void compute_c_plus_minus_i(int plus_or_minus,
			    std::vector<double> const &vec_0,
			    std::vector<double> const &vec_1,
   			    std::vector<int> const &maps_0,
   			    std::vector<int> const &inverse_map,
			    std::vector<double> &c_i);

void compute_si_sj(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj);

void compute_si_plus_sj_minus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_p_sj_m);

void compute_si_minus_sj_plus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_m_sj_p);

void compute_si_z_sj_z(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
		   Matrix &si_z_sj_z);

void compute_si_minus_sj_plus_sk_z(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   std::vector<double> &three_pt);

void compute_si_minus_sj_minus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj);

#endif
