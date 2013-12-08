#ifndef ED_HEADER
#define ED_HEADER

#include"hamiltonian.h"
#include"global.h"
#include"matrix_functions.h"
#include"number_functions.h"
#include"printing_functions.h"
#include"bethe_lapack_interface.h"
#include"math_utilities.h"
#include"hamiltonian_spin_functions.h"

double measure_sz(std::vector<int> const &config);

double compute_spin(double s_s_plus_one);

void ed_get_eigs(Ham &h, std::vector<double> &eigs);

void lanczos_no_sym_get_eigs(Ham &h, int it, std::vector<double> &eigs);

void lanczos_spin_sym(Ham &h,int iterations,
		      std::vector<double> &eigs,
		      std::vector<double> &szs,
		      std::vector<double> &spins,
                      bool measure_s=false, bool ipr=true, bool only_sz_non_negative=false);

void lanczos_given_map(Ham &h,
                      int iterations,
                      std::vector<int> const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &spins,
                      bool measure_s=false, bool ipr=true);

void lanczos_with_hints_given(int iterations,int how_may_vecs, 
		      	std::vector< std::vector<int> > const &map,
		      	std::vector< std::vector<double> > const &hints,
                      	std::vector<double> &eigs,
			Matrix &eigenvecs, bool ipr=true);

void ed_with_hints_given(std::vector< std::vector<int> > 		const &map,
		      	 std::vector< std::vector<double> > 	        const &hints,
                      	 std::vector<double> 			        &eigs,
			 Matrix 					&eigenvecs,
			 bool 					        ipri=true);

void lanczos_requested_sz(Ham &h,
                      int iterations, 
                      std::vector<double> &eigs,
		      double s_z,
		      std::vector<double> &spins,
                      bool measure_s);

#endif
