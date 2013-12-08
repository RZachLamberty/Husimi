#ifndef LANCZOS_MORE_HEADER
#define LANCZOS_MORE_HEADER

#include"hamiltonian.h"
#include"global.h"
#include"matrix_functions.h"
#include"number_functions.h"
#include"printing_functions.h"
#include"bethe_lapack_interface.h"
#include"math_utilities.h"
#include"hamiltonian_spin_functions.h"


void lanczos_given_map_evecs(Ham &h,
                      int iterations,
                      std::vector<int>  const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &spins,
                      std::vector< std::vector<double> > &lowest_evecs,
		      int states_requested=2,
		      bool ipr=false);

#endif
