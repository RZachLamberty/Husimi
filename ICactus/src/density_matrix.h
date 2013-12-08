#ifndef DENSITY_MATRIX_HEADER
#define DENSITY_MATRIX_HEADER

#include"global.h"
#include"matrix.h"
#include"matrix_functions.h"
#include"printing_functions.h"

void density_matrix_calc_with_map(std::vector<double> const &evec, 
			 std::vector<int> const &map,
			 std::vector<int> const &inverse_map,
			 std::vector<int> const &all_states,
			 std::vector<int> const &which_indices_to_keep,
			 Matrix &den_mat);

void off_diag_density_matrix_calc_with_map(std::vector<double> const &evec_1, 
			std::vector<int> const &map_1,
			std::vector<double> const &evec_2, 
			std::vector<int> const &map_2,
			std::vector<int> const &inverse_map,
			std::vector<int> const &all_states,
			std::vector<int> const &which_indices_to_keep,
			Matrix &den_mat);

void density_matrix_calc(std::vector<double> const &evec, 
			 std::vector<int> const &all_states,
			 std::vector<int> const &which_indices_to_keep,
			 Matrix &den_mat);


#endif
