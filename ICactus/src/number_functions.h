#ifndef NUMBER_FUNCTIONS_HEADER
#define NUMBER_FUNCTIONS_HEADER

#include"global.h"

//////////////////////////////////////////////////////////////////////////////
// Useful number functions  
//////////////////////////////////////////////////////////////////////////////

int n_choose_k(int n,int k);

void constrained_dets(int num_sites,int num_ones,std::vector<int> &dets);

std::vector<int> convert_ind_to_vec(int index, 
				    std::vector<int> nstates);
int convert_vec_to_ind(std::vector<int> const &vec, 
		       std::vector<int> nstates);
void get_adj_list_from_pairs(std::vector< std::vector<int> > const &pairs,
			     std::vector< std::vector<int> > &adj_list);
void convert_num_to_vec(int num, int base, int num_bits, 
			std::vector<int> &config);
int convert_vec_to_num(std::vector<int> &config, int base);
int str_to_int(std::string str);
double str_to_d(std::string str);
bool str_to_bool(std::string str);
template <class T>
std::string to_string (T const &t);
std::vector<int> convert_string_to_vec(std::string const &s);
std::vector<double> convert_string_to_vec_double(std::string const &s);
std::vector< std::vector<int> > convert_string_to_vec_of_vec(std::string const &s);
std::string dtos(double dbl);

#endif

