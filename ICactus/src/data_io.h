#ifndef DATA_IO_HEADER
#define DATA_IO_HEADER

#include"global.h"

template <class T>
void dump_function_auto(std::vector<T> const &y, 
                        std::string filename);

template <class T>
void dump_function(std::vector<T> const &x,
                   std::vector<T> const &y, 
                   std::string filename);

void dump_multiple_functions(std::vector<double> const &x,
                             std::vector< std::vector<double> > const &funcs, 
                             std::string filename);
void dump_function_double(std::vector<double> const &x,
                   std::vector<double> const &y, 
                   std::string filename);

#endif

