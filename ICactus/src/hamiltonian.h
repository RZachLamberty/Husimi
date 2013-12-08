#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER

#include"global.h"
using namespace std;

class Ham{

    public:
    int num_sites;
    std::vector< std::vector<int> > pairs_list;
    virtual void operator()
                       (std::vector<int> const &config, 
                        std::vector< std::vector<int> > &touched_sites_list, 
                        std::vector< std::vector<int> > &vals_on_touched_list, 
                        std::vector< complex<double> > &hints_list) {};
    
    virtual void init(){};
    
    virtual Ham* clone() const=0;     
            
};

#endif
