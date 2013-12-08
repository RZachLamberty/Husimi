#ifndef HENLEY_1D_MODEL_HEADER
#define HENLEY_1D_MODEL_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 XXZ MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////

#include"hamiltonian.h"

class Henley_1d: public Ham
{
    public:
    std::vector<double> j;
    //std::vector< std::vector<int> > pairs_list;

    public:
    void init(std::vector<double> j,
              std::vector< std::vector<int> > list_of_pairs)
    {
        int i,k,max;

        max=0;
        this->j=j;
        this->pairs_list=list_of_pairs;
    
        for (i=0;i<list_of_pairs.size();i++)
        {
            for (k=0;k<2;k++)
            {
                if (list_of_pairs[i][k]>max) {max=list_of_pairs[i][k];}
            }
        }
        this->num_sites=max+1;
    }
    
    
    void operator()(std::vector<int> const &config,
                    std::vector< std::vector<int> > &touched_sites_list, 
                    std::vector< std::vector<int> > &vals_on_touched_list,
                    std::vector< complex<double> > &hints_list);
    
    Ham* clone() const
    {return new Henley_1d(*this);}    
};

void henley_1d_setup(std::string filename, 
               Henley_1d &henley_1d);


#endif
