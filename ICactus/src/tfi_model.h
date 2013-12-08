#ifndef TFI_MODEL_HEADER
#define TFI_MODEL_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 XXZ MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////

#include"hamiltonian.h"

class Spin_Half_TFIM: public Ham
{
    private:
    double h_x;
    double j_z;
    //std::vector< std::vector<int> >  pairs_list;
    int sz;

    public:
    void init(int num_sites,
              double j_z,
              double h_x,
              std::vector< std::vector<int> >  list_of_pairs)
    {
        this->num_sites=num_sites;
        this->h_x=h_x;
        this->j_z=j_z;
        this->pairs_list=list_of_pairs;
    }
    
    void init(double j_z,
              double h_x,
              std::vector< std::vector<int> > list_of_pairs)
    {
        int i,j,max;

        max=0;
        this->h_x=h_x;
        this->j_z=j_z;
        this->pairs_list=list_of_pairs;
    
        for (i=0;i<list_of_pairs.size();i++)
        {
            for (j=0;j<2;j++)
            {
                if (list_of_pairs[i][j]>max) {max=list_of_pairs[i][j];}
            }
        }
        this->num_sites=max+1;
    }
    
    void operator()(std::vector<int> const &config,
                    std::vector< std::vector<int> > &touched_sites_list, 
                    std::vector< std::vector<int> > &vals_on_touched_list,
                    std::vector< complex<double> > &hints_list);
    
    Ham* clone() const
    {
        return new Spin_Half_TFIM(*this);
    }    
};

void tfim_setup(std::string filename, 
               Spin_Half_TFIM &xxz);


#endif
