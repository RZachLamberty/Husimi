#ifndef XXZ_MODEL_HEADER
#define XXZ_MODEL_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 XXZ MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////

#include"hamiltonian.h"
std::vector<int> get_eta(std::vector< std::vector<int> > list_of_pairs);

class Spin_Half_XXZ: public Ham
{
    public:
    double           j_x;
    double           j_z;
    double           d;
    double           hstag;
    double           spin;
    int              sz;
    std::vector<int> eta;
    std::vector<int> sp_sites;
    int              gen;
    public:
    void init(int num_sites,
              double j_x,
              double j_z,
              std::vector< std::vector<int> > list_of_pairs)
    {
        this->num_sites=num_sites;
        this->j_x=j_x;
        this->j_z=j_z;
        this->hstag=0.0;
        this->pairs_list=list_of_pairs;
        this->eta=get_eta(this->pairs_list);
        this->sp_sites.clear();
        this->spin=0.5;
        this->gen=0;
    }
    
    void init(double j_x,
              double j_z,
              std::vector< std::vector<int> > list_of_pairs)
    {
        int i,j,max;

        max=0;
        this->j_x=j_x;
        this->j_z=j_z;
        this->hstag=0.0;
        this->pairs_list=list_of_pairs;
    
        for (i=0;i<list_of_pairs.size();i++)
        {
            for (j=0;j<2;j++)
            {
                if (list_of_pairs[i][j]>max) {max=list_of_pairs[i][j];}
            }
        }
        this->num_sites=max+1;
        this->eta=get_eta(this->pairs_list);
        this->sp_sites.clear();
        this->spin=0.5;
        this->gen=0;
    }
    
    void init(double j_x,
              double j_z,
              double hstag,
              std::vector< std::vector<int> > list_of_pairs)
    {
        int i,j,max;

        max=0;
        this->j_x=j_x;
        this->j_z=j_z;
        this->hstag=hstag;
        this->pairs_list=list_of_pairs;
    
        for (i=0;i<list_of_pairs.size();i++)
        {
            for (j=0;j<2;j++)
            {
                if (list_of_pairs[i][j]>max) {max=list_of_pairs[i][j];}
            }
        }
        this->num_sites=max+1;
        this->eta=get_eta(this->pairs_list);
        this->sp_sites.clear();
        this->spin=0.5;
        this->gen=0;
    }
    
    void init(double j_x,double j_z,double spin,int gen)
    {

        this->j_x=j_x;
        this->j_z=j_z;
        this->hstag=0.0;
        this->spin=spin;
        this->gen=gen;
    }
    
    void init(double j_x,double j_z,double d,double spin,int gen)
    {

        this->j_x=j_x;
        this->j_z=j_z;
        this->d=d;
        this->hstag=0.0;
        this->spin=spin;
        this->gen=gen;
    }
    void operator()(std::vector<int> const &config,
                    std::vector< std::vector<int> > &touched_sites_list, 
                    std::vector< std::vector<int> > &vals_on_touched_list,
                    std::vector< complex<double> > &hints_list);
    
    Ham* clone() const
    {
        return new Spin_Half_XXZ(*this);
    }    
};

void xxz_setup(std::string filename, 
               Spin_Half_XXZ &xxz);


#endif
