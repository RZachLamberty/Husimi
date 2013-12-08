#include"global.h"
#include"tfi_model.h"
#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
#include"search_for.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//                                 TRANSVERSE FIELD ISING MODEL
//
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Spin_Half_TFIM::operator()
                    (std::vector<int> const &config,
                     std::vector< std::vector<int> > &touched_sites_list, 
                     std::vector< std::vector<int> > &vals_on_touched_list,
                     std::vector< complex<double> > &hints_list)
{
    std::vector< std::vector<int> >:: iterator p;
    int site;

    calc_hints_szsz_all(this->j_z,
                        this->pairs_list,
                        config,
                        touched_sites_list,
                        vals_on_touched_list,
                        hints_list);

    for (site=0;site<config.size();site++)
    {
        calc_hints_sx(this->h_x,
                      site, 
                      config,
                      touched_sites_list,
                      vals_on_touched_list,
                      hints_list);
    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////

void tfim_setup(string filename, 
               Spin_Half_TFIM &tfim)

{    
    ///////////////////////////////////////////
    // Set the TFIM Hamiltonian
    ///////////////////////////////////////////

    bool found=true;
    string str,str_ret;
    double h_x,j_z;
    std::vector< std::vector<int> > pairs;
    int nsites;

    // coupling parameters search
    str="h_x";
    search_for(str,filename,str_ret,found);
    
    if (found)
    {
        h_x=str_to_d(str_ret);
    }
    else
    {
        h_x=0.0;
    }

    str="j_z";
    search_for(str,filename,str_ret,found);
    
    if (found)
    {
        j_z=str_to_d(str_ret);
    }
    else
    {
        j_z=1.0;
    }
    
    str="pairs";
    search_for(str,filename,str_ret,found);
    if (found)
    {
          if (str_ret.substr(0,1)==string("["))
          {
            pairs=convert_string_to_vec_of_vec(str_ret);
          }
    }    
    
    /*str="nsites";
    search_for(str,filename,str_ret,found);
    if (found)
    {
        nsites=str_to_int(str_ret);
    }    
    
    tfim.init(nsites,j_z,h_x,pairs);*/
    
    tfim.init(j_z,h_x,pairs);

    
}

