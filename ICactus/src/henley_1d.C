#include"global.h"
#include"henley_1d.h"
#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
#include"search_for.h"
#include"printing_functions.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//                                 XXZ MODEL
//
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Henley_1d::operator()
                    (std::vector<int> const &config,
                     std::vector< std::vector<int> > &touched_sites_list, 
                     std::vector< std::vector<int> > &vals_on_touched_list,
                     std::vector< complex<double> > &hints_list)
{
    std::vector< std::vector<int> >:: iterator p;
    int first,second;
    int ind;

    ind=0;

    for (p=this->pairs_list.begin();p!=this->pairs_list.end();p++)
    {
        //cout<<"j="<<(this->j)[ind];
	first=(*p)[0];second=(*p)[1];
        calc_hints_sxsx_sysy((this->j)[ind],
                             first,second,config,touched_sites_list,
                             vals_on_touched_list,hints_list);
    	calc_hints_szsz((this->j)[ind],first,second,config,touched_sites_list,
			vals_on_touched_list,hints_list);
	ind=ind+1;
    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////

void henley_1d_setup(string filename, Henley_1d &henley_1d)

{    
    ///////////////////////////////////////////
    // Set the Henley_1d Hamiltonian
    ///////////////////////////////////////////

    bool found=true;
    string str,str_ret;
    std::vector<double> j;
    std::vector< std::vector<int> > pairs;

    // coupling parameters search
    str="j";
    search_for(str,filename,str_ret,found);
    
    if (found){j=convert_string_to_vec_double(str_ret);}
    else {return;}

    str="pairs";
    search_for(str,filename,str_ret,found);
    if (found)
    {
          if (str_ret.substr(0,1)==string("["))
          {pairs=convert_string_to_vec_of_vec(str_ret);}
    }
    else {return;}
    
    if (j.size()==pairs.size()) {henley_1d.init(j,pairs);}
    else {cout<<"Expect an error message because j size doesnt match with pair size"<<endl;}
    
    //print_vec(j);
    //print_mat_int(pairs);
    //print_mathematica_pairs(pairs);
    //print_mathematica_pairs(henley_1d.pairs_list);

}


