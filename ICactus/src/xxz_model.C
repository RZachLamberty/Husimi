#include"global.h"
#include"xxz_model.h"
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

void Spin_Half_XXZ::operator()
                    (std::vector<int> const &config,
                     std::vector< std::vector<int> > &touched_sites_list, 
                     std::vector< std::vector<int> > &vals_on_touched_list,
                     std::vector< complex<double> > &hints_list)
{
    std::vector< std::vector<int> >:: iterator p;
    int first,second;
    
    if (this->j_x !=0)
    {
	    for (p=this->pairs_list.begin();p!=this->pairs_list.end();p++)
	    {
		first=(*p)[0];second=(*p)[1];
		calc_hints_sxsx_sysy(this->j_x,
				     first,
				     second,
				     config,
				     touched_sites_list,
				     vals_on_touched_list,
				     hints_list);
	    }
    }

    if (this->j_z !=0) 
    {
     	calc_hints_szsz_all(this->j_z,
                        this->pairs_list,
                        config,
                        touched_sites_list,
                        vals_on_touched_list,
                        hints_list);
     }

     //if (this->hstag !=0) 
     //{calc_hints_stag_sz(this->hstag,this->eta,config,touched_sites_list,vals_on_touched_list,hints_list);}
     if (this->hstag !=0) 
     {
	calc_hints_stag_sz(this->hstag,this->eta,config,touched_sites_list,vals_on_touched_list,hints_list);
	hints_list[hints_list.size()-1]=0.0;
        for (int i=0;i<this->sp_sites.size();i++)
        {
		hints_list[hints_list.size()-1]+=(this->eta[this->sp_sites[i]]*this->hstag*(double(config[this->sp_sites[i]])-0.5));
	}
     }
 
     
}
////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> get_eta(std::vector< std::vector<int> > list_of_pairs)
{    
   int max=0;
   int num_sites;
   std::vector<int> eta;

   for (int i=0;i<list_of_pairs.size();i++)
   {
            for (int j=0;j<2;j++)
            {
                if (list_of_pairs[i][j]>max) {max=list_of_pairs[i][j];}
            }
    }
    
    num_sites=max+1;

    //cout<<"Num sites="<<num_sites;
    eta.resize(num_sites);         
    eta[0]=1;

    for (int i=0;i<num_sites;i++)
    {
	//cout<<"i="<<i<<endl;
	if (eta[i]==0)
	{
		for (int j=0;j<list_of_pairs.size();j++)
		{
			if (list_of_pairs[j][0]==i and eta[list_of_pairs[j][1]]!=0) 
			{
			      eta[i]=-eta[list_of_pairs[j][1]];
			}
			if (list_of_pairs[j][1]==i and eta[list_of_pairs[j][0]]!=0) 
			{
			      eta[i]=-eta[list_of_pairs[j][0]];
			}
		}	
	}
    }
    return eta;			
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void xxz_setup(string filename, 
               Spin_Half_XXZ &xxz)

{    
    ///////////////////////////////////////////
    // Set the XXZ Hamiltonian
    ///////////////////////////////////////////

    bool found=true;
    string str,str_ret;
    double d,j_x,j_z,hst,spin;
    std::vector< std::vector<int> > pairs;
    int gen;

    // coupling parameters search
    str="j_x";
    search_for(str,filename,str_ret,found);
    if (found){j_x=str_to_d(str_ret);}
    else {j_x=0.0;}

    str="j_z";
    search_for(str,filename,str_ret,found);
    if (found){j_z=str_to_d(str_ret);}
    else{j_z=1.0;}
    
    str="d";  // Anisotropy
    search_for(str,filename,str_ret,found);
    if (found){d=str_to_d(str_ret);}
    else{d=0.0;}
    
    str="spin";
    search_for(str,filename,str_ret,found);
    if (found){spin=str_to_d(str_ret);}
    else{spin=0.5;}
    
    str="hst";
    search_for(str,filename,str_ret,found);
    if (found){hst=str_to_d(str_ret);}
    else{hst=0.0;}
    
    str="gen";
    search_for(str,filename,str_ret,found);
    if (found){gen=str_to_int(str_ret);}
    else {gen=0;}
    
    xxz.init(j_x,j_z,d,spin,gen);
    
    /*str="pairs";
    search_for(str,filename,str_ret,found);
    if (found)
    {
          if (str_ret.substr(0,1)==string("[")) pairs=convert_string_to_vec_of_vec(str_ret);
    }*/    
    
    /*str="nsites";
    search_for(str,filename,str_ret,found);
    if (found){nsites=str_to_int(str_ret);}    
    xxz.init(nsites,j_x,j_z,pairs);*/
    
    /*xxz.init(j_x,j_z,hst,pairs);
    cout<<"H stag="<<hst<<endl;*/

    /*xxz.sp_sites.push_back(0);
    xxz.sp_sites.push_back(2);
    xxz.sp_sites.push_back(5);
    xxz.sp_sites.push_back(6);
    xxz.sp_sites.push_back(9);
    xxz.sp_sites.push_back(10);
    xxz.sp_sites.push_back(12);
    xxz.sp_sites.push_back(13);*/
    
    /*xxz.sp_sites.push_back(1);
    xxz.sp_sites.push_back(3);
    xxz.sp_sites.push_back(4);
    xxz.sp_sites.push_back(7);
    xxz.sp_sites.push_back(8);
    xxz.sp_sites.push_back(11);*/
    
    /*xxz.sp_sites.push_back(0);
    xxz.sp_sites.push_back(2);
    xxz.sp_sites.push_back(9);
    xxz.sp_sites.push_back(10);
    
    
    cout<<"Special sites"<<endl;
    print_vec(xxz.sp_sites);*/
    
    /*cout<<"J_x="<<xxz.j_x<<endl;
    cout<<"J_z="<<xxz.j_z<<endl;
    cout<<"Hst="<<xxz.hstag<<endl;
    cout<<"Eta"<<endl;
    print_vec(xxz.eta);*/
}


