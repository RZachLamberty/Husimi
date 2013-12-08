#include"global.h"

// Search
#include"search_for.h"

// Hamiltonian
#include"hamiltonian.h"

// Heisenberg
#include"tfi_model.h"
#include"xxz_model.h"
#include"henley_1d.h"

// ED
#include"ed.h"
#include"lanczos_more.h"

// NRG and DMRG
#include"nrg_spin.h"
#include"density_matrix.h"

//GET AI FROM CI
#include"get_ai.h"

using namespace std;

int main()
{
    time_t start,end;
    double dif;

    bool found;
    string filename="./input_files/input.txt";
    string str_ret;

 /////////////////////////////////////////////////////////////////////////////
 //                    HAMILTONIAN READ AND SETUP calls
 /////////////////////////////////////////////////////////////////////////////
    string hamiltonian;
    Ham *ham=NULL;
    Spin_Half_TFIM tfim;
    Spin_Half_XXZ xxz;
    Henley_1d henley_1d;

    bool ham_found;

    cout<<endl;
    cout<<"========================================="<<endl;
    cout<<"I am reading the Hamiltonian information "<<endl;
    cout<<"========================================="<<endl;

    search_for(string("hamiltonian"),filename,hamiltonian,ham_found);
    if (ham_found) 
    {
        if (hamiltonian.compare("tfim")==0)
        {    
            tfim_setup(filename,tfim);cout<<endl;
            cout<<"Done with Spin Half Transverse Field Ising model (TFIM) setup"<<endl;
            ham=tfim.clone();
        }
        else if (hamiltonian.compare("xxz")==0)
        {
            xxz_setup(filename,xxz);cout<<endl;
            cout<<"Done with XXZ model setup"<<endl;
            ham=xxz.clone();
            cout<<"Done with XXZ cloning"<<endl;
        }
        else if (hamiltonian.compare("henley_1d")==0)
        {
            henley_1d_setup(filename,henley_1d);cout<<endl;
            cout<<"Done with Henley's 1D model setup"<<endl;
            ham=henley_1d.clone();
            //cout<<"Cloned Henley's 1D model setup"<<endl;
        }
	else
        {
            cout<<endl;cout<<"I could not find the requested hamiltonian"<<endl;
            //return 0;
        }
    }
 
/////////////////////////////////////////////////////////////////////////////
//                              EXACT DIAGONALIZATION
/////////////////////////////////////////////////////////////////////////////
    std::vector<double> eigs,szs,spins;

    if (ham_found)
    {
	    //print_mathematica_pairs((*ham).pairs_list);
    	    //cout<<"FINISHED Printing pairs list"<<endl;
	    
            search_for(string("diagonalize_all"),filename,str_ret,found);
	    if (str_to_bool(str_ret))
	    {		
        	if (hamiltonian.compare("tfim")==0)
		{
			lanczos_no_sym_get_eigs(*ham,512,eigs);
			cout<<"ED energies (Lanczos)"<<endl;
			print_vec(eigs,true);
		}
        	
		if (hamiltonian.compare("xxz")==0 or hamiltonian.compare("henley_1d")==0)
		{
			bool iprNow = false;
			lanczos_spin_sym(*ham,200,eigs,szs,spins,false,iprNow,true);
			cout<<"ED energies and Szs (Lanczos)"<<endl;
			//print_vec(eigs,true);
			for (int n=0;n<eigs.size();n++)
			{cout<<boost::format("%+.10f") %eigs[n]<<"  "<<boost::format("%+.3f") %szs[n]<<endl;}
			cout<<"ED energies -Mathematica List (Lanczos)"<<endl;
			for (int n=0;n<min(64,int(eigs.size()));n++)
			{cout<<"{"<<xxz.j_z<<","<<eigs[n]-eigs[0]<<"},";}
		        cout<<endl;	
			cout<<"ED Szs -Mathematica List (Lanczos)"<<endl;
			for (int n=0;n<min(64,int(eigs.size()));n++)
			{cout<<"{"<<xxz.j_z<<","<<szs[n]<<"},";}
		        cout<<endl;	
		}

		cout<<"Ham sites"<<(*ham).num_sites<<endl;

	    	if ((*ham).num_sites<=13)
		{ ed_get_eigs(*ham,eigs);
		  cout<<"ED energies (full)"<<endl;
		  print_vec(eigs,true);
		}
	    }   
	    search_for(string("diagonalize"),filename,str_ret,found);
	    if (str_to_bool(str_ret))
	    {		
	    	//if ((*ham).num_sites<=11){ ed_get_eigs(*ham,eigs);}
	    	//else{lanczos_no_sym_get_eigs(*ham,512,eigs);}
        	if (hamiltonian.compare("xxz")==0)
		{
			lanczos_no_sym_get_eigs(*ham,512,eigs);
			cout<<"ED energies (Lanczos)"<<endl;
			print_vec(eigs,true);
		}
        	
		if (hamiltonian.compare("xxz")==0 or hamiltonian.compare("henley_1d")==0)
		{
			//lanczos_spin_sym(*ham,64,eigs,szs,spins,false);
			//cout<<"ED energies (Lanczos)"<<endl;
			//print_vec(eigs,true);
			
			if ((*ham).num_sites<24) 
			{
				lanczos_requested_sz(*ham,128,eigs,0.0,spins,true);
				cout<<"ED energies (Lanczos requested S_z=0)"<<endl;
				for (int n=0;n<eigs.size();n++)
				{cout<<eigs[n]<<"  "<<spins[n]<<endl;}
				print_vec(eigs,true);
				
				lanczos_requested_sz(*ham,128,eigs,1.0,spins,true);
				cout<<"ED energies (Lanczos requested S_z=1)"<<endl;
				for (int n=0;n<eigs.size();n++)
				{cout<<eigs[n]<<"  "<<spins[n]<<endl;}
				/*print_vec(eigs,true);*/
			
				lanczos_requested_sz(*ham,128,eigs,2.0,spins,true);
				cout<<"ED energies (Lanczos requested S_z=2)"<<endl;
				for (int n=0;n<eigs.size();n++)
				{cout<<eigs[n]<<"  "<<spins[n]<<endl;}
				/*print_vec(eigs,true);*/
			}
			else
			{
				lanczos_requested_sz(*ham,256,eigs,0.0,spins,false);
				cout<<"ED energies (Lanczos requested S_z=0)"<<endl;

				print_vec(eigs,true);
				
				//cout<<"Spins "<<endl;
				//print_vec(spins,true);
				//cout<<endl;

				/*lanczos_requested_sz(*ham,256,eigs,1.0,spins,false);
				cout<<"ED energies (Lanczos requested S_z=1)"<<endl;
				print_vec(eigs,true);

				
				lanczos_requested_sz(*ham,256,eigs,2.0,spins,false);
				cout<<"ED energies (Lanczos requested S_z=2)"<<endl;
				print_vec(eigs,true);*/
				
				/*lanczos_requested_sz(*ham,1024,eigs,3.0,spins,false);
				cout<<"ED energies (Lanczos requested S_z=3)"<<endl;
				print_vec(eigs,true);*/

			}
			/*lanczos_requested_sz(*ham,512,eigs,1.0,spins,false);
			cout<<"ED energies (Lanczos requested S_z=1)"<<endl;
			print_vec(eigs,true);
			
			lanczos_requested_sz(*ham,512,eigs,-1.0,spins,false);
			cout<<"ED energies (Lanczos requested S_z=-1)"<<endl;
			print_vec(eigs,true);*/
		}


	    	if ((*ham).num_sites<=11)
		{ ed_get_eigs(*ham,eigs);
		  cout<<"ED energies (full)"<<endl;
		  print_vec(eigs,true);
		}
   	    	
	    }   
    	     
    }

/////////////////////////////////////////////////////////////////////////////
//                             Infinite DMRG (White/Wilson)
/////////////////////////////////////////////////////////////////////////////
    int                                                 max_states;
    int                                                 size;
    std::vector<double> 		                target_szs;
    std::vector<int>  			                target_numbers;
    std::vector< std::vector< Matrix > >                sp_tr,sp_hams,sp_s_p,sp_s_m,sp_s_z,sp_s_z2;
    std::vector< std::vector<int> >                     inverse_subspace_map,inverse_sz_map;
    std::vector< std::vector<double> > 		   	eigs_by_spin;
    cout<<"Ham_found="<<ham_found<<endl;	
    if (ham_found)
    {
	    search_for(string("nrg"),filename,str_ret,found);
	    if (str_to_bool(str_ret))
	    {
    	    	search_for(string("max_states"),filename,str_ret,found);
		if (found) {max_states=str_to_int(str_ret);}
		else {max_states=4;}
		cout<<"max_states = "<<max_states<<endl;
    		search_for(string("target_szs"),filename,str_ret,found);
    		if (found) if (str_ret.substr(0,1)==string("[")) target_szs=convert_string_to_vec_double(str_ret);
		else{cout<<"Could not find target s_z for DMRG. I will chose only Ground state from NRG warmup"<<endl;}    
    		search_for(string("target_numbers"),filename,str_ret,found);
    		if (found){if (str_ret.substr(0,1)==string("[")) target_numbers=convert_string_to_vec(str_ret);}
		else{cout<<"Could not find target numbers for DMRG. I will chose only Ground state from NRG warmup"<<endl;}    
		time(&start);
		nrg_spin(xxz,target_szs,target_numbers,eigs_by_spin,spins,inverse_subspace_map,inverse_sz_map,sp_tr,sp_s_p,sp_s_m,sp_s_z,sp_s_z2,sp_hams,max_states,false);
		time(&end);
		dif=difftime(end,start);

		cout<<"==================================================================="<<endl;
        	cout<<"Total time to do NRG was "<<dif<<" seconds"<<endl;
        	cout<<"==================================================================="<<endl;
	    }   
     }

/////////////////////////////////////////////////////////////////////////////
//                              CLEAN UP 
/////////////////////////////////////////////////////////////////////////////

    if (ham!=NULL) {delete ham;ham=NULL;}
    return 0;

}
//////////////////////////////////////////////////////////////////////////////
