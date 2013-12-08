#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////////
void calc_hints_sxsx_sysy(double coupling, 
                         int first, 
                         int second, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        // Apply S+ S- 
        if (config[first] != config[second])
        {
            touched_sites.push_back(first);
            touched_sites.push_back(second);
            
            vals_on_touched.push_back(config[second]);
            vals_on_touched.push_back(config[first]);

            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);
            
            hints_list.push_back(coupling*0.5);
        }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_szsz(double coupling, 
                     int first, 
                     int second, 
                     std::vector<int> const &config,
                     std::vector< std::vector<int> > &touched_sites_list,
                     std::vector< std::vector<int> > &vals_on_touched_list,
                     std::vector< complex<double> > &hints_list)

{
            std::vector<int> touched_sites, vals_on_touched;
            
            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);
            hints_list.push_back(coupling*(double(config[first])-0.5)*(double(config[second])-0.5));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void calc_hints_szsz_all(double coupling, 
                         std::vector< std::vector<int> > const &pairs, 
                         std::vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)

{
            //szsz is diagonal in sz basis
            int first,second;
            double hint=0.0;
            std::vector< std::vector<int> >:: const_iterator p;
            std::vector<int> touched_sites, vals_on_touched;
            
            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);

            for (p=pairs.begin();p!=pairs.end();p++)
            {
                first=(*p)[0];second=(*p)[1];
                hint+=((double(config[first])-0.5)*(double(config[second])-0.5));
            }
            hints_list.push_back( complex<double> (coupling*hint) );
}
     
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_stag_sz(double hstag, 
                         std::vector<int> const &eta, 
                         std::vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
            //szsz is diagonal in sz basis
            int i;
            double hint=0.0;
            std::vector<int> touched_sites, vals_on_touched;
            
            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);

            for (i=0;i<eta.size();i++)
            {
                hint+=((double(config[i])-0.5)*(double(eta[i])));
            }
            hints_list.push_back( complex<double> (hstag*hint) );
}
     

///////////////////////////////////////////////////////////////////////////////

void calc_hints_sx(double coupling, 
                   int site, 
                   std:: vector<int> const &config,
                   std::vector< std::vector<int> > &touched_sites_list,
                   std::vector< std::vector<int> > &vals_on_touched_list,
                   std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        // Apply S+ + S - / 2 --> Flip a spin! 
        
        touched_sites.push_back(site);
        
        vals_on_touched.push_back((config[site]+1)%2); //Flipped spin half! 

        touched_sites_list.push_back(touched_sites);
        vals_on_touched_list.push_back(vals_on_touched);
        
        hints_list.push_back(coupling*0.5);
}
///////////////////////////////////////////////////////////////////////////////

void compute_c_plus_minus_i(int plus_or_minus,
			    std::vector<double> const &vec_0,
			    std::vector<double> const &vec_1,
   			    std::vector<int> const &maps_0,
   			    std::vector<int> const &inverse_map,
			    std::vector<double> &c_i)
{
      int i,j;
      int bra_config,ket_config;
      int loc;
      int num_sites=c_i.size();
      std::vector<int> config;

      //overlap=0.0;
      for (i=0;i<num_sites;i++) c_i[i]=0.0;
      
      //double norm=0.0;
      // check norms
      //for (int i=0;i<vec_0.size();i++) norm=norm+vec_0[i]*vec_0[i];
      //cout<<"Norm of 0 state"<<norm<<endl;
     
      //norm=0.0; 
      //for (int i=0;i<vec_1.size();i++) norm=norm+vec_1[i]*vec_1[i];
      //cout<<"Norm of 1 state"<<norm<<endl;
      
      if (plus_or_minus==1)
      {
              //cout<<"c_i for S=1 S_z = +1  S=0 S_z=0   FORMULA c_i=<exc|S_i_+|GS>"<<endl;
              for (i=0;i<num_sites;i++)
	      {
		   for (j=0;j<maps_0.size();j++)
		   {
			ket_config=maps_0[j];
			convert_num_to_vec(ket_config,2,num_sites,config);
			if (config[i]==0)
			{
				bra_config=ket_config+pow(2,num_sites-1-i);
				loc=inverse_map[bra_config];
				c_i[i]+=(vec_1[loc]*vec_0[j]);
			}   
		   }
	      }
      }
      
      if (plus_or_minus==-1)
      {
              //cout<<"c_i for S=1 S_z = -1  S=0 S_z=0   FORMULA c_i=<exc|S_i_-|GS>"<<endl;
              for (i=0;i<num_sites;i++)
	      {
		   for (j=0;j<maps_0.size();j++)
		   {
			ket_config=maps_0[j];
			convert_num_to_vec(ket_config,2,num_sites,config);
			if (config[i]==1)
			{
				bra_config=ket_config-pow(2,num_sites-1-i);
				loc=inverse_map[bra_config];
				c_i[i]+=(vec_1[loc]*vec_0[j]);
			}   
		   }
	      }
      }
      
      if (plus_or_minus==0)
      {
             //cout<<"c_i for S=1 S_z =  0  S=0 S_z=0   FORMULA c_i=<exc|S_i_z|GS>"<<endl;
             for (i=0;i<num_sites;i++)
	     {
		   for (j=0;j<maps_0.size();j++)
		   {
			ket_config=maps_0[j];
			convert_num_to_vec(ket_config,2,num_sites,config);
			c_i[i]+=(vec_0[j]*vec_1[j]*(double(config[i])-0.5));
			//overlap+=(vec_0[j]*vec_1[j]);
		   }   
	     }
	   
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_sj(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=m;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    
                            if (m!=n)
                            {
				    if (config[m]==0 and config[n]==1)
				    {
					bra_config=ket_config+pow(2,num_sites-1-m);
					bra_config=bra_config-pow(2,num_sites-1-n);
					loc=inverse_map[bra_config];
					si_sj(m,n)+=(0.5*(vec_1[loc]*vec_0[j]));
				    }

				    if (config[n]==0 and config[m]==1)
				    {
					bra_config=ket_config-pow(2,num_sites-1-m);
					bra_config=bra_config+pow(2,num_sites-1-n);
					loc=inverse_map[bra_config];
					si_sj(m,n)+=(0.5*(vec_1[loc]*vec_0[j]));
				    }
			    }
                            else
                            {si_sj(m,n)+=(0.5*(vec_1[j]*vec_0[j]));}

			    si_sj(m,n)+=(vec_1[j]*vec_0[j]*(double(config[m])-0.5)*(double(config[n])-0.5));
			    si_sj(n,m)=si_sj(m,n);
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_plus_sj_minus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    if ( (config[m]==0 and config[n]==1) and (m!=n))
			    {
				bra_config=ket_config+pow(2,num_sites-1-m);
				bra_config=bra_config-pow(2,num_sites-1-n);
				loc=inverse_map[bra_config];
				si_sj(m,n)+=((vec_1[loc]*vec_0[j]));
			    }
                            if (m==n and config[m]==1)
		            {
				si_sj(m,n)+=((vec_1[j]*vec_0[j]));
                            }
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_minus_sj_minus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    if ( (config[m]==1 and config[n]==1) and (m!=n))
			    {
				bra_config=ket_config-pow(2,num_sites-1-m);
				bra_config=bra_config-pow(2,num_sites-1-n);
				loc=inverse_map[bra_config];
				si_sj(m,n)+=((vec_1[loc]*vec_0[j]));
			    }
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_minus_sj_plus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    if ( (config[m]==1 and config[n]==0) and (m!=n))
			    {
				bra_config=ket_config-pow(2,num_sites-1-m);
				bra_config=bra_config+pow(2,num_sites-1-n);
				loc=inverse_map[bra_config];
				si_sj(m,n)+=((vec_1[loc]*vec_0[j]));
			    }
                            if (m==n and config[m]==0)
		            {
				si_sj(m,n)+=((vec_1[j]*vec_0[j]));
                            }
			}
	      }
      }
}
////////////////////////////////////////////////////////////////////////////////
void compute_si_minus_sj_plus_sk_z(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   std::vector<double> &three_pt)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      three_pt.clear();
      three_pt.resize(num_sites*num_sites*num_sites);
      three_pt.assign(num_sites*num_sites*num_sites,0.0);
 
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int l=0;l<num_sites;l++)
		       {
			       for (int j=0;j<vec_0.size();j++)
			       {
				    ket_config=maps_0[j];
				    convert_num_to_vec(ket_config,2,num_sites,config);
				    if ((config[m]==1 and config[n]==0) and (m!=n))
				    {
					bra_config=ket_config-pow(2,num_sites-1-m);
					bra_config=bra_config+pow(2,num_sites-1-n);
					loc=inverse_map[bra_config];
					three_pt[num_sites*num_sites*m+num_sites*n+l]+=((vec_1[loc]*vec_0[j]*(double(config[l])-0.5)));
				    }
				    if (m==n and config[m]==0)
				    {
					three_pt[num_sites*num_sites*m+num_sites*n+l]+=((vec_1[j]*vec_0[j])*(double(config[l])-0.5));
				    }
				}
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_z_sj_z(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
		   Matrix &si_sj)
{
      int 	        ket_config;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    si_sj(m,n)+=((vec_1[j]*vec_0[j]*(double(config[m])-0.5)*(double(config[n])-0.5)));
		       }
	      }
      }
}

