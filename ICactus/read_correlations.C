#include <fstream>
#include <iostream>
#include <string>
#include"./src/matrix.h"
#include<stdio.h>
#include<fstream>
#include<sstream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<wchar.h>
#include<cwchar>
#include<boost/lexical_cast.hpp>
#include<boost/function.hpp>
#include<boost/format.hpp>
#include<complex>
#include<iomanip>
#include<ctime>
#include"./bethe_lapack_interface.h"

using namespace std;

void split(string &str, string sep_char, string &str_bef, string &str_aft)
{
    // Separate the string into two at the specified separating character
    size_t found;
    size_t length;
    int loc;

    char buffer[5000];
    // I am assuming the user wont put an entry greater than 100 charaters long! 
    // (50 was kinda not enough -given long name of files)

    found=str.find(sep_char);
    if (found!=string::npos)
    {
        loc=int(found);
        
        length=str.copy(buffer,loc,0);
        buffer[length]='\0';
        str_bef=buffer;
        
        length=str.copy(buffer,str.size()-loc,loc+1);
        buffer[length]='\0';
        str_aft=buffer;
    }
}

void search_for(string str_in,string filename, string &str_ret, bool &found)
{
    // The input file has been formatted in a certain way 
    // We assume this form - For eg.
    // //////////////////////////////////////////////
    // hamiltonian=heisenberg
    // j_x=1.0
    // j_z=0.5
    ///...
    // //////////////////////////////////////////////
    // i.e. Format is
    // string_in=string_ret
    
    str_ret="";
    string str_read,str_bef,str_aft;
    
    found=false;
    int line=0;

    // Open file
    //cout<<"Defining ifstream"<<endl;
    ifstream fl(filename.c_str(),std::ios::in);
    //cout<<"Defined ifstream"<<endl;
    
    if (!fl)
    {   cerr<<"File "<<filename<<" could not be opened"<<endl;
        exit(1);
    }
    
    if (fl.is_open())
    {
        // Begin to read lines
        while ((! fl.eof()) and (not found))
        {
            getline(fl,str_read);
            line++;
            //cout<<"Line number I just read = "<<line<<endl;
            //cout<<"String I just read ="<<str_read<<endl;
            // Take input string and separator at "=" sign
            split(str_read,string("="),str_bef,str_aft);
            if (str_bef.compare(str_in)==0)
            {   
                fl.close();
                str_ret=str_aft;
                found=true;
            }
        }
        //if (not found) {cout<<" File "<<filename<<" was opened successfully but requested string ' "<< str_in<<" ' not found "<<endl;}
    }
}

template <class T>
std::string to_string (T const &t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}
template std::string to_string<int>(int const &);
template std::string  to_string<double>(double const &);

int str_to_int(string str)
{ return boost::lexical_cast<int>(str);}

double str_to_d(string str)
{ return boost::lexical_cast<long double>(str);}

std::vector<int> convert_string_to_vec(std::string const &s)
{
    // Eg s=[1,0] would get converted to vec[1,0]
    std::vector<int> v;
    int num;
    
    std::string temp="";
    for (int i=1;i<s.size();i++)
    {
        if (s.substr(i,1)!=string(",") and s.substr(i,1)!=string("]") 
            and s.substr(i,1)!=string(")")
            and s.substr(i,1)!=string("}"))
        {   temp+=s.substr(i,1);}    
        else
        {   num=str_to_int(temp);
            v.push_back(num);
            temp="";}
    }
    return v;
}

std::vector< std::vector<int> > convert_string_to_vec_of_vec(std::string const &s)
{
    // Eg s=[[1,0],[1,1]] would get converted to vec[1,0]+vec[1,1]
    std::vector<int> temp_vec;
    std::vector< vector<int> > v;

    string temp="";
    for (int i=1;i<s.size();i++)
    {
        if (s.substr(i,1)!=string("]") and i!=s.size())
        {   temp+=s.substr(i,1);
        }    
        else
        {
            temp+=s.substr(i,1);
            temp_vec=convert_string_to_vec(temp);
            v.push_back(temp_vec);
            temp="";
            i++;
        }
    }
    return v;
}

void print_vec_acc(std::vector<double> &vec,bool vert,int size)
{
    if (size==0) {size=vec.size();}
    
    size=min(size,int(vec.size()));

    if (not vert)
    {		
    	for (int i=0;i<size;i++){cout<<boost::format("%+.10f") %vec[i]<<"  ";}
    	cout<<endl;
    }
    else
    {   for (int i=0;i<size;i++){cout<<boost::format("%+.10f") %vec[i]<<endl;} }
}

void print_real_mat(Matrix &mat)
{
    int i,j;
    for (i=0;i<mat.NRows();i++)
    {
        for (j=0;j<mat.NCols();j++) cout<<boost::format("%+.2f") % mat(i,j)<<" ";
        cout<<endl;
    }
    cout<<endl;
}

void read_correlations_and_ci(int size,Matrix &correlations,std::vector<double> &ci)
{
   int ind_1,ind_2;
   string name_s,name2_s;
   
   correlations.resize(size,size);
   ci.resize(size);

   name_s="./"+to_string(size)+"pure/"+to_string(size)+"_site_corr.txt";
   char *name;
   name=new char [name_s.size()+1];
   strcpy (name,name_s.c_str());
  
   ifstream infilecorr(name);
   while (infilecorr)
   {
	    string s;
	    if (!getline( infilecorr, s )) break;

	    istringstream ss( s );
	    vector <string> record;

	    while (ss)
	    {
	      string s;
	      if (!getline( ss, s, ',' )) break;
	      record.push_back( s );
	    }
            ind_1=str_to_int(record[0]);
            ind_2=str_to_int(record[1]);
            correlations(ind_1,ind_2)=str_to_d(record[2]);
            correlations(ind_2,ind_1)=correlations(ind_1,ind_2);
    }
   if (!infilecorr.eof()){cerr << "Fooey!\n";}

   //print_real_mat(correlations);
   
   name2_s="./"+to_string(size)+"pure/"+to_string(size)+"_ci.txt";
   char *name2;
   name2=new char [name2_s.size()+1];
   strcpy (name2,name2_s.c_str());
   ifstream infileci( name2 );
   ind_1=0;
   while (infileci and ind_1<size)
   {
	    string s;
            infileci>>s; 
            infileci>>s;
            ci[ind_1]=str_to_d(s);
            ind_1=ind_1+1;
    }
   //print_vec_acc(ci,true,size);

}

void quantities_of_interest (std::vector< std::vector<int> > &pairs,Matrix &correlations, std::vector<double> &ci, double &gap_num,double &gap_denom,double &sma_gap,double &overlap)
{
   std::vector<double> ui;
   Matrix corr_eigenvecs,correlation_matrix;
   int nsites;
   double norm,orth;
   nsites=ci.size();
   std::vector<int> ipiv;
   ipiv.resize(nsites);
   int info;

   correlation_matrix=correlations;
   corr_eigenvecs=correlations;
   ui=ci;
   for (int i=0;i<correlations.NRows()*correlations.NCols();i++) correlation_matrix[i]=(1.0/3.0)*correlations[i];
   for (int i=0;i<correlations.NRows()*correlations.NCols();i++) corr_eigenvecs[i]=(4.0/3.0)*correlations[i];
   dgesv(nsites,1,&*corr_eigenvecs.begin(),nsites,&*ipiv.begin(),&*ui.begin(),nsites,info);

   //cout<<"Correlation matrix Siz Sjz"<<endl; 
   //print_real_mat(correlation_matrix);	
   
   //cout<<"C_i in first and ground"<<endl;
   //print_vec_acc(ci,true,nsites);
   
   norm=0.0;overlap=0.0;orth=0.0;

  cout<<"=============================================================================================================="<<endl;
  cout <<"                                  Optimal ui ansatz for "<<nsites<<" sites"<<endl;
  cout<<"=============================================================================================================="<<endl;
   
  for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		if (i==j) {norm=norm+(0.5*ui[i]*ui[i]);}
                else {norm=norm+(2.0*correlation_matrix(i,j)*ui[i]*ui[j]);}
	}
		overlap=overlap+(ui[i]*ci[i]);
   }

   for (int i=0;i<nsites;i++) ui[i]=ui[i]/sqrt(norm);

   //cout<<"ui in first and ground (after normalizing)"<<endl;
   //print_vec_acc(ui,true,nsites);
   
   cout<<"Overlap numerator                              (taking optimal ui)  "<<overlap<<" for "<<nsites<<" sites"<<endl;
   cout<<"Norm                                           (taking optimal ui)  "<<norm<<" for "<<nsites<<" sites"<<endl;
   cout<<"Overlap of mode with true 1st excited state is (taking optimal ui)  "<<overlap/sqrt(norm)<<" for "<<nsites<<" sites"<<endl;

   double sum_ai=0.0;
 
   for (int i=0;i<nsites;i++) {sum_ai=sum_ai+ui[i];}

   sum_ai=sum_ai/double(nsites);

   for (int i=0;i<nsites;i++) {ui[i]=-sum_ai+ui[i];}
   
   //cout<<"ui in first and ground (after subtracting)"<<endl;
   //print_vec_acc(ui,true,nsites);

   //cout<<"Checking orthogonality constraint for ui"<<endl;
   orth=0.0;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		orth=orth+(0.5*correlation_matrix(i,j)*ui[i]);
   	}
   }
   //cout<<"Orth ="<<orth<<endl;
 
   gap_num=0.0;
   gap_denom=0.0;
   for (int i=0;i<ui.size();i++) ui[i]=ui[i]/abs(ui[nsites-1]);

   for (int i=0;i<pairs.size();i++)
   {
	gap_num+=(pow((ui[pairs[i][0]]-ui[pairs[i][1]]),2.0)*(-2.0)*(correlation_matrix(pairs[i][0],pairs[i][1])));
	//gap_num+=(ui[pairs[i][0]]*ui[pairs[i][1]]*(-2.0/3.0)*(correlations(pairs[i][0],pairs[i][1])));
   }
   cout<<"Gap numerator                                  (taking optimal ui)  "<<gap_num<<" for "<<nsites<<" sites"<<endl;

   for (int i=0;i<nsites;i++)
   {
 	for (int j=0;j<nsites;j++)
	{
		gap_denom+=(ui[i]*ui[j]*(2.0)*correlation_matrix(i,j));
	}
   }
   cout<<"Gap denominator                                (taking optimal ui) "<<gap_denom<<" for "<<nsites<<" sites"<<endl;
   cout<<"SMA Gap                                        (taking optimal ui) "<<gap_num/gap_denom<<" for "<<nsites<<" sites"<<endl;
   cout<<endl;
 
  // Try the LR ansatz
  cout<<"=============================================================================================================="<<endl;
  cout <<"                                  LR ansatz for "<<nsites<<" sites"<<endl;
  cout<<"=============================================================================================================="<<endl;

  for (int i=0;i<nsites;i++) ui[i]=ui[i]/abs(ui[i]);
 
  //print_vec_acc(ui,true,nsites);

  norm=0.0;overlap=0.0;orth=0.0;
  cout<<endl;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		if (i==j) {norm=norm+(0.5*ui[i]*ui[i]);}
                else {norm=norm+(2.0*correlation_matrix(i,j)*ui[i]*ui[j]);}
	}
		overlap=overlap+(ui[i]*ci[i]);
   }

   cout<<"Overlap numerator                              (taking LR)         "<<overlap<<" for "<<nsites<<" sites"<<endl;
   cout<<"Norm                                           (taking LR)         "<<norm<<" for "<<nsites<<" sites"<<endl;
   cout<<"Overlap of mode with true 1st excited state is (taking LR)         "<<overlap/sqrt(norm)<<" for "<<nsites<<" sites"<<endl;

   //cout<<"Checking orthogonality constraint for ui"<<endl;
   orth=0.0;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		orth=orth+(0.5*correlation_matrix(i,j)*ui[i]);
   	}
   }
   //cout<<"Orth ="<<orth<<endl;
 
   gap_num=0.0;
   gap_denom=0.0;

   for (int i=0;i<pairs.size();i++) gap_num+=(pow((ui[pairs[i][0]]-ui[pairs[i][1]]),2.0)*(-2.0)*(correlation_matrix(pairs[i][0],pairs[i][1])));
   cout<<"Gap numerator                                  (taking LR)        "<<gap_num<<" for "<<nsites<<" sites"<<endl;

   for (int i=0;i<nsites;i++)
   {
 	for (int j=0;j<nsites;j++)
	{
		gap_denom+=(ui[i]*ui[j]*(2.0)*correlation_matrix(i,j));
	}
   }
   cout<<"Gap denominator                                (taking LR)        "<<gap_denom<<" for "<<nsites<<" sites"<<endl;
   cout<<"SMA Gap                                        (taking LR)        "<<gap_num/gap_denom<<" for "<<nsites<<" sites"<<endl;
   cout<<endl;

  // Try the staggered ansatz
  cout<<"=============================================================================================================="<<endl;
  cout <<"                                  Staggered ansatz for "<<nsites<<" sites"<<endl;
  cout<<"=============================================================================================================="<<endl;

  for (int i=0;i<nsites;i++) ui[i]=ci[i]/abs(ci[i]);
 
  //print_vec_acc(ui,true,nsites);

  norm=0.0;overlap=0.0;orth=0.0;
  cout<<endl;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		if (i==j) {norm=norm+(0.5*ui[i]*ui[i]);}
                else {norm=norm+(2.0*correlation_matrix(i,j)*ui[i]*ui[j]);}
	}
		overlap=overlap+(ui[i]*ci[i]);
   }

   cout<<"Overlap numerator                              (taking stag)         "<<overlap<<" for "<<nsites<<" sites"<<endl;
   cout<<"Norm                                           (taking stag)         "<<norm<<" for "<<nsites<<" sites"<<endl;
   cout<<"Overlap of mode with true 1st excited state is (taking stag)         "<<overlap/sqrt(norm)<<" for "<<nsites<<" sites"<<endl;

   //cout<<"Checking orthogonality constraint for ui"<<endl;
   orth=0.0;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		orth=orth+(0.5*correlation_matrix(i,j)*ui[i]);
   	}
   }
   //cout<<"Orth ="<<orth<<endl;
 
   gap_num=0.0;
   gap_denom=0.0;

   for (int i=0;i<pairs.size();i++) gap_num+=(pow((ui[pairs[i][0]]-ui[pairs[i][1]]),2.0)*(-2.0)*(correlation_matrix(pairs[i][0],pairs[i][1])));
   cout<<"Gap numerator                                  (taking stag)        "<<gap_num<<" for "<<nsites<<" sites"<<endl;

   for (int i=0;i<nsites;i++)
   {
 	for (int j=0;j<nsites;j++)
	{
		gap_denom+=(ui[i]*ui[j]*(2.0)*correlation_matrix(i,j));
	}
   }
   cout<<"Gap denominator                                (taking stag)        "<<gap_denom<<" for "<<nsites<<" sites"<<endl;
   cout<<"SMA Gap                                        (taking stag)        "<<gap_num/gap_denom<<" for "<<nsites<<" sites"<<endl;
   cout<<endl;

}



int main()
{
   Matrix correlations_30(30,30),correlations_62(62,62),correlations_126(126,126);
   std::vector<double> c_30(30),c_62(62),c_126(126);
   std::vector< vector<int> > pairs_30,pairs_62,pairs_126;
   int nsites;
   double gap_num,gap_denom,sma_gap,overlap;
   string str_ret;
   bool found;

   search_for(string("pairs_30"),string("./pure_bethe_adj.txt"),str_ret,found);
   pairs_30=convert_string_to_vec_of_vec(str_ret);
   search_for(string("pairs_62"),string("./pure_bethe_adj.txt"),str_ret,found);
   pairs_62=convert_string_to_vec_of_vec(str_ret);
   search_for(string("pairs_126"),string("./pure_bethe_adj.txt"),str_ret,found);
   pairs_126=convert_string_to_vec_of_vec(str_ret);

   read_correlations_and_ci(30,correlations_30,c_30);
   read_correlations_and_ci(62,correlations_62,c_62);
   read_correlations_and_ci(126,correlations_126,c_126);
   quantities_of_interest (pairs_30,correlations_30, c_30, gap_num,gap_denom,sma_gap,overlap);
   quantities_of_interest (pairs_62,correlations_62, c_62, gap_num,gap_denom,sma_gap,overlap);
   quantities_of_interest (pairs_126,correlations_126, c_126, gap_num,gap_denom,sma_gap,overlap);
   return 0;
}
