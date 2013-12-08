#include"get_ai.h"
#include"global.h"
#include"printing_functions.h"
#include"matrix_functions.h"
#include"matrix.h"

using namespace std;

void get_ai(std::string file_corrs,std::vector<double> &c_i_10,std::vector<double> &a_i_10)
{
   Matrix correlation_matrix,corr_eigenvecs;
   std::vector<int> ipiv;
   int info;
   int nsites=c_i_10.size();
   

   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
	  correlation_matrix(i,j)=correlation_matrix(i,j)/3.0;
	}
    }
   a_i_10.clear();
   a_i_10=c_i_10;
   corr_eigenvecs=correlation_matrix;
   dgesv(nsites,1,&*corr_eigenvecs.begin(),nsites,&*ipiv.begin(),&*a_i_10.begin(),nsites,info);

   cout<<"Correlation matrix Siz Sjz"<<endl; 
   print_real_mat(correlation_matrix);	
   
   cout<<"C_i in first and ground"<<endl;
   print_vec(c_i_10);
   
   double norm=0.0;
   double overlap=0.0;
   double orth=0.0;

   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		if (i==j) {norm=norm+(0.5*a_i_10[i]*a_i_10[i]);}
                else {norm=norm+(2.0*correlation_matrix(i,j)*a_i_10[i]*a_i_10[j]);}
	}
		overlap=overlap+(a_i_10[i]*c_i_10[i]);
   }

   for (int i=0;i<nsites;i++) a_i_10[i]=a_i_10[i]/sqrt(norm);

   cout<<"a_i in first and ground (after normalizing)"<<endl;
   print_vec(a_i_10);
   
   cout<<"Overlap numerator                              (taking a_i)  "<<overlap<<endl;
   cout<<"Norm                                           (taking a_i)  "<<norm<<endl;
   cout<<"Overlap of mode with true 1st excited state is ( taking a_i) "<<overlap/sqrt(norm)<<endl;

   double sum_ai=0.0;
 
   for (int i=0;i<nsites;i++) {sum_ai=sum_ai+a_i_10[i];}

   sum_ai=sum_ai/double(nsites);

   for (int i=0;i<nsites;i++) {a_i_10[i]=-sum_ai+a_i_10[i];}
   
   cout<<"a_i in first and ground (after subtracting)"<<endl;
   print_vec(a_i_10);

   cout<<"Checking orthogonality constraint for a_i"<<endl;
   orth=0.0;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
		orth=orth+(0.5*correlation_matrix(i,j)*a_i_10[i]);
   	}
   }

}
