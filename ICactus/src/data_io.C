#include"data_io.h"
#include"number_functions.h"

using namespace std;

template <class T>
//////////////////////////////////////////////////////////////////////////////
void dump_function_auto(std::vector<T> const &y, 
                        std::string filename)
{             
      std::string tmp_string;
      char *filechar;
      
      filechar=new char [filename.size()+1];
      strcpy (filechar,filename.c_str());
      ofstream data_file(filechar);

      if (data_file.is_open())
      {
            for (int j=0;j<y.size();j++)
            {
                tmp_string=to_string(j)+string("   ")+to_string(y[j]);
                data_file<<tmp_string;data_file<<string("\n");
            }
            data_file.close();
      }        
      else
      {     cout<<"oops! I could not dump the data for some reason! "<<endl;}

      delete[] filechar;
}

template void dump_function_auto<int>(std::vector<int> const &,
			              std::string);

template void dump_function_auto<double>(std::vector<double> const &,
					 std::string);


/////////////////////////////////////////////////////////////////////////////
template <class T>
void dump_function(std::vector<T> const &x,
                   std::vector<T> const &y, 
                   std::string filename)
{             
      std::string tmp_string;
      char *filechar;
      
      filechar=new char [filename.size()+1];
      strcpy (filechar,filename.c_str());
      ofstream data_file(filechar);

      if (data_file.is_open())
      {
            for (int j=0;j<x.size();j++)
            {
                tmp_string=to_string(x[j])+string("   ")+to_string(y[j]);
                data_file<<tmp_string;data_file<<string("\n");
            }
            data_file.close();
      }        
      else
      {     cout<<"oops! I could not dump the data for some reason! "<<endl;}

      delete[] filechar;
}

void dump_function_double(std::vector<double> const &x,
                   std::vector<double> const &y, 
                   std::string filename)
{             
      std::string tmp_string;
      char *filechar;
      
      filechar=new char [filename.size()+1];
      strcpy (filechar,filename.c_str());
      ofstream data_file(filechar);

      if (data_file.is_open())
      {
            for (int j=0;j<x.size();j++)
            {
                tmp_string=dtos(x[j])+string("   ")+dtos(y[j]);
                data_file<<tmp_string;data_file<<string("\n");
            }
            data_file.close();
      }        
      else
      {     cout<<"oops! I could not dump the data for some reason! "<<endl;}

      delete[] filechar;
}

template void dump_function<int>(std::vector<int> const &,
				 std::vector<int> const &,	
				 std::string);

template void dump_function<double>(std::vector<double> const &,
				    std::vector<double> const &,
                                    std::string);

//////////////////////////////////////////////////////////////////////////////

void dump_multiple_functions(std::vector<double> const &x,
                             std::vector< std::vector<double> > const &funcs, 
                             std::string filename)

{              
      std::string tmp_string;
      char *filechar;
      int i,j;

      filechar=new char [filename.size()+1];
      strcpy (filechar,filename.c_str());
      ofstream data_file(filechar);

      if (data_file.is_open())
      {
        for (i=0;i<funcs[0].size();i++)
        {
            tmp_string=string("");
            for (j=0;j<funcs.size();j++)
            {tmp_string+=to_string(x[j])+string("   ")+to_string(funcs[j][i])+string("   ");}
            data_file<<tmp_string;data_file<<string("\n");
        }  
        data_file.close();
      }        
      else
      {cout<<"oops! I could not dump the data for some reason! "<<endl;}
      delete[] filechar;
}

/////////////////////////////////////////////////////////////////////////////

