#include <iostream>
#include<vector>
#include<string>
#include"number_functions.h"

using namespace std;
    
int main()
{
    std::string                         tmp_string,tmp_string_2;
    std::vector< std::vector<int> >     pairs;
    int                                 gen=7;
    std::vector<int>                    pair,active_nodes,new_active_nodes;
    int                                 nsites;

    pair.push_back(0);
    pair.push_back(1);
    pairs.push_back(pair);
    active_nodes.push_back(0);
    active_nodes.push_back(1);
    nsites=2;

    for (int i=0;i<gen;i++)
    {
	for (int j=0;j<active_nodes.size();j++)
        {
		pair.clear();
		pair.push_back(active_nodes[j]);
                pair.push_back(nsites);
                new_active_nodes.push_back(nsites);
		pairs.push_back(pair);
                nsites=nsites+1;
		pair.clear();
		pair.push_back(active_nodes[j]);
                pair.push_back(nsites);
                new_active_nodes.push_back(nsites);
		pairs.push_back(pair);
                nsites=nsites+1;
        }
        active_nodes=new_active_nodes;
        new_active_nodes.clear();
    }

    for (int l=0;l<pairs.size();l++)
    {
					//tmp_string=string("pairs_cluster_")+to_string(ctr)+string("={");    
					tmp_string+=to_string(pairs[l][0])+string("->")+to_string(pairs[l][1]);
					tmp_string_2+=string("[")+to_string(pairs[l][0])+string(",")+to_string(pairs[l][1])+string("]");
					if (l!=pairs.size()-1)
					{tmp_string+=string(",");tmp_string_2+=string(",");}
     }
     tmp_string+=string("}");
     tmp_string+=string("\n");
     tmp_string_2+=string("]");
     cout<<tmp_string_2<<endl;
     cout<<tmp_string<<endl;

     return 0;
}
