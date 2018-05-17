#include "iostream"
#include "vector"
#include "ctime"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

#include <algorithm> 
// #include <limits>
// #include <set>
#include <queue>


#define mp make_pair
#define pb push_back

// #include "omp.h"

using namespace std;


int main(int argc, char* argv[])
{   
    string inf = argv[1];
    string of = argv[2];
    int k = atoi(argv[3]);  
    ifstream infile(inf.c_str());
    int nodes, edges;
    infile>>nodes>>edges;
    string line;
    getline(infile,line);
 
    int i;
    vector<vector<int> > adj;
    // input taking

    adj.resize(nodes+1);
    for (i=1;i<=nodes;i++)
    {
    	// cout<<"********************"<<i<<"****************************\n";
    	string line;
    	getline(infile,line);
    	// cout<<line<<endl;
    	stringstream ss(line);
    	int num;
    	while(ss>>num)
    	{
    		// cout<<num<<endl;
    		adj[i].push_back(num);
    	}
    }


    //caraseining
    // splitting
    // uncoarsening
    // small_bisect(adj);
    //initialising result

    vector<int> res;
    res.resize(nodes+1);
    for(i=1;i<=nodes;i++)
    {
    	res[i] = k-1;
    }

    int size = (int)(nodes/k);
    for(i=0;i<(k-1);i++)
    {
        int count = 0;
        for(j=1;j<=nodes;j++)
        {
            if(count==size)
                break;
            if(res[j]==(k-1))
            {
                // res[j] = i;
                queue<int> q;
                q.push(j);
                while((!q.empty())&&(count<size))
                {
                    count = count + 1;
                    int x = q.front();
                    q.pop();
                    res[x] = i;
                    for (int l = 0;l<adj[x].size();l++)
                    {
                        if(res[adj[x][l]]==(k-1))
                            q.push(adj[x][l]);
                    }
                }
            }
        }
    }
    // writing the result

    // coarsen(adj);
    ofstream outfile(of.c_str());

    for(i=1;i<nodes;i++)
    {
    	outfile<<res[i]<<' ';
    }
    outfile<<res[i]<<' ';
    outfile.close();
    return 0;
}

