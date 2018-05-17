#include "iostream"
#include "vector"
#include "ctime"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

#include <algorithm> 
#include <limits>
#include <set>


#define mp make_pair
#define pb push_back

// #include "omp.h"

using namespace std;

vector <vector <pair<int,int> > >  adj;
void printAdjList(vector <vector <pair<int,int> > > &adj){
    int ii,b,c=0;
    for(ii=0;ii<adj.size();ii++){
        c=0;
        for(b=0;b<adj[ii].size();b++){
            if(c==0){
                cout<<ii<<" - ";
            }
            c=1;
            cout<< adj[ii][b].first<< "(" <<adj[ii][b].second<<") ";
        }
        if(c==1){
            cout<<endl;
        } 
    }
}



vector<int> small_bisect(vector <vector <pair<int,int> > >  &adj){
    int s = 1;
    // printAdjList(adj);
    // cout<<"  "<<endl;

    int i;
    set <pair<int,int> > vertices;
    vertices.insert(make_pair(0,s)); 
    vector <int> distance(adj.size(),90000000);
    distance[s] = 0;
    pair<int,int> curr;

    // cout<<adj.size()<< " -- "<<endl;

    while(vertices.size()!=0){
        curr = *vertices.begin();
        vertices.erase(curr);

        // cout<< curr.second << " ";
        int c1,c2;
        for(i=0;i<adj[curr.second].size();i++){
            c1 =  distance[curr.second] + adj[curr.second][i].second;
            c2 = distance[adj[curr.second][i].first];
            if( c1 < c2){
                vertices.erase(make_pair(c2,adj[curr.second][i].first));
                vertices.insert(make_pair(c1,adj[curr.second][i].first));
                distance[adj[curr.second][i].first] = c1;
                // before[adj[curr.second][i].first] = curr.second;
            }
        }
    }
    distance[0] = numeric_limits<int>::max();
    vector<int> yy(distance.size());
    size_t n(0);
    // cout<< "distances are  "<<endl;
    // for (auto v : distance)
    //     std::cout << v << ' ';

    generate(begin(yy),end(yy), [&]{ return n++; });

    sort(  begin(yy), 
            end(yy),
                [&](int i1, int i2) { return distance[i1] < distance[i2]; } );
    // cout<<"y is "<<endl;
    // for (auto v : yy)
    //     std::cout << v << ' ';
    vector<int> out(adj.size(),0);
    int done = 0;
    while(done<adj.size()/2){
        // if(yy[done]!=0){
            out[yy[done]] = 1;
            done++;
        // }
    }
    out[0] = -1;    
    // zeroth element is garbage 

    //return vector starting with -1 and consisting of half zeros and half ones
    // corresponding to a cut
    
    // cout<<"out is "<<endl;
    // for (auto v : out)
    //     std::cout << v << ' ';

    return out;
}


// pair<vector<pair<int,int> > , vector<bool> > match(vector<vector<int,int> > adj)
// {
//      n = adj.size();
//      vector<bool> ism(n,false);
//      vector<int> visit(n,0);

//      for(i=0;i<n;i++)
//         visit[i] = i;

//      for(i=0;i<n;i++)
//      {
//         int m = rand()%i;
//         visit[i] = visit[m];
//         visit[m] = i;
//      }

//      vector<pair<int,int> > m;


//      for (i=0;i<n;i++)
//      {
//         int u = visit[i];
//         if(!ism[u])
//         {
//             int mx_wt = 0;
//             int mtch = -1;
//             for (int j = 0;j<adj[u].size();j++)
//             {
//                 if(!ism[adj[u][j].first])
//                 {
//                     int tmp = adj[u][j].second;
//                     if(tmp>mx_wt)
//                     {
//                         tmp = mx_wt;
//                         mtch = adj[u][j].first;
//                     }
//                 }
//             }
//             if(mtch!=-1)
//             {
//                 m.push_back(make_pair(u,mtch));
//                 ism[u] = true;
//                 ism[mtch] = true;
//             }
//         }
//      }

//     return(make_pair(m,ism));
// }



int coarsen_single_step(vector <vector <pair<int,int> > >  &adj){
    
    pair<vector<pair<int,int> > , vector<bool> >  res;

    // res = match(adj);

    vector<pair<int,int> >  matching = res.first;
    vector<bool> isMatched = res.second;
    vector<int> transformation(adj.size());
    
    int i=0,counter = 1;
    for(i=0;i<matching.size();i++){
        transformation[matching[i].first] = counter;
        transformation[matching[i].second] = counter;
        counter++;
    }
    for(i=1;i<isMatched.size();i++){
        if(!isMatched[i]){
            transformation[i] = counter;
            counter++;
        }
    }
    vector <vector <pair<int,int> > >  newadj(counter); //counter is 1 more than no of new nodes
    
    int ii,b,c=0,exists=0;
    for(ii=1;ii<adj.size();ii++){
        for(b=0;b<adj[ii].size();b++){
           // ii to ad[ii][b].first
            exists=0;
            for(c=0;c<newadj[transformation[ii]].size();c++){
                if(newadj[transformation[ii]][c].first == transformation[adj[ii][b].first] ){
                    newadj[transformation[ii]][c].second += 1;
                    exists = 1;
                }
            }
            if(exists!=1){
                newadj[transformation[ii]].pb(mp(transformation[adj[ii][b].first],1));
            }
        }
    }
    return 0;
}



vector<int> KWayCut_rec(vector <vector <pair<int,int> > >  &adj, int k,int offset){
    vector<int> bisection = small_bisect(adj);
    int i,j;
    if(k==2){
        if(offset == -1){
            return bisection;
        }
        else{
            for(i=1;i<bisection.size();i++){
                bisection[i] += 2*offset;
            }
            return bisection;
        }
    }
    vector <vector <pair<int,int> > >  adj0(adj.size());
    vector <vector <pair<int,int> > >  adj1(adj.size());
    vector<int> bisect0;
    vector<int> bisect1;
    vector<int> ret(bisection.size());
    ret[0] = -1;
    int index_0=1, index_1=1;
    if(offset==-1){
        // call with 0 and 1
        // for(i=1;i<bisection.size();i++){
        //     if(bisection[i]==0){
        //         adj0.pb(adj[i]);
        //     }
        //     else if(bisection[i]==1){
        //         adj1.pb(adj[i]);                
        //     }
        //     else{
        //         cout <<"error see k way"<<endl;
        //     }
        //}
        for(i=1;i<adj.size();i++)
        {
            for(j=0;j<adj[i].size();j++)
            {
                pair<int,int> p = adj[i][j];
                if(bisection[i]==bisection[adj[i][j].first])
                {
                    if(bisection[i]==1)
                        adj1[i].pb(p);
                    else
                        adj0[i].pb(p);
                }
            }
        }
        
        bisect0 = KWayCut_rec(adj0,k/2,0);
        bisect1 = KWayCut_rec(adj1,k/2,1);
        for(i=1;i<bisection.size();i++){
            if(bisection[i]==0){
                ret[i] = bisect0[index_0];
                index_0++;
            }
            else if(bisection[i]==1){
                ret[i] = bisect1[index_1];
                index_1++;
            }
            else{
                cout <<"error see k way"<<endl;
            }
        }

        return ret;

    }
    else{
        // call with 2*off and 2*off+1
        for(i=1;i<adj.size();i++)
        {
            for(j=0;j<adj[i].size();j++)
            {
                pair<int,int> p = adj[i][j];
                if(bisection[i]==bisection[adj[i][j].first])
                {
                    if(bisection[i]==1)
                        adj1[i].pb(p);
                    else
                        adj0[i].pb(p);
                }
            }
        }
        bisect0 = KWayCut_rec(adj0,k/2,2*offset);
        bisect1 = KWayCut_rec(adj1,k/2,2*offset+1);
        for(i=1;i<bisection.size();i++){
            if(bisection[i]==0){
                ret[i] = bisect0[index_0];
                index_0++;
            }
            else if(bisection[i]==1){
                ret[i] = bisect1[index_1];
                index_1++;
            }
            else{
                cout <<"error see k way"<<endl;
            }
        }

        return ret;
    }

}


vector<int> KWayCut(vector <vector <pair<int,int> > >  &adj, int k){
    return KWayCut_rec(adj,k,-1);
}


int main(int argc, char* argv[])
{   
    string inf = argv[1];
    string of = argv[2];
    int k = atoi(argv[3]);  
    ifstream infile(inf);
    int nodes, edges;
    infile>>nodes>>edges;
    string line;
    getline(infile,line);
 
    int i;

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
    		adj[i].push_back(mp(num,1));
    	}
    }


    //caraseining
    // splitting
    // uncoarsening
    // small_bisect(adj);
    //initialising result

    vector<int> res = KWayCut(adj,k);
    
    
    // res.resize(nodes+1);
    // for(i=1;i<=nodes;i++)
    // {
    // 	res[i] = i%k;
    // }

    // writing the result

    // coarsen(adj);
    ofstream outfile(of);

    for(i=1;i<nodes;i++)
    {
    	outfile<<(i%k)<<' ';
    }
    outfile<<(i%k)<<' ';
    outfile.close();
    return 0;
}


