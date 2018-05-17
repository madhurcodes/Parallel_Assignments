#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>        
#include <algorithm> 
#include <limits>
#include <set>
#include <cmath>


#define mp make_pair
#define pb push_back

// #include "omp.h"

using namespace std;

// vector <vector <pair<int,int> > >  adj;
void printAdjList(vector <vector <pair<int,int> > > &adj){
    int ii,b,c=0;
    for(ii=0;ii<adj.size();ii++){
        c=0;
        for(b=0;b<adj[ii].size();b++){
            if(c==0){
                // cout<<ii<<" - ";
            }
            c=1;
            // cout<< adj[ii][b].first<< "(" <<adj[ii][b].second<<") ";
        }
        if(c==1){
            // cout<<endl;
        } 
    }
}

vector<int> iter_refine1(vector<int> lbl,vector<vector<pair<int,int> > > &adj)
{
    int z_cnt = 0;
    int o_cnt = 0;

    int i,n;

    n = lbl.size();

    vector<int> d(n,0);

    for(i=1;i<adj.size();i++)
    {
        for(int j=0;j<adj[i].size();j++)
        {
            int v = adj[i][j].first;
            int x = adj[i][j].second;
            if(lbl[i]==lbl[v])
            {
                d[i] = d[i] - x;
                d[v] = d[v] - x;
            }
            else
            {
                d[i] = d[i] + x;
                d[v] = d[v] + x;
            }
        }
    }
    // cout<<"d done one"<<endl;
    vector <pair<int,int> > d_0;
    vector <pair<int,int> > d_1;
    for(i=1;i<n;i++)
    {
        if(lbl[i]==0)
        {
            z_cnt++ ;
            d_0.pb(mp(-1*d[i],i));
        }
        else
        {
            o_cnt++ ;
            d_1.pb(mp(-1*d[i],i));
        }
    }
    if(z_cnt>o_cnt)
    {
        sort(d_0.begin(),d_0.end());
        int diff = int(0.95*(z_cnt - o_cnt)/2);
        for(i=0;i<diff;i++)
            lbl[d_0[i].second] = 1;
    }
    else if(z_cnt < o_cnt)
    {
        sort(d_1.begin(),d_1.end());
        int diff = int(0.95*(o_cnt - z_cnt)/2);
        for(i=0;i<diff;i++)
            lbl[d_1[i].second] = 0;
    }
    // cout<<"size equalization done"<<endl;
    int cnt = 0;
    int gmax = 1;
    while((gmax>0)&&(cnt<1)&&false)
    {
        cnt = cnt + 1;
        vector<int> d(n,0);
        // cout<<"************"<<cnt<<"**************"<<endl;
        for(i=1;i<adj.size();i++)
        {
            for(int j=0;j<adj[i].size();j++)
            {
                int v = adj[i][j].first;
                int x = adj[i][j].second;
                if(lbl[i]==lbl[v])
                {
                    d[i] = d[i] - x;
                    d[v] = d[v] - x;
                }
                else
                {
                    d[i] = d[i] + x;
                    d[v] = d[v] + x;
                }
            }
        }
        // cout<<"d done"<<endl;
        vector<int> g,a,b;
        // int limit = (int)(log(n/2)/log(2));
        int limit = 0;
        // for(i=0;i<(n/2);i++)
        for(i=0;i<limit;i++)
        {
            int temp_a = 0;
            int temp_b = 0;
            int temp_g = 0;

            for(int a = 1;a<n;a++)
            {
                for(int b = a+1;b<n;b++)
                {
                    if(lbl[a]!=lbl[b])
                    {
                        int cst = d[a] + d[b];
                        // for(int j = 0;j<adj[a].size();j++)
                        // {
                        //     if(adj[a][j].first==b)
                        //     {
                        //         cst = cst - 2*adj[a][j].second;
                        //     }
                        // }
                        if(cst>temp_g)
                        {
                            temp_g = cst;
                            if(lbl[a]==0)
                            {
                                temp_a = a;
                                temp_b = b;
                            }
                            else
                            {
                                temp_a = b;
                                temp_b = a;
                            }
                        }
                    }
                }
            }

            if(temp_a!=0)
            {
                lbl[temp_a] = -1;
                lbl[temp_b] = -1;
                g.pb(temp_g);
                a.pb(temp_a);
                b.pb(temp_b);
                
            }
        }

        for(i=0;i<a.size();i++)
        {
            lbl[a[i]] = 0;
            lbl[b[i]] = 1;
        }

        gmax = 0;
        int kmax = -1;
        int cursum = 0;

        for (int k1 = 0;k1<g.size();k1++)
        {
            cursum = cursum + g[k1];
            if(cursum>gmax)
            {
                kmax = k1;
                gmax = cursum;
            }
        }

        for(i=0;i<=kmax;i++)
        {
            lbl[a[i]] = 1;
            lbl[b[i]] = 0;
        }
    }
    // cout<<"exiting iter refine"<<endl;
    return lbl;
}






vector<int> iter_refine(vector<int> lbl,vector<vector<pair<int,int> > > &adj)
{
    int z_cnt = 0;
    int o_cnt = 0;

    int i,n;

    n = lbl.size();

    for(i=1;i<n;i++)
    {
        if(lbl[i]==0)
        {
            z_cnt++ ;
            // d_0.pb(mp(-1*d[i],i));
        }
        else
        {
            o_cnt++ ;
            // d_1.pb(mp(-1*d[i],i));
        }
    }
    int cnt = 0;
    int diff = abs(z_cnt - o_cnt);
    while(cnt<60)
    {
        vector<int> d(n,0);

        for(i=1;i<adj.size();i++)
        {
            for(int j=0;j<adj[i].size();j++)
            {
                int v = adj[i][j].first;
                int x = adj[i][j].second;
                if(lbl[i]==lbl[v])
                {
                    d[i] = d[i] - x;
                    // d[v] = d[v] - x;
                }
                else
                {
                    d[i] = d[i] + x;
                    // d[v] = d[v] + x;
                }
            }
        }

        int mx_cng = 0;
        int cng_lbl = -1;
        int cng_val;
        if(z_cnt>o_cnt)
        {
            cng_val = 1;
            for(i=1;i<n;i++)
            {
                if(lbl[i]==0)
                {
                    if(d[i]>mx_cng)
                    {
                        mx_cng = d[i];
                        cng_lbl = i;
                    }
                }
            }
        }
        else if(z_cnt<o_cnt)
        {
            cng_val = 0;
            for(i=1;i<n;i++)
            {
                if(lbl[i]==1)
                {
                    if(d[i]>mx_cng)
                    {
                        mx_cng = d[i];
                        cng_lbl = i;
                    }
                }
            }
        }
        // else{
        //     break;
        // }
        if(mx_cng <= 0){
            cnt++;
            // break;
        }
        else{
            lbl[cng_lbl] = cng_val;
            cnt = 0;
        }

        diff = abs(z_cnt -o_cnt);
    }

    return lbl;
}




vector<int> small_bisect(vector <vector <pair<int,int> > >  &adj){
    // cout<<"**********inside small bisect*************"<<endl;
    int s = 1;
    // printAdjList(adj);
    // cout<<"  "<<endl;

    int i;
    set <pair<int,int> > vertices;
    vertices.insert(make_pair(0,s)); 
    vector <int> distance(adj.size(),numeric_limits<int>::max());
    // vector <int> before(adj.size(), -1);
    // before[s] = s;
    distance[s] = 0;
    pair<int,int> curr;
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
    // cout<<adj.size()<<endl;
    while(done<adj.size()/2){
        // cout<<"hello world"<<endl;
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



pair<vector<pair<int,int> >,vector<bool> > match(vector<vector<pair<int,int> > > &adj)
{
     int i;
     int n = adj.size();
     vector<bool> ism(n,false);
     vector<int> visit(n,0);
     ism[0] = true;
     
     for(i=0;i<n;i++)
        visit[i] = i;
     
    //  std::srand ( unsigned ( std::time(0) ) );
    //  random_shuffle(visit.begin()+1,visit.end());
     // cout<<n<<endl;
     // for(i=1;i<n;i++)
     // {
     //    int x = rand();
     //    // cout<<x<<endl;
     //    int m = (x)%i;
     //    // cout<<i<<" "<<m<<endl;
     //    visit[i] = visit[m];
     //    visit[m] = i;
     // }

     vector<pair<int,int> > m;

     for (i=1;i<n;i++)
     {
        int u = visit[i];
        if(!ism[u])
        {
            int mx_wt = 1000000;
            int mtch = -1;
            for (int j = 0;j<adj[u].size();j++)
            {
                if(!ism[adj[u][j].first])
                {
                    int tmp = adj[u][j].second;
                    if(tmp<mx_wt)
                    {
                        mx_wt = tmp;
                        mtch = adj[u][j].first;
                    }
                }
            }
            if(mtch!=-1)
            {
                m.push_back(make_pair(u,mtch));
                ism[u] = true;
                ism[mtch] = true;
            }
        }
     }
    // vector <vector <pair<int,int> > >  adj1;
    // adj1.pb(m);
    // printAdjList(adj1);
    return(make_pair(m,ism));
}


pair<vector<vector<pair<int,int> > > , vector<int> > coarsen_single_step(vector <vector <pair<int,int> > >  &adj){
    
    pair<vector<pair<int,int> > , vector<bool> >  res;

    res = match(adj);

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
    pair<vector<vector<pair<int,int> > > , vector<int> > ret1;
    ret1 = mp(newadj,transformation);
    return ret1;
}

vector<int> fullbisect(vector <vector <pair<int,int> > >  &adj)
{
    pair<vector<vector<pair<int,int> > > , vector<int> > ret1 ;

    ret1 = coarsen_single_step(adj);
    int step = 0;
    vector <  pair<vector<vector<pair<int,int> > > , vector<int> > > store;
    store.pb(ret1);

    while(ret1.first.size()>200 && (step<50))
    {
        // cout<<step<<" "<<ret1.first.size()<<endl;
        step++;
        ret1 = coarsen_single_step(ret1.first);
        store.pb(ret1);
    }

    vector<int> small_bisect_ret = small_bisect(ret1.first);
    // cout<<"small bisect done\n";
    // small_bisect_ret = iter_refine(small_bisect_ret,ret1.first);

    vector<int> label ; 
    vector<int> temp;
    temp = small_bisect_ret;

    for (int i=store.size()-1;i>=0;i--)
    {
        // cout<<"********"<<i<<"*****************"<<endl;
        label.resize(0);
        label.resize(store[i].second.size());
        label[0] = -1;
        for(int j = 1;j<label.size();j++)
        {
            label[j] = temp[store[i].second[j]];
        }
        // cout<<"beginning iter refine "<<endl;
        if(i>0)
            label = iter_refine(label,store[i-1].first);
        temp = label;
    }
    label = iter_refine(label,adj);
    label = iter_refine1(label,adj);    
    return label;
}


vector<int> KWayCut_rec(vector <vector <pair<int,int> > >  &adj, int k,int offset){
    // cout<<"**********inside KWayCut_rec **************"<<offset<<"***************"<<endl;
    vector<int> bisection = fullbisect(adj);

    // bisection = iter_refine(bisection,adj);
    int i,j;
    if(k==2){
        for(i=1;i<bisection.size();i++){
                bisection[i] += 2*offset;
            }
        return bisection;
    }
    vector<int> newnum(bisection.size());
    int o_cnt = 0;
    int z_cnt = 0;
    for(i=1;i<bisection.size();i++)
    {
        if(bisection[i]==0)
        {
            o_cnt++;
            newnum[i] = o_cnt;
        }
        else
        {
            z_cnt++;
            newnum[i] = z_cnt;
        }
    }

    vector <vector <pair<int,int> > >  adj0(o_cnt+1);
    vector <vector <pair<int,int> > >  adj1(z_cnt+1);
    vector<int> bisect0;
    vector<int> bisect1;
    // vector<int> ret(bisection.size());
    // ret[0] = -1;
    int index_0=1, index_1=1;
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
                p.first = newnum[p.first];
                if(bisection[i]==bisection[adj[i][j].first])
                {
                    if(bisection[i]==1)
                        adj1[newnum[i]].pb(p);
                    else
                        adj0[newnum[i]].pb(p);
                }
            }
        }
        
        bisect0 = KWayCut_rec(adj0,k/2,2*offset);
        bisect1 = KWayCut_rec(adj1,k/2,1+2*offset);

        vector<int> bisection_new(bisection.size(),-1);
        for(i=1;i<bisection.size();i++)
        {
            if(bisection[i]==0)
                bisection_new[i] = bisect0[newnum[i]];
            else
                bisection_new[i] = bisect1[newnum[i]];
        }

        // for(i=1;i<bisection.size();i++){
        //     if(bisection[i]==0){
        //         ret[i] = bisect0[index_0];
        //         index_0++;
        //     }
        //     else if(bisection[i]==1){
        //         ret[i] = bisect1[index_1];
        //         index_1++;
        //     }
        //     else{
        //         cout <<"error see k way"<<endl;
        //     }
        // }

        return bisection_new;

}


vector<int> KWayCut(vector <vector <pair<int,int> > >  &adj, int k){
    return KWayCut_rec(adj,k,0);
}


int main(int argc, char* argv[])
{   
    // cout<<"at beginning"<<endl;
    string inf = argv[1];
    string of = argv[2];
    int k = atoi(argv[3]);  
    // cout<<"read arguments"<<endl;
    ifstream infile(inf.c_str());
    // cout<<"opened file"<<endl;
    int nodes, edges;
    infile>>nodes>>edges;
    string line;
    getline(infile,line);
    vector <vector <pair<int,int> > >  adj;
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
    // cout<<"before bisect"<<endl;
    // small_bisect(adj);
    // cout<<"before matching"<<endl;
    // match(adj);
    // return 0;
    //initialising result

    vector<int> res;
    res = KWayCut(adj,k);
    // res.resize(nodes+1);
    // for(i=1;i<=nodes;i++)
    // {
    // 	res[i] = i%k;
    // }
    // writing the result

    // coarsen(adj);
    ofstream outfile(of.c_str());

    for(i=1;i<nodes;i++)
    {
    	outfile<<res[i]<<' ';
    }
    outfile<<res[i]<<"\r\n";
    int cost = 0;
    for(i=0;i<adj.size();i++)
    {
        for(int j = 0;j<adj[i].size();j++)
        {
            if(res[i]!=res[adj[i][j].first])
                cost++;
        }
    }
    cost = cost/2;
    cout<<cost<<endl;
    outfile.close();
    cout<<"We did our stuff"<<endl;
    vector<int> score(k,0);
    for(i=1;i<adj.size();i++)
        score[res[i]]++;
    int mn = score[0];
    int mx = score[0];
    for(i=0;i<k;i++)
    {
        mn = min(mn,score[i]);
        mx = max(mx,score[i]);
    }
    float corr = (float)(mx-mn)/float(adj.size()-1);
    cout<<corr<<endl;
    return 0;
}

