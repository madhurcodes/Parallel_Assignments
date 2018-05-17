#include "iostream"
#include "vector"
#include "omp.h"
#include <set>

using namespace std;

#define point pair<int, int>

int findSide(point p1, point p2, point p)
{
    int val = (p.second - p1.second) * (p2.first - p1.first) - (p2.second - p1.second) * (p.first - p1.first);
 
    if (val > 0)
        return 1;
    if (val < 0)
        return -1;
    return 0;
}
 
int distance(point p, point q)
{
    return (p.second - q.second) * (p.second - q.second) + (p.first - q.first) * (p.first - q.first);
}

int linedistance(point p1, point p2, point p)
{
    return abs ((p.second - p1.second)*(p2.first - p1.first)-(p2.second - p1.second) * (p.first - p1.first));
}
 
void quickHull(vector<point> &a,set<point> &hull,  int p1, int p2, int side,int nthr_tot)
{
    int ind = -3;
    int max_distance = -1;
    
    // #pragma omp parallel for  num_threads(nthr) reduction(max:max_distance)
    for (unsigned i=0; i<a.size(); i++)
    {
        int temp = linedistance(a[p1], a[p2], a[i]);
        if (findSide(a[p1], a[p2], a[i]) == side && temp > max_distance)
        {
            max_distance = temp;
            ind = i;
        }
    }

    if (ind == -3)
    {   
        #pragma omp critical
        {
        hull.insert(a[p1]);
        hull.insert(a[p2]);
        }
    }
    else{
        #pragma omp parallel num_threads(nthr_tot)
        {
        #pragma omp sections
        {
        #pragma omp section
        {
            // cout<< "tid "<<omp_get_thread_num()<< " total "<<omp_get_num_threads()<<endl;
        quickHull(a,hull, ind, p1, -findSide(a[ind], a[p1], a[p2]),nthr_tot);
        }
        #pragma omp section
        {
            // cout<< "tid "<<omp_get_thread_num()<< " total "<<omp_get_num_threads()<<endl;            
        quickHull(a,hull, ind, p2, -findSide(a[ind], a[p2], a[p1]),nthr_tot);
        }
        }
        }
    }
}
 

vector< point > calcConvexHull(vector< vector<int> > image, int nthr) {
// vector< point > calcConvexHull(int image, int nthr) {
    unsigned i,j;
    vector< point> inp;
    for(i=0;i<image.size();i++){
        for(j=0;j<image[0].size();j++){
            if(image[i][j]==1){
                inp.push_back(make_pair(i,j));
            }
        }
    }

    set<point> hull;
    // vector< point> test;

    // test.push_back(make_pair(1,1));
    // test.push_back(make_pair(2,2));
    // test.push_back(make_pair(1,2));
    // test.push_back(make_pair(0,0));
    // test.push_back(make_pair(0,3));
    // test.push_back(make_pair(4,4));
    // test.push_back(make_pair(3,3));
    // test.push_back(make_pair(3,1));
    // test.push_back(make_pair(-4,3));
    // test.push_back(make_pair(7,3));
    // test.push_back(make_pair(5,-3));
    // test.push_back(make_pair(11,13));
    // test.push_back(make_pair(4,-9));
    // test.push_back(make_pair(4,0));
    // test.push_back(make_pair(-3,22));
    // test.push_back(make_pair(-2,8));

    // vector<point> a = inp;

    if (inp.size() <= 2)
    {
        // cout << "Error";
        return inp;
    }
 
    int minx = 0, maxx = 0;
    for (unsigned i=0; i<inp.size(); i++)
    {
        if (inp[i].first < inp[minx].first)
            minx = i;
        if (inp[i].first > inp[maxx].first)
            maxx = i;
    }
 
    omp_set_nested(1); 
    quickHull(inp, hull,  minx, maxx, 1,nthr);
    quickHull(inp, hull, minx, maxx, -1,nthr);

    // cout << "The points are:\n";

    // while (!hull.empty())
    // {
    //     cout << "(" <<( *hull.begin()).first << ", "
    //          << (*hull.begin()).second << ") "<<endl;
    //     hull.erase(hull.begin());
    // }


    vector<point> ret( hull.begin(), hull.end() );
    return ret;
}


// int main(){
//     calcConvexHull(1,3);
//     return 0;
// }
