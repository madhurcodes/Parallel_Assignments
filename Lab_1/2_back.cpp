#include "iostream"
#include "vector"
#include "omp.h"
#include<set>

using namespace std;

#define point pair<int, int>

set<point> hull;

int findSide(point p1, point p2, point p)
{
    int val = (p.second - p1.second) * (p2.first - p1.first) -
              (p2.second - p1.second) * (p.first - p1.first);
 
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
    return abs ((p.second - p1.second) * (p2.first - p1.first) -
               (p2.second - p1.second) * (p.first - p1.first));
}
 

void quickHull(vector<point> &a, int n, point p1, point p2, int side)
{
    int ind = -1;
    int max_distance = 0;
 
    for (int i=0; i<n; i++)
    {
        int temp = linedistance(p1, p2, a[i]);
        if (findSide(p1, p2, a[i]) == side && temp > max_distance)
        {
            ind = i;
            max_distance = temp;
        }
    }

    if (ind == -1)
    {
        hull.insert(p1);
        hull.insert(p2);
        return;
    }
 
    quickHull(a, n, a[ind], p1, -findSide(a[ind], p1, p2));
    quickHull(a, n, a[ind], p2, -findSide(a[ind], p2, p1));
}
 
void printHull(vector<point> &a, int n)
{
    if (n < 3)
    {
        cout << "Convex hull not possible\n";
        return;
    }
 
    int min_x = 0, max_x = 0;
    for (int i=0; i<n; i++)
    {
        if (a[i].first < a[min_x].first)
            min_x = i;
        if (a[i].first > a[max_x].first)
            max_x = i;
    }
 
    quickHull(a, n, a[min_x], a[max_x], 1);
 
    quickHull(a, n, a[min_x], a[max_x], -1);
 
    cout << "The points in Convex Hull are:\n";
    while (!hull.empty())
    {
        cout << "(" <<( *hull.begin()).first << ", "
             << (*hull.begin()).second << ") "<<endl;
        hull.erase(hull.begin());
    }
}
// vector< point > calcConvexHull(vector< vector<int> > image, int num_threads) {
vector< point > calcConvexHull(int image, int num_threads) {
    
    vector< point> inp;
    
    vector< point> test;
    test.push_back(make_pair(0,0));
    test.push_back(make_pair(0,3));
    test.push_back(make_pair(1,1));
    test.push_back(make_pair(2,2));
    test.push_back(make_pair(1,2));
    test.push_back(make_pair(4,4));
    test.push_back(make_pair(3,3));
    test.push_back(make_pair(3,1));
    printHull(test,8);
    // #omp pragma parallel 
    // {


    // }
    vector<point> ret( hull.begin(), hull.end() );
    return ret;
}


int main(){
    calcConvexHull(1,1);
    return 0;
}
