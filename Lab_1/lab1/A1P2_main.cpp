#include "iostream"
#include "fstream"
#include "vector"
#include "ctime"
#include "omp.h"
#include "algorithm"

using namespace std;

extern vector< pair<int, int> > calcConvexHull(vector< vector<int> > image, int num_threads);

bool sortBy(const pair<int, int> &a, const pair<int, int> &b){
    if (a.first == a.second){
        return (a.second < b.second);
    }
    else{
        return (a.first < b.first);
    }
}

int main(int argc, char* argv[]) {
    int num_threads = atoi(argv[2]);
    int m, n;

    ifstream in;
    in.open(argv[1]);

    in >> m >> n;
    vector< vector<int> > image;
    image.resize(m);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int a;
            in >> a;
            image[i].push_back(a);
        }
    }

    in.close();

    vector< pair<int, int> > convexHull;
    double start_time = omp_get_wtime();
    convexHull = calcConvexHull(image, num_threads);
    double time_taken = omp_get_wtime() - start_time;

    // Printing stats and results
    cout<< time_taken << endl;

    sort(convexHull.begin(), convexHull.end(), sortBy);

    ofstream out;
    out.open(argv[3]);
    out<< convexHull.size() << endl;

    for (int i = 0; i < convexHull.size(); i++){
        out << convexHull[i].first << " " << convexHull[i].second << endl;
    }
    out.close();

    return 0;
}