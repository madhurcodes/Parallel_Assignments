#include "iostream"
#include "vector"
#include "ctime"
#include "omp.h"

using namespace std;

extern vector<int> calcPrefixSum(vector<int> input, int num_threads);

int main() {
    int num_threads;
    cin >> num_threads;
    int n;
    cin >> n;
    vector<int> input;
    input.resize(n);
    vector<int> des_out;
    des_out.resize(n);
    for (int i = 0; i < n; i++) {
        input[i]  = -9 + (rand() % static_cast<int>(18 + 1));
        // cout<<input[i]<<endl;
    }

    des_out[0] = input[0];
    for (int i = 1; i < n; i++) {
        des_out[i]  = des_out[i-1]+input[i];
    }

    vector<int> prefixSum;
    double start_time = omp_get_wtime();
    prefixSum = calcPrefixSum(input, num_threads);
    double time_taken = omp_get_wtime() - start_time;

    // Printing stats and results
    cout<< time_taken << endl;
    cout<< prefixSum.size() << endl;

    for (unsigned i = 0; i < prefixSum.size(); i++){
        if(prefixSum[i]!=des_out[i]){
            cout<<"ERROR"<<endl;
            return -1;
        }
        // cout << prefixSum[i] << " " ;
    }
    cout << endl;

    return 0;
}