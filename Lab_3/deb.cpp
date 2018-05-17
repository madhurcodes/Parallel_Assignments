#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <cstring>

using namespace std;

int main(int argc , char * argv[]){

    MPI_Init(&argc, &argv); 
    int num_columns = atoi(argv[1]);
    char* base_filename = argv[2];

    // char fname[80];
    // char num[10];
    // sprintf(num,"%d",num_columns);
    // strcpy(fname,base_filename);
    // strcat(fname,num);
    // cout<<"fname si "<<fname<<endl;
    int prank; MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    int psize; MPI_Comm_size(MPI_COMM_WORLD, &psize);



    vector<vector<pair<float,char *> > > mycolumns;
    int n_excess_columns = num_columns % psize;
    int n_my_columns;
    if (prank < n_excess_columns){
        n_my_columns = (num_columns/psize) + 1;
    }
    else{
        n_my_columns = (num_columns/psize);
    }
    int i = 0;

    vector <int> col_to_proc_rank(num_columns);
    int j = 0;
    int k = 0;
    for (i=0;i<num_columns;i++){

        col_to_proc_rank[i] = j;
        if (j<n_excess_columns){
            if(k>=(num_columns/psize)){
                j+=1;
                k = 0;
            }
            else{
                k++;
            }
        }
        else{
            if(k>=(num_columns/psize)-1){
                j+=1;
                k = 0;
            }
            else{
                k++;
            }

        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(prank==0){
    for(i=0;i<num_columns;i++){
        cout<<"Col - "<< i<< "Rank " << col_to_proc_rank[i]<<endl;
    }
    }
    MPI_Finalize();
    
    return 0;
}