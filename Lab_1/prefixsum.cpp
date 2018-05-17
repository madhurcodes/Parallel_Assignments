#include "iostream"
#include "vector"
#include "ctime"
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include "omp.h"

using namespace std;

 vector<int> calcPrefixSum(vector<int> input, int num_threads)
{
   int size = input.size();
   int nthr,*tt,*pp;
   vector<int> x(size,0);

  #pragma omp parallel num_threads(num_threads)
  {
     int i;
    #pragma omp single
    {
      nthr = omp_get_num_threads();
    //   printf("Nthr - %d\n",nthr);
      tt = (int *) malloc(sizeof( int)*nthr);    
      pp = (int *) malloc(sizeof( int)*nthr);    
    }
     int tid = omp_get_thread_num();
     int sum = 0;

    #pragma omp for schedule(static) 
    for(i=0; i<size; i++) {
        // printf("Tid - %lld, i = %lld\n",tid,i);
      sum += input[i];
      x[i] = sum;
    }
    pp[tid] = sum;
    #pragma omp barrier
    // printf("tid- %lld,sum - %lld\n",tid,sum);

    for(i=1; i<nthr; i*=2) {
        if(tid>=i){
            tt[tid] = pp[tid] + pp[tid-i];
        }   
        #pragma omp barrier
        #pragma omp single
        memcpy(pp + 1, tt + 1, sizeof( int) * (nthr - 1));
    }

    // #pragma omp single
    // {
    //      int k;
    //     for(k=0;k<nthr;k++){
    //         printf("tid - %lld , tt - %lld\n",k,pp[k]);
    //     }
    // }
    #pragma omp for schedule(static)
    for(i=0; i<size; i++) {
    if(tid>=1){
      x[i] += pp[tid-1];
    }
    //   printf("i - %lld,sum - %lld \n",i,x[i*8]);
    }
  }
  free(tt);
  free(pp);
  return x;
}
