#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "omp.h"
#define N 1499336 //33554432 // 2 ^ 25

void comp(long long int *input_array, long long int *out, unsigned long size)
{
  long long int nthr,*tt,*pp, *x = out;
  #pragma omp parallel num_threads(3)
  {
    long long int i;
    #pragma omp single
    {
      nthr = omp_get_num_threads();
      printf("Nthr - %lld\n",nthr);
      tt = malloc(sizeof(long long int)*nthr);    
      pp = malloc(sizeof(long long int)*nthr);    
    }
    long long int tid = omp_get_thread_num();
    long long int sum = 0;

    #pragma omp for schedule(static) 
    for(i=0; i<size; i++) {
        // printf("Tid - %lld, i = %lld\n",tid,i);
      sum += input_array[i];
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
        memcpy(pp + 1, tt + 1, sizeof(long long int) * (nthr - 1));
    }

    #pragma omp single
    {
        long long int k;
        for(k=0;k<nthr;k++){
            // printf("tid - %lld , tt - %lld\n",k,pp[k]);
        }
    }    
    #pragma omp for schedule(static)
    for(i=0; i<size; i++) {
    if(tid>=1){
      x[i] += pp[tid-1];
    }
    //   printf("i - %lld,sum - %lld \n",i,x[i]);
    }
  }
  free(tt);
  free(pp);
  
}

int main(void ) {

  long long int *input_array, *pprefixsum ;

  input_array = (long long int*) malloc(sizeof(long long int) * N);
  pprefixsum = (long long int*) malloc(sizeof(long long int) * N);

  for(long long int i=0; i<N; i++) input_array[i] = i+1;
//   for(long long int i=0; i<N; i++) printf("%d ", input_array[i]); printf("\n");
    double start = omp_get_wtime();

  comp(input_array, pprefixsum, N);
  double end = omp_get_wtime();
  printf("Time Taken - %f s.\n",end-start);    
//   for(long long int i=0; i<N; i++) printf("%d ", pprefixsum[i]); printf("\n");
//   for(long long int i=0; i<N; i++) printf("%d ", (i+1)*(i+2)/2); printf("\n");
  for(long long int i=0; i<N; i++) {
      if(pprefixsum[i]!=(i+1)*(i+2)/2 ){
          printf("EEEEEEEEError");
          break;
      }
  }
  free(input_array);
  free(pprefixsum);
  return 0;
}