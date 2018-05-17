#include <stdio.h>
#include<omp.h>
int main ()  
{
 #pragma omp parallel
 {
 int i = omp_get_thread_num();
 
 printf("Hello %d ",i);
 printf("world! %d \n",i);
 } 
 }
