#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int cmpfunc (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}


int main()
{
  int m = 10000;
  int d = 1000;
  int k = 300;
  int nonzerosPerRow = 10;  
  int nonzeros = nonzerosPerRow * m;

  //This is the sparse matrix
  int rows[nonzeros];
  int cols[nonzeros];
  float vals[nonzeros];

  for (int t=0;t<nonzeros;t++){
    rows[t] = rand() % m;
    cols[t] = rand() % d;
    vals[t] = ((float)rand()/(float)(RAND_MAX));
  }
  qsort(rows, nonzeros, sizeof(int), cmpfunc);


  float denseMat[d][k];
  float product[d][k];
  float temp[k]; 
  int headPtr = 0;
  int rowIndex = rows[headPtr];
  int ptr = 0;

  clock_t start,end;
  double cpu_time_used;


  for(int i=0; i < d; i++)
    for(int j=0; j<k ; j++)
      denseMat[i][j] = ((float)rand()/(float)(RAND_MAX));

  for(int i=0; i < d; i++)
    for(int j=0; j<k ; j++)
      product[i][j] = 0;
  

  start = clock();
  while (ptr != nonzeros){
    headPtr = ptr;
    rowIndex = rows[headPtr];
    for (int j=0; j<k ; j++)
      temp[j] = 0;

    while (ptr != nonzeros && rows[ptr] == rowIndex){
      for (int t=0; t<k ; t++)
	temp[t] += denseMat [cols[ptr]] [t] * vals[ptr];
      ptr ++;
    }
    for (int j=headPtr; j<ptr ; j++){
      for (int t=0; t<k ; t++)
	product [cols[j]][t] += temp[t] * vals[j];
    }
  }

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("Done in %f seconds\n",cpu_time_used);

}
