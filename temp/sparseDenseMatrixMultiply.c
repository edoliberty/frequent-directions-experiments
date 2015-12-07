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
  
  int nonzeros = nonzerosPerRow*m;

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


  float X[d][k];
  float X2[d][k];
  
  float temp[d];

  clock_t start, end;
  double cpu_time_used;
  start = clock();
  int rounds = 50;
  for (int round = 0; round < rounds; round++){
    int t1=0,t2;
    for (int row=0; row<m; row++){ // for every row in A                                                                                                                                                                                                                                                                                                                     
      for(int j=0;j<d;j++) temp[j]=0; // zeroing out temp                                                                                                                                                                                                                                                                                                                    
      while (rows[t1]<row && t1<nonzeros) t1++;
      if (rows[t1]>row) continue;
      // at this point rows[t1]==row                                                                                                                                                                                                                                                                                                                                         
      t2=t1;
      while (rows[t2]==row && t2<nonzeros) {
	t2++;
      }
      // at this point t2-t1 is the number of non zeros in the row
      
      // computing the temp vector
      for (int t=t1; t<t2; t++){
	int col = cols[t];
	float val = vals[t];
	float* x = X[col];
	for (int j=0;j<k;j++){
	  temp[j]+=val*x[j];
	}
      }
      
      for (int t=t1; t<t2; t++){
	int col = cols[t];
	float val = vals[t];
	float* x2 = X2[col];
	for (int j=0;j<k;j++){
	  x2[j]+=val*temp[j];
	}
      }
    }
  }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("Done in %f seconds\n",cpu_time_used);
}

