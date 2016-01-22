#include <stdio.h>
#include "SparseMatrix.h"
#include <stdlib.h>
#include <time.h>

int main(){

  SparseMatrix A;
  int dim = 1000;
  int n = 10000;
  int ell = 300;
  int i;
  double start, end, cpu_time_used;
  
  // QR test:begin
  double* G = (double*) malloc(sizeof(double) * dim * ell);
  for(i=0; i < dim*ell; i++)
    G[i] = ((float)rand()/(float)(RAND_MAX));
  //print_two_dim("G before QR=",G, dim, ell);
  //printline();
  //qrDecomp(G, dim, ell);
  //print_two_dim("G after QR=",G, dim, ell);
  //printline();
  // QR test: end
  

  start = clock();
  init_matrix(&A, dim);
  SparseVector arr[n];
  for (i=0; i < n; i++){
    random_init(&arr[i], dim);
  }

  for (i=0; i < n; i++)
    append(&A, &arr[i]);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("for append in %f seconds\n",cpu_time_used);

  start = clock(); 
  //double* B = sparseShrink(&A, ell);
  for(i=0; i< 100; i++)
    covarianceMult(&A, dim, ell, &G);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Done in %f seconds\n",cpu_time_used);
 
}
