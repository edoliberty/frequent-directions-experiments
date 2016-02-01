#include <stdio.h>
//#include "sparseMatrix.h"
#include <stdlib.h>
#include <time.h>
#include "sparseSketcher.h"

int main(){
 
  int n = 10000;
  int dim = 1000;
  int ell = 30;
  int i;
  double start, end, cpu_time_used;
  
  SparseSketcher A;
  init_sketcher(&A, ell, dim);
  
  SparseVector arr[n];
  for (i=0; i < n; i++){
    random_init(&arr[i], dim);
  }
  start = clock();
  for (i=0; i < n; i++)
    append_row(&A, &arr[i]);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Done in %f seconds\n",cpu_time_used);

  get_sketch(&A);
  print_two_dim("sketch = ", A.sketch, ell, dim);

  /*
  // QR test:begin
  double* G = (double*) malloc(sizeof(double) * dim * ell);
  for(i=0; i < dim*ell; i++)
    G[i] = ((float)rand()/(float)(RAND_MAX));
  //print_two_dim("G",G,ell, dim);
  //printline();

 
  init_matrix(&A, dim);
  SparseVector arr[n];
  for (i=0; i < n; i++){
    random_init(&arr[i], dim);
  }

  for (i=0; i < n; i++)
    append(&A, &arr[i]);
 
  
  double* p = (double*) malloc(sizeof(double) * ell * A.nextRow);
  double* temp_vec = (double*) malloc(sizeof(double) * ell);
  double* temp_mat = (double*) malloc(sizeof(double) * ell * A.dimension);


  start = clock(); 
  //double* B = sparseShrink(&A, ell);
  for(i=0; i< rounds; i++)
    ;
    //covarianceMult(&A, dim, ell, &G, temp, &product);
    //blockPowerMethod(&A, ell, 0.25, G, p, temp_vec, temp_mat);
    //temp_mat = sparseShrink(&A, ell);
  end = clock();

  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Done in %f seconds\n",cpu_time_used);
  
  //print_two_dim("cov = ",G, dim, ell);
  //print_two_dim("block power = ", p, A.nextRow, ell);
  */
}
