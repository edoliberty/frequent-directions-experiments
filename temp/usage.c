#include <stdio.h>
//#include "SparseVector.h"
#include "SparseMatrix.h"
#include <stdlib.h>

int main(){

  SparseMatrix sp;
  int dimension = 30;
  int nnz = 15;
  init(&sp , nnz, dimension);
  
  SparseVector arr[10];
  for (int i=0; i < 10; i++){
    random_init(&arr[i], dimension);
  }
  
  for (int i=0; i < 10; i++){
//print(&arr[i]);
    append(&sp, &arr[i]);
  }

  int ell = 10;
  double** denseMat = (double**) malloc(sizeof(double*) * dimension);

  for (int i=0; i < dimension; i++){
    denseMat[i] = (double*) malloc(sizeof(double) * ell);
  }

  for(int i=0; i < dimension; i++)
    for(int j=0; j<ell ; j++)
      denseMat[i][j] = ((float)rand()/(float)(RAND_MAX));

  double** t = covarianceMult(&sp, dimension, ell, denseMat);
 
  for(int i=0; i < dimension; i++){
    for(int j=0; j<ell ; j++)
      printf("%.3f ", t[i][j]);
    printf("\n");
  }

  /*
  SparseVector svarr[10];
  int cols[3] = {1,4,5};
  double vals[3] = {0.1,0.2,0.3};
  init(&sv, 3, dimension, cols , vals);

  SparseVector sv2;
  int cols2[4] = {0,5,9,8};
  double vals2[4] = {0.4,0.21,1.3,4};
  init(&sv2, 4, dimension, cols2 , vals2);

  append(&sp, &sv);
  append(&sp, &sv2);
  printf("%d %d %d ",sp.columnDim, sp.nextRow, sp.pointer);
  */
}
