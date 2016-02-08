#include "sparseVector.h"


void init_sparseVector(SparseVector* self, int dim, int cols[], double vals[]){
  self-> nnz = sizeof(cols) / sizeof(cols[0]);
  self-> dimension = dim;
  self-> cols = (int*) malloc(sizeof(int) * self-> nnz);
  self-> values = (double*) malloc(sizeof(double) * self-> nnz);
  self-> squaredNorm = 0;
  int i;

  for (i=0; i < self-> nnz; i++){
    self-> cols[i] = cols[i];
    self-> values[i] = vals[i];
    self-> squaredNorm += pow(vals[i] , 2);
  }
}

void random_init_sparseVector(SparseVector* self, int dim, int nnz){
  self-> dimension = dim;  
  self-> nnz = nnz;
  self-> cols = (int*) malloc(sizeof(int) * self-> nnz);
  self-> values = (double*) malloc(sizeof(double) * self-> nnz);
  self-> squaredNorm = 0;
  int i;

  for (i=0; i < self-> nnz; i++){
    double newly_gen = rand() % dim;
    int flag = 1;
    int j= 0;

    while (flag == 1){
      for (j=0; j < i; j++)
	if (newly_gen == self-> cols[j])
	  break;
      if (j == i)
	flag = 0; 
      else
	newly_gen = rand() % dim;
    }
    self-> cols[i] = newly_gen;
    self-> values[i] = (int)ceil( ((double)rand()/(double)(RAND_MAX)) * 10);
    self-> squaredNorm += pow(self-> values[i] , 2);
  }
}

void print_sparseVector(SparseVector* self){
  int i;
  printf("nnz=%d",self-> nnz);
  for (i=0; i< self-> nnz; i++)
    printf("(%d, %.2f)", self-> cols[i], self-> values[i]  );
  printf("\n");
}

/*
void printDense(SparseVector* self){
  //sort first
  int i, j=0, k;
  printf("[");
  for (i=0; i < self-> nnz; i++){
    for (k=j; k < self-> cols[i]; k++)
      printf("0 ");
    printf("%f ", self-> values[i]);
    j = self-> cols[i];
  }
  printf("]\n");
}
*/
