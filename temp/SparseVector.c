#include "SparseVector.h"
#include <stdio.h>
#include <stdlib.h>

void init(SparseVector* sv, int nnz, int dim, int cols[], double vals[]){
  sv->nnz = nnz;
  sv->dimension = dim;
  sv->cols = (int*) malloc(sizeof(int) * nnz);
  sv->values = (double*) malloc(sizeof(double) * nnz);

  for (int i=0; i < nnz; i++){
    //assert(cols[i]>=0);// && cols[i]<dim);
    sv->cols[i] = cols[i];
    sv->values[i] = vals[i];
  }
}

void random_init(SparseVector* sv, int dim){
  sv->dimension = dim;
  
  sv->nnz = rand() % dim;
  sv->cols = (int*) malloc(sizeof(int) * sv->nnz);
  sv->values = (double*) malloc(sizeof(double) * sv->nnz);

  for (int i=0; i < sv->nnz; i++){
    double newly_gen = rand() % dim;
    int flag = 1;
    int j= 0;

    while (flag == 1){
      for (j=0; j < i; j++)
	if (newly_gen == sv->cols[j])
	  break;
      if (j == i)
	flag = 0; 
      else
	newly_gen = rand() % dim;
    }
    sv->cols[i] = newly_gen;
    sv->values[i] = ((double)rand()/(double)(RAND_MAX)); ;
  }
}

void print(SparseVector* sv){
  for (int i=0; i< sv->nnz; i++)
    printf("(%d, %.2f)", sv->cols[i], sv->values[i]  );
  printf("\n");
}
