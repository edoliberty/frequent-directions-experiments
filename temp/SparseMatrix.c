#include "SparseMatrix.h"
#include <stdio.h>
#include <stdlib.h>

void init (SparseMatrix* sp, int nnz, int dim){
  sp->current_nnz = 0;
  sp->columnDim = dim;
  sp->nextRow = 0;
  sp->pointer = 0;

  sp->rows = (int*) malloc(sizeof(int) * nnz);
  sp->cols = (int*) malloc(sizeof(int) * nnz);
  sp->values = (double*) malloc(sizeof(double) * nnz);

}

void append (SparseMatrix* sp, SparseVector* sv){
  for (int i=0; i < sv->nnz ; i++){
    sp->rows[sp->pointer] = sp->nextRow;
    sp->cols[sp->pointer] = sv->cols[i];
    sp->values[sp->pointer] = sv->values[i];
    sp->pointer ++;
  }
  sp->nextRow ++;
  sp->current_nnz += sv->nnz;
}



double** covarianceMult (SparseMatrix *sp, int dimension, int ell, double** denseMat){
  int headPtr = 0;
  int ptr = 0;
  int rowIndex = sp->rows[headPtr];


  
  double* temp = (double*) malloc(sizeof(double) * ell);
  double** product = (double**) malloc(sizeof(double*) * dimension);
  for (int i=0; i < dimension; i++){
    product[i] = (double*) malloc(sizeof(double) * ell);
  }
  
  //double temp[ell];
  //double product[dimension][ell];

  for(int i=0; i < dimension; i++)
    for(int j=0; j<ell ; j++)
      product[i][j] = 0;

  // computation of A^TA * DenseMat
  while (ptr != sp->pointer){
    headPtr = ptr;
    rowIndex = sp->rows[headPtr];
    for (int j=0; j<ell ; j++)
      temp[j] = 0;

    while (ptr != sp->pointer && sp->rows[ptr] == rowIndex){
      for (int t=0; t<ell ; t++)
	temp[t] += denseMat [sp->cols[ptr]] [t] * sp->values[ptr];
      ptr ++;
    }
    for (int j=headPtr; j<ptr ; j++){
      for (int t=0; t<ell ; t++)
	product [sp->cols[j]][t] += temp[t] * sp->values[j];
    }
  }
  
  return product;
}

