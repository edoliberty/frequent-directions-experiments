#include <stdio.h>
#include <stdlib.h>
#include "SparseVector.h"

//typedef struct SparseMatrix SparseMatrix;
typedef struct {
  int* rows;
  int* cols;
  double* values;
  int columnDim;
  int pointer;
  int nextRow;
  int current_nnz;
  
} SparseMatrix;


void init (SparseMatrix* sp, int nnz, int dim);
void append (SparseMatrix *sp, SparseVector *sv);
//int* getShape (SparseMatrix *self);
double** covarianceMult (SparseMatrix *sp, int dimension, int ell, double** denseMat);
