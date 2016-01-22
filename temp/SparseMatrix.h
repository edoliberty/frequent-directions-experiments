#include <stdio.h>
#include <stdlib.h>
#include "SparseVector.h"

//typedef struct SparseMatrix SparseMatrix;
typedef struct {
  int* rows;
  int* cols;
  double* values;
  int dimension;
  int pointer;
  int nextRow;
  int current_nnz;
  int arrays_length;
  
} SparseMatrix;


void init_matrix (SparseMatrix* sp, int dim);
void append (SparseMatrix *sp, SparseVector *sv);
void covarianceMult (SparseMatrix *sp, int dimension, int ell, double** G);

double* blockPowerMethod(SparseMatrix *sp, int ell, double epsilon);
double* matrixMult (SparseMatrix *sp, int ell, double* G);

double* transposeRightMult (SparseMatrix *sp, int ell, double* G);
double* rightMult (SparseMatrix *sp, int ell, double* G);
double* leftMult (SparseMatrix *sp, int ell, double* G);
double* sparseShrink(SparseMatrix *sp, int ell);
