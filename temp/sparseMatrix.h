#include <stdio.h>
#include <stdlib.h>
#include "sparseVector.h"

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
void covarianceMult (SparseMatrix *sp, int dimension, int ell, double** G, double* temp, double** product);

void blockPowerMethod(SparseMatrix *sp, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat);
double* matrixMult (SparseMatrix *sp, int ell, double* G);

void transposeRightMult (SparseMatrix *sp, int ell, double* G, double* product);
double* rightMult (SparseMatrix *sp, int ell, double* G);
void leftMult (SparseMatrix *sp, int ell, double* G, double* product);
//double* sparseShrink(SparseMatrix *sp, int ell);
