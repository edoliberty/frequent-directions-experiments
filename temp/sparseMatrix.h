#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

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
  double squaredFrob;

} SparseMatrix;


void init_sparseMatrix (SparseMatrix* self, int dim, int init_len);
void extend_sparseMatrix(SparseMatrix* self);
void append_to_sparseMatrix (SparseMatrix *self, SparseVector *sv);
void print_sparseMatrix(SparseMatrix* self);
void qrDecomp(double* G, lapack_int d, lapack_int ell);
void covarianceMult (SparseMatrix *self, int dimension, int ell, double** G, double* temp, double** product);
void leftMult (SparseMatrix *self, int ell, double* G, double* product);
double* rightMult (SparseMatrix *self, int ell, double* G);
void transposeRightMult (SparseMatrix *self, int ell, double* G, double* product);
void blockPowerMethod(SparseMatrix *self, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat);
double* getCovariance_sparseMatrix(SparseMatrix *self);
void densify_sparseMatrix(SparseMatrix* self, double* output);

#endif
