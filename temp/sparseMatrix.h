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
  SparseVector* vectors;
  int nextRow;
  int dimension;
  int current_nnz;
  double squaredFrob;

} SparseMatrix;


void init_sparseMatrix (SparseMatrix* self, int dim, int len);
void append_to_sparseMatrix (SparseMatrix *self, SparseVector *sv);
void print_sparseMatrix(SparseMatrix* self);
void covMultiply_sparseMatrix (SparseMatrix *self, int dimension, int ell, double** G, double* temp, double** product);
void leftMult (SparseMatrix *self, int ell, double* G, double* product);
void transposeRightMult (SparseMatrix *self, int ell, double* G, double* product);
void blockPowerMethod(SparseMatrix *self, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat);
double* getCovariance_sparseMatrix(SparseMatrix *self);
void densify_sparseMatrix(SparseMatrix* self, double* output);

#endif
