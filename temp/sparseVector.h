#ifndef SPARSEVEC_H
#define SPARSEVEC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct{
  double squaredNorm;
  double* values;
  int dimension;
  int* cols;
  int nnz;

} SparseVector;

void init_sparseVector(SparseVector* self, int dim, int cols[], double vals[]);
void random_init_sparseVector(SparseVector* self, int dim);
void print_sparseVector(SparseVector* self);
//void printDense(SparseVector* self);

#endif
