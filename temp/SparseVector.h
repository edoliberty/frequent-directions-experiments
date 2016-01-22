#include <stdio.h>
#include <stdlib.h>

typedef struct{
  int dimension;
  int* cols;
  double* values;
  int nnz;

} SparseVector;

void init(SparseVector* sv, int dim, int cols[], double vals[]);

void random_init(SparseVector* sv, int dim);

void print(SparseVector* sv);
