#ifndef SPARSESKETCHER_H
#define SPARSESKETCHER_H

#include "sparseMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>

typedef struct {
  char* class_name;
  int dimension;
  int ell;
  int m;
  double* sketch;
  SparseMatrix buffer;
  int nnz_threshold;

} SparseSketcher;


void init_sparseSketcher(SparseSketcher* self, int ell, int dim );
void append_to_sparseSketcher(SparseSketcher* self, SparseVector* sv);
void sparseShrink(SparseSketcher* self);
void denseShrink(SparseSketcher* self);
void rotate_sparseSketcher(SparseSketcher* self);
void get_sparseSketch(SparseSketcher* self);

#endif
