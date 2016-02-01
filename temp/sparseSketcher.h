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
  double* sketch;
  SparseMatrix buffer;
  int nnz_threshold;

} SparseSketcher;


void init_sketcher(SparseSketcher* ss, int ell, int dim );
void append_row(SparseSketcher* ss, SparseVector* sv);
void sparseShrink(SparseSketcher* ss);
void denseShrink(SparseSketcher* ss);
void rotate(SparseSketcher* ss);
void get_sketch(SparseSketcher* ss);
