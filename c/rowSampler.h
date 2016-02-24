#ifndef ROWSAMPLER_H
#define ROWSAMPLER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "sparseVector.h"
#include "singleItemSampler.h"

typedef struct {
  char* class_name;
  int dimension;
  int ell;
  double* sketch;
  SingleItemSampler* samplers;

} RowSampler;


void init_rowSampler(RowSampler* self, int ell, int dim );
void append_to_rowSampler(RowSampler* self, SparseVector* sv);
void get_rowSamplerSketch(RowSampler* self);

#endif
