#ifndef RANDSUM_H
#define RANDSUM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "sparseVector.h"

typedef struct {
  char* class_name;
  double* sketch;
  int dimension;
  int ell;

} RandomSum;


void init_randomSum(RandomSum* self, int ell, int dim );
void append_to_randomSum(RandomSum* self, SparseVector* sv);

#endif
