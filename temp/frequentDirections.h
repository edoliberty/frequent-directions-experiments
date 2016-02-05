#ifndef FD_H
#define FD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "sparseVector.h"
#include <lapacke.h>

typedef struct {
  char* class_name;
  int dimension;
  int ell;
  int m;
  int nextRow;
  double* sketch;

} FrequentDirections;


void init_fd(FrequentDirections* self, int ell, int dim );
void append_to_fd(FrequentDirections* self, SparseVector* sv);
void get_fdSketch(FrequentDirections* self);
void rotate_fd(FrequentDirections* self);

#endif
