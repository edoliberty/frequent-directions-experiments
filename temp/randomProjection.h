#ifndef RANDPROJ_H
#define RANDPROJ_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "sparseVector.h"

typedef struct {
  char* class_name;
  int dimension;
  int ell;
  double* sketch;

} RandomProjection;


void init_randomProj(RandomProjection* self, int ell, int dim );
void append_to_randomProj(RandomProjection* self, SparseVector* sv);

#endif
