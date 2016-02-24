#include "randomSum.h"

void init_randomSum(RandomSum* self, int ell, int dim ){
  self->class_name = "RandomSum";
  self->dimension = dim;
  self->ell = ell;
  self->sketch = (double*) malloc(dim * ell * sizeof(double));
  memset(self->sketch, 0, sizeof(double) * ell * dim);
  srand(time(NULL));
}


void append_to_randomSum(RandomSum* self, SparseVector* sv){
  int rid = rand() % (self->ell);
  int sign = (-2) * (rand() % 2) + 1;
  int index;
  for(int i=0; i<sv->nnz; i++){
    index = rid * self->dimension + sv->cols[i];
    self->sketch[index] += sign * (sv->values[i]);
  }
}

