#include "randomProjection.h"

void init_randomProj(RandomProjection* self, int ell, int dim ){
  self->class_name = "RandomProjection";
  self->dimension = dim;
  self->ell = ell;
  self->sketch = (double*) malloc(dim * ell * sizeof(double));
  memset(self->sketch, 0, sizeof(double) * ell * dim);
  srand(time(NULL));
}


void append_to_randomProj(RandomProjection* self, SparseVector* sv){
  int sign, i, j, index;

  for(i=0; i < self->ell; i++){
    sign = (-2) * (rand() % 2) + 1;
    for(j=0; j < sv->nnz; j++){
      index = i * self->dimension + sv->cols[j];
      self->sketch[index] += (sign/sqrt(self->ell)) * (sv->values[j]);
    }
  }

}

