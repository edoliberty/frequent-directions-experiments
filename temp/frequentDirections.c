#include "frequentDirections.h"

void init_fd(FrequentDirections* self, int ell, int dim ){
  self->class_name = "FrequentDirections";
  self->dimension = dim;
  self->ell = ell;
  self->m = 2*ell;
  self->sketch = (double*) malloc(dim * (self->m) * sizeof(double));
  memset(self->sketch, 0, sizeof(double) * self->m * dim);
  self->nextRow = 0;
}


void append_to_fd(FrequentDirections* self, SparseVector* sv){
  //if (sv->nnz == 0)
  //  return;

  if (self->nextRow == self->m)
    rotate_fd(self);
  

  int i, index, j = 0;
  for(i = 0; i < sv->dimension; i++){ 
    index = (self->nextRow) * (self->dimension) + i;
    if(j < sv->dimension && sv->cols[j] == i ){
      self->sketch[index] = sv->values[j];
      j++;
    }else{
      self->sketch[index] = 0;
    }
  }
  
  self->nextRow ++;
}


void rotate_fd(FrequentDirections* self){
  //self->nextRow = 0;
  
  double S[self-> m], U[(self-> m) * (self-> m)], Vt[(self->m) * self->dimension];
  int i, j, info;

  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', self->m, self->dimension, self->sketch, self->dimension, S, U, self->m, Vt, self->dimension);

  // shrink S
  for(i=0; i < self->ell; i++){
    S[i] = sqrt( pow(S[i],2) - pow(S[self->ell - 1],2) );
  }

  // compute S*Vt
  for(i=0; i < self->ell; i++)
    for(j=0; j < self->dimension; j++)
      self->sketch[i * self->dimension + j] = Vt[i * self->dimension + j] * S[i] ;

  self->nextRow = self->ell;  
}

void get_fdSketch(FrequentDirections* self){
  memset(&(self->sketch[self->nextRow * self->dimension]), 0, (self->m - self->nextRow) * self->dimension * sizeof(double));
  rotate_fd(self);
}
