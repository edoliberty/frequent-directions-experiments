#include "frequentDirections.h"

void init_fd(FrequentDirections* self, int ell, int dim ){
  self->class_name = "FrequentDirections";
  self->dimension = dim;
  self->ell = ell;
  self->m = 2*ell;
  self->sketch = (double*) malloc(dim * (self->m) * sizeof(double));
  self->nextRow = 0;
}


void append_to_fd(FrequentDirections* self, SparseVector* sv){

  if (self->nextRow == self->m)
    rotate_fd(self);
  

  int j = 0;
  int rid = (self->nextRow) * (self->dimension);

  double* vec = densify_sparseVector(sv);
  

  for(int i = 0; i < sv->dimension; i++) 
    self->sketch[rid + i] = vec[i];

  self->nextRow ++;
  free(vec);  
}


void rotate_fd(FrequentDirections* self){
  double* S = (double*) malloc(sizeof(double) * self->m);
  double* U = (double*) malloc(sizeof(double) * self->m * self->m);
  double* Vt = (double*) malloc(sizeof(double) * self->m * self->dimension);


  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', self->m, self->dimension, self->sketch, self->dimension, S, U, self->m, Vt, self->dimension);


  // compute S*Vt
  for(int i=0; i < self->ell; i++){
    S[i] = sqrt( pow(S[i],2) - pow(S[self->ell - 1],2) );
    for(int j=0; j < self->dimension; j++)
      self->sketch[i * self->dimension + j] = Vt[i * self->dimension + j] * S[i] ;
  }

  memset(&self->sketch[self->ell * self->dimension], 0, self->ell * self->dimension * sizeof(double));


  self->nextRow = self->ell;  
  free(S); free(U); free(Vt); 
}

void get_fdSketch(FrequentDirections* self) {}
