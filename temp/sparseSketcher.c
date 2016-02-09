#include "sparseSketcher.h"

void init_sparseSketcher(SparseSketcher* self, int ell, int dim ){
  self->class_name = "sparseSketcher";
  self->dimension = dim;
  self->ell = ell;
  self->m = 2*ell;
  self->sketch = (double*) malloc(sizeof(double) * (self->m) * dim);
  memset(self->sketch, 0, sizeof(double) * (self->m) * dim);
  init_sparseMatrix(&(self->buffer), dim, dim);
  self->nnz_threshold = ell * dim;
}  


void append_to_sparseSketcher(SparseSketcher* self, SparseVector* sv){
  if((self->buffer).current_nnz >= self->nnz_threshold || (self->buffer).nextRow >= self->dimension)
    rotate_sparseSketcher(self);
  append_to_sparseMatrix(&(self->buffer), sv);
}

void rotate_sparseSketcher(SparseSketcher *self){
  sparseShrink(self);
  denseShrink(self);
}

void get_sparseSketch(SparseSketcher *self){
  rotate_sparseSketcher(self);
}

void sparseShrink(SparseSketcher *self){
  if((self->buffer).nextRow > self->ell){
    double* temp_vec = (double*) malloc(sizeof(double) * self->ell);
    double* temp_mat = (double*) malloc(sizeof(double) * self->ell * self->dimension);
    double* G = (double*) malloc(self->ell * self->dimension * sizeof(double));
    double* Z = (double*) malloc(self->ell * (self->buffer).nextRow * sizeof(double));

    for(int i=0; i < self->ell * self->dimension; i++)
      G[i] = ( (float)rand() / (float)(RAND_MAX) );

    blockPowerMethod(&(self->buffer), self->ell, 1, G, Z, temp_vec, temp_mat);
    free(temp_vec); 
    free(G);

    //computing P = ZtA, temp_mat is P
    transposeRightMult(&(self->buffer), self->ell, Z, temp_mat);
    free(Z);

    // svd(ZtA)
    double S[self->ell], U[(self->ell) * (self->ell)], Vt[self->dimension * (self->ell)];
  
    //double* S = (double*) malloc(sizeof(double) * self->ell);
    //double* U = (double*) malloc(sizeof(double) * self->ell * self->ell);
    //double* Vt = (double*) malloc(sizeof(double) * self->dimension * self->ell);

    int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', self->ell, self->dimension, temp_mat, self->dimension, S, U, self->ell, Vt, self->dimension);
    free(temp_mat);


    // shrink S and compute S*Vt
    for(int i=0; i < self->ell; i++){
      S[i] = sqrt( pow(S[i],2) - pow(S[self->ell-1],2) );
      for(int j=0; j < self->dimension; j++)
	self->sketch[(self->ell + i) * self->dimension + j] = Vt[i * self->dimension + j] * S[i] ;
    }
  }else{ // self->buffer has atmost ell rows
    SparseVector temp;
    int itr = (self->buffer).nextRow;

    for(int i=0; i < itr; i++){
      temp = (self->buffer).vectors[i];
      for(int j=0; j < temp.nnz; j++)
	self->sketch[(self->ell + i) * self->dimension + temp.cols[j]] = temp.values[j];
    }
  }

  // reset buffer
  (self->buffer).current_nnz = 0;
  (self->buffer).nextRow = 0;
  (self->buffer).squaredFrob = 0;
}


void denseShrink(SparseSketcher* self){
  double S[2*self->ell], U[(2*self->ell) * (2*self->ell)], Vt[self->dimension * (2*self->ell)];

  //double* S = (double*) malloc(sizeof(double) * 2 * self->ell);
  //double* U = (double*) malloc(sizeof(double) * 4 * self->ell * self->ell);
  //double* Vt = (double*) malloc(sizeof(double) * 2 * self->dimension * self->ell);

  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', 2*self->ell, self->dimension, self->sketch, self->dimension, S, U, 2*self->ell, Vt, self->dimension);

  for(int i=0; i < self->ell; i++){
    S[i] = sqrt( pow(S[i],2) - pow(S[self->ell-1],2) );
    for(int j=0; j < self->dimension; j++)
      self->sketch[i * self->dimension + j] = Vt[i * self->dimension + j] * S[i] ;
  }

  //free(S); free(U); free(Vt);
  memset(&self->sketch[self->ell * self->dimension], 0, self->ell * self->dimension * sizeof(double));
}
