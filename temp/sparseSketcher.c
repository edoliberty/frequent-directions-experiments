#include "sparseSketcher.h"

void init_sparseSketcher(SparseSketcher* self, int ell, int dim ){
  self->class_name = "sparseSketcher";
  self->dimension = dim;
  self->ell = ell;
  self->m = 2*ell;
  self->sketch = (double*) malloc(sizeof(double) * (self->m) * dim);
  memset(self->sketch, 0, sizeof(double) * (self->m) * dim);
  init_sparseMatrix(&(self->buffer), dim, (ell+1)*dim );
  self->nnz_threshold = ell * dim;
}  


void append_to_sparseSketcher(SparseSketcher* self, SparseVector* sv){
  if((self->buffer).current_nnz >= self->nnz_threshold || (self->buffer).nextRow >= self->dimension)
    rotate_sparseSketcher(self);
  append_to_sparseMatrix(&(self->buffer), sv);
}

void rotate_sparseSketcher(SparseSketcher *self){
  /*
  (self->buffer).current_nnz = 0;
  (self->buffer).nextRow = 0;
  (self->buffer).pointer = 0;    
  */
  sparseShrink(self);
  denseShrink(self);
}

void get_sparseSketch(SparseSketcher *self){
  rotate_sparseSketcher(self);
}

void sparseShrink(SparseSketcher *self){
  if((self->buffer).nextRow > self->ell){
    //printf("IF \n");
    double* temp_vec = (double*) malloc(sizeof(double) * self->ell);
    double* temp_mat = (double*) malloc(sizeof(double) * self->ell * self->dimension);
    double* G = (double*) malloc(self->ell * self->dimension * sizeof(double));
    double* Z = (double*) malloc(self->ell * (self->buffer).nextRow * sizeof(double));
    int i,j;

    for(i=0; i < self->ell * self->dimension; i++)
      G[i] = ( (float)rand() / (float)(RAND_MAX) );

    blockPowerMethod(&(self->buffer), self->ell, 1, G, Z, temp_vec, temp_mat);
    free(temp_vec); 
    free(G);

    //computing P = ZtA, temp_mat is P
    transposeRightMult(&(self->buffer), self->ell, Z, temp_mat);
    free(Z);

    // svd(ZtA)
    double S[self->ell], U[self->ell * self->ell], Vt[self->dimension * self->ell];
    int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', self->ell, self->dimension, temp_mat, self->dimension, S, U, self->ell, Vt, self->dimension);
    free(temp_mat);


    // shrink S and compute S*Vt
    for(i=0; i < self->ell; i++){
      S[i] = sqrt( pow(S[i],2) - pow(S[self->ell-1],2) );
      for(j=0; j < self->dimension; j++)
	self->sketch[(self->ell + i) * self->dimension + j] = Vt[i * self->dimension + j] * S[i] ;
    }
  }else{ // self->buffer has atmost ell rows
    //printf("has atmost ell rows");
    int headptr = 0, ptr = 0;
    int rowIndex = (self->buffer).rows[headptr];
    int i = 0;

    while(ptr != (self->buffer).pointer){
      headptr = ptr;
      rowIndex = (self->buffer).rows[headptr];

      while(ptr != (self->buffer).pointer && (self->buffer).rows[ptr] == rowIndex){
	self->sketch[(self->ell + i) * self->dimension + (self->buffer).cols[ptr]] = (self->buffer).values[ptr];
	ptr ++;
      }
      i++;
    }
  }

  // reset buffer
  (self->buffer).current_nnz = 0;
  (self->buffer).nextRow = 0;
  (self->buffer).pointer = 0;    
}


void denseShrink(SparseSketcher* self){
  double S[2*self->ell], U[(2*self->ell) * (2*self->ell)], Vt[self->dimension * (2*self->ell)];
  int i, j, info;

  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', 2*self->ell, self->dimension, self->sketch, self->dimension, S, U, 2*self->ell, Vt, self->dimension);

  //for(i=0; i < self->ell; i++)
  
  for(i=0; i < self->ell; i++){
    S[i] = sqrt( pow(S[i],2) - pow(S[self->ell-1],2) );
    for(j=0; j < self->dimension; j++)
      self->sketch[i * self->dimension + j] = Vt[i * self->dimension + j] * S[i] ;
  }

  memset(&self->sketch[self->ell * self->dimension], 0, self->ell * self->dimension * sizeof(double));
}
