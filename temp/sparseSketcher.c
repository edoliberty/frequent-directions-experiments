#include "sparseSketcher.h"

void init_sketcher(SparseSketcher* ss, int ell, int dim ){
  ss->class_name = "sparseSketcher";
  ss->dimension = dim;
  ss->ell = ell;
  ss->sketch = (double*) malloc(sizeof(double) * 2 * ell * dim);
  init_matrix(&(ss->buffer), dim );
  ss->nnz_threshold = ell * dim;
  memset(ss->sketch, 0, sizeof(double)*2*ell*dim);
}

void rotate(SparseSketcher *ss){
  sparseShrink(ss);
  denseShrink(ss);
}

void get_sketch(SparseSketcher *ss){
  rotate(ss);
}


void append_row(SparseSketcher* ss, SparseVector* sv){
  if((ss->buffer).current_nnz >= ss->nnz_threshold || (ss->buffer).nextRow >= ss->dimension)
    rotate(ss);
  append(&(ss->buffer), sv);
}

void sparseShrink(SparseSketcher *ss){

  double* temp_vec = (double*) malloc(sizeof(double) * ss->ell);
  double* temp_mat = (double*) malloc(sizeof(double) * ss->ell * ss->dimension);
  double* G = (double*) malloc(ss->ell * ss->dimension * sizeof(double));
  double* Z = (double*) malloc(ss->ell * (ss->buffer).nextRow * sizeof(double));
  int i,j;

  for(i=0; i < ss->ell * ss->dimension; i++)
    G[i] = ( (float)rand() / (float)(RAND_MAX) );

  blockPowerMethod(&(ss->buffer), ss->ell, 0.25, G, Z, temp_vec, temp_mat);
  free(temp_vec); 
  free(G);

  //temp_mat becomes ZtA
  transposeRightMult(&(ss->buffer), ss->ell, Z, temp_mat);
  free(Z);

  // svd(ZtA)
  double S[ss->ell], U[ss->ell * ss->ell], Vt[ss->dimension * ss->ell];
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', ss->ell, ss->dimension, temp_mat, ss->dimension, S, U, ss->ell, Vt, ss->dimension);
  free(temp_mat);

  // shrink S
  for(i=0; i < ss->ell; i++){
    S[i] = sqrt( pow(S[i],2) - pow(S[ss->ell-1],2) );
    //printf("S: %f ,",S[i]);
  }
  //printf("\n");
  // compute S*Vt
  for(i=0; i < ss->ell; i++)
    for(j=0; j < ss->dimension; j++)
      ss->sketch[(ss->ell + i) * ss->dimension + j] = Vt[i * ss->dimension + j] * S[i] ;

  // reset buffer
  (ss->buffer).current_nnz = 0;
  (ss->buffer).nextRow = 0;
  (ss->buffer).pointer = 0;    
}


void denseShrink(SparseSketcher* ss){
  double S[2*ss->ell], U[(2*ss->ell) * (2*ss->ell)], Vt[ss->dimension * (2*ss->ell)];
  int i, j, info;

  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', 2*ss->ell, ss->dimension, ss->sketch, ss->dimension, S, U, 2*ss->ell, Vt, ss->dimension);
  for(i=0; i < ss->ell; i++)
    S[i] = sqrt( pow(S[i],2) - pow(S[ss->ell-1],2) );
  
  for(i=0; i < ss->ell; i++)
    for(j=0; j < ss->dimension; j++)
      ss->sketch[i * ss->dimension + j] = Vt[i * ss->dimension + j] * S[i] ;

  memset(&ss->sketch[ss->ell * ss->dimension], 0, ss->ell * ss->dimension * sizeof(double));
}
