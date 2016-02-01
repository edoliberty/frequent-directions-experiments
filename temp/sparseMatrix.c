#include "sparseMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
//#include "common.c"

void init_matrix (SparseMatrix* sp, int dim){
  sp->current_nnz = 0;
  sp->dimension = dim;
  sp->nextRow = 0;
  sp->pointer = 0;
  sp->arrays_length = dim;
  sp->rows = (int*) malloc(sizeof(int) * dim);
  sp->cols = (int*) malloc(sizeof(int) * dim);
  sp->values = (double*) malloc(sizeof(double) * dim);
}

void print_sparseMatrix(SparseMatrix* sp){
  int ptr = 0;
  int rowIndex = sp->rows[ptr];

  while (ptr != sp->pointer){
    if (sp->rows[ptr] != rowIndex){
      rowIndex = sp->rows[ptr];
      printf("\n");
    }
    printf("(%d,%d: %f)",sp->rows[ptr],sp->cols[ptr],sp->values[ptr]);
    ptr ++;
  }
  printf("\n");    
}


void extend_matrix(SparseMatrix* sp){
  int* info = realloc( sp->rows, 2 * sizeof(int) * (sp->arrays_length) );
  if(info != NULL)
    sp->rows = info;

  info = realloc( sp->cols, 2 * sizeof(int) * (sp->arrays_length)  );
  if(info)
    sp->cols = info;

  double* info2 = (double*) realloc( sp->values, 2 * sizeof(double) * (sp->arrays_length)  );
  if(info)
    sp->values = info2;

  sp->arrays_length *= 2;
}

void append (SparseMatrix* sp, SparseVector* sv){
  // check if sv->dimension == sp->dimension
  if (sv->dimension != sp->dimension){
    printf("Vector and matrix have different dimensions %d and %d.\n", sv->dimension, sp->dimension);
    return;
  }

  if(sp->pointer + sv->nnz > sp->arrays_length){// we need to double the size
    extend_matrix(sp);
  }
  int i;

  for (i=0; i < sv->nnz ; i++){
    sp->rows[sp->pointer] = sp->nextRow;
    sp->cols[sp->pointer] = sv->cols[i];
    sp->values[sp->pointer] = sv->values[i];
    sp->pointer ++;
  }
  sp->nextRow ++;
  sp->current_nnz += sv->nnz;
}



// QR decomposition
// Q is stored in G in row-wise format
// R is not returned
void qrDecomp(double* G, lapack_int d, lapack_int ell) {
  int i,j;
  double *tau  = (double *) malloc( ell * sizeof( double ) );
  memset(tau, 0, ell * sizeof(double));

  lapack_int x = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, d, ell, G, ell, tau);
  x = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, d, ell, ell, G, ell, tau);
  free(tau);  
}


// computes A^TA.G 
// stores the result in G
void covarianceMult (SparseMatrix* sp, int d, int ell, double** G, double* temp, double** product){
  int headptr = 0;
  int ptr = 0;
  int rowIndex = sp->rows[headptr];
  int j,t;

  for(j=0; j < d * ell; j++)
    (*product)[j] = 0;
  
  // computation of A^TA.G
  while (ptr != sp->pointer){
    headptr = ptr;
    rowIndex = sp->rows[headptr];
    for (j=0; j<ell ; j++)
      temp[j] = 0;

    //A*G for one row of A
    while (ptr != sp->pointer && sp->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++)
	temp[t] += (*G) [sp->cols[ptr]*ell + t] * sp->values[ptr];
      ptr ++;
    }
    for (j=headptr; j<ptr ; j++){
      for (t=0; t<ell ; t++)
	(*product) [sp->cols[j] * ell + t] += temp[t] * sp->values[j];
    }
  }
  //free(temp);
  //free(*G);
  double* G_addr = *G;
  *G = *product;
  *product = G_addr;
  
}


// computes A*G, product is the output
void leftMult (SparseMatrix *sp, int ell, double* G, double* product){
  int ptr = 0;
  int j,t;
  int rowIndex = sp->rows[ptr];
  double* temp = (double*) malloc(sizeof(double) * ell);  

  for(j=0; j < (sp->nextRow) * ell; j++)
    product[j] = 0;

  while (ptr != sp->pointer){
    rowIndex = sp->rows[ptr];
    for (j=0; j<ell ; j++)
      temp[j] = 0;

    //A*G for one row of A
    while (ptr != sp->pointer && sp->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++)
	temp[t] += G [sp->cols[ptr]*ell + t] * sp->values[ptr];
      ptr ++;
    }
    //product = that row of A*G
    for (t=0; t<ell ; t++)
      product[rowIndex * ell + t] = temp[t];
  }
  free(temp);
}


// computes G*A
double* rightMult (SparseMatrix *sp, int ell, double* G){
  int ptr = 0;
  int rowIndex = sp->rows[ptr];
  int j,t;

  // product is ell * d
  double* product = (double*) malloc( (sp->dimension) * ell * sizeof(double));
  for(j=0; j < (sp->dimension) * ell; j++)
    product[j] = 0;


  while (ptr != sp->pointer){
    rowIndex = sp->rows[ptr];
    while (ptr != sp->pointer && sp->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++)
	product[t * (sp->dimension) +rowIndex] += G[t*(sp->nextRow) + rowIndex] * sp->values[ptr];

      ptr ++;
    }
  }
  return product;
}

// computes G^T*A
void transposeRightMult (SparseMatrix *sp, int ell, double* G, double* product){
  int ptr = 0;
  int rowIndex = sp->rows[ptr];
  int j,t;

  // product is ell * d
  for(j=0; j < (sp->dimension) * ell; j++)
    product[j] = 0;
  
  while (ptr != sp->pointer){
    rowIndex = sp->rows[ptr];
    while (ptr != sp->pointer && sp->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++){
	product[ t * (sp->dimension) +sp->cols[ptr] ] += G[rowIndex * ell + t] * sp->values[ptr];
      }
      ptr ++;
    }
  }
}


void blockPowerMethod(SparseMatrix *sp, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat){

  int iterations = (int) ceil(10 * (log(sp->dimension / epsilon) / epsilon));
  int i;

  for(i=0; i < iterations; i++){
    qrDecomp(G, sp->dimension, ell);
    covarianceMult(sp, sp->dimension, ell, &G, temp_vec, &temp_mat); 
  }
  // approx right singular vectors
  leftMult (sp, ell, G, lsv);
  qrDecomp(lsv, sp->nextRow, ell);

}


/*
double* sparseShrink(SparseMatrix *sp, int ell){

  double* temp_vec = (double*) malloc(sizeof(double) * ell);
  double* temp_mat = (double*) malloc(sizeof(double) * ell * sp->dimension);
  double* G = (double*) malloc(ell * sp->dimension * sizeof(double));
  double* Z = (double*) malloc(ell * sp->nextRow * sizeof(double));
  int i,j;

  for(i=0; i < ell * sp->dimension; i++)
    G[i] = ( (float)rand() / (float)(RAND_MAX) );

  blockPowerMethod(sp, ell, 0.25, G, Z, temp_vec, temp_mat);
  free(temp_vec); 
  free(G);

  //temp_mat becomes ZtA
  transposeRightMult(sp, ell, Z, temp_mat);
  free(Z);

  // svd(ZtA)
  double S[ell], U[ell * ell], Vt[sp->dimension * ell];
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S',  ell, sp->dimension, temp_mat, sp->dimension, S, U, ell, Vt, sp->dimension);

  // shrink S
  for(i=0; i<ell; i++)
    S[i] = sqrt( pow(S[i],2) - pow(S[ell-1],2) );
  
  // compute S*Vt
  for(i=0; i<ell; i++)
    for(j=0; j<sp->dimension; j++)
      temp_mat[i*sp->dimension + j] = Vt[i*sp->dimension + j] * S[i] ;


  return temp_mat;
}
*/
