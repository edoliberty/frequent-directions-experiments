#include "sparseMatrix.h"


void init_sparseMatrix (SparseMatrix* self, int dim, int len){
  self->current_nnz = 0;
  self->dimension = dim;
  self->nextRow = 0;
  self->squaredFrob = 0;

  self->vectors = (SparseVector*) malloc(sizeof(SparseVector) * len);
}


void append_to_sparseMatrix (SparseMatrix* self, SparseVector* sv){
  if (sv->dimension != self->dimension){
    printf("Vector and matrix have different dimensions %d and %d.\n", sv->dimension, self->dimension);
    return;
  }

  self->vectors[self->nextRow] = *sv;
  self->nextRow ++;
  self->squaredFrob += sv->squaredNorm;
  self->current_nnz += sv->nnz;
}


void print_sparseMatrix(SparseMatrix* self){
  int i;
  for(i=0; i < self->nextRow; i++)
    print_sparseVector(&(self->vectors[i]));
  
}


/*
 * computes A^TA.G where G is d*ell matrix
 * returns the result in G 
 * temp is 1*ell working memory
 * product is d*ell working memory
*/
void covMultiply_sparseMatrix (SparseMatrix* self, int d, int ell, double** G, double* temp, double** product){
  int i,j,t;
  SparseVector vec;

  for(j=0; j < d * ell; j++)
    (*product)[j] = 0;
  
  for(i = 0; i < self-> nextRow; i++){
    memset(temp, 0, sizeof(double) * ell);
    vec = self-> vectors[i];

    /* A*G for one row of A */
    for (t=0; t<ell ; t++)
      for(j=0; j < vec.nnz; j++)
	temp[t] += vec.values[j] * (*G)[vec.cols[j]*ell + t];

    /* At*temp for corresponding column of At */
    for(j=0; j < vec.nnz; j++)
      for (t=0; t<ell ; t++)
	(*product) [vec.cols[j] * ell + t] += temp[t] * (vec.values[j]);
  }

  double* G_addr = *G;
  *G = *product;
  *product = G_addr;

}


/* computes A*G for G being d*ell matrix
 * output is stored in double* product
 */
void leftMult (SparseMatrix* self, int ell, double* G, double* product){
  int i,j,t;
  memset(product, 0, sizeof(double) * (self->nextRow) * ell);
  SparseVector vec;

  for(i=0; i < self-> nextRow; i++){
    vec = self-> vectors[i];
    for (t=0; t<ell ; t++)
      for(j=0; j < vec.nnz; j++)
	product[i * ell + t] += vec.values[j] * G[vec.cols[j]*ell + t];
  } 
}


/* computes Gt*A
 * G has ell columns
 * output is returned in double* product
 */
void transposeRightMult (SparseMatrix* self, int ell, double* G, double* product){
  int j,i,t;

  // product is ell * d
  memset(product, 0, sizeof(double) * (self->dimension) * ell);
  SparseVector temp;

  for(i=0; i < self-> nextRow; i++){
    temp = self-> vectors[i];
    for (t=0; t<ell ; t++)
      for(j=0; j < temp.nnz; j++)
	product[ t * (self->dimension) + temp.cols[j] ] += G[i * ell + t] * temp.values[j];     
  }
}


void blockPowerMethod(SparseMatrix *self, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat){
  int iterations = (int) ceil(1 * (log(self->dimension / epsilon) / epsilon));
  int i;
  for(i=0; i < iterations; i++){
    if(i % 3 == 0)
      qrDecomp(G, self->dimension, ell);
    covMultiply_sparseMatrix(self, self->dimension, ell, &G, temp_vec, &temp_mat); 
  }

  // approx left singular vectors
  leftMult (self, ell, G, lsv);
  qrDecomp(lsv, self->nextRow, ell);

}


/* returns covariance matrix, i.e. AtA
 */
double* getCovariance_sparseMatrix(SparseMatrix* self){

  double* cov = (double*) malloc(sizeof(double) * self->dimension * self-> dimension);
  memset(cov, 0 , self->dimension * self-> dimension * sizeof(double));

  int t, i, j, elemIndex;
  double val;
  SparseVector temp;

  for(t=0; t < self->nextRow; t++)
    temp = self->vectors[t];
    for(i=0; i< temp.nnz; i++)
      for(j=0; j< temp.nnz; j++){
	elemIndex = temp.cols[i] * self-> dimension + temp.cols[j];
	val = (temp.values[i]) * (temp.values[j]);
	cov[elemIndex] += val;
      }

  return cov;
}

void densify_sparseMatrix(SparseMatrix* self, double* output){

  int t,i;
  memset(output, 0, sizeof(double) * self->nextRow * self->dimension);
  SparseVector temp;
  
  for(t=0; t < self->nextRow; t++){
    temp = self->vectors[t];
    for(i=0; i<temp.nnz; i++)
      output[ t * self->dimension + temp.cols[i] ] = temp.values[i];
  }
}

