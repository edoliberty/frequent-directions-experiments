#include "sparseMatrix.h"


void init_sparseMatrix (SparseMatrix* self, int dim, int init_len){
  self->current_nnz = 0;
  self->dimension = dim;
  self->nextRow = 0;
  self->pointer = 0;
  self->arrays_length = init_len;
  self->rows = (int*) malloc(sizeof(int) * init_len);
  self->cols = (int*) malloc(sizeof(int) * init_len);
  self->values = (double*) malloc(sizeof(double) * init_len);
  self->squaredFrob = 0;
}

void extend_sparseMatrix(SparseMatrix* self){
  printf("in extend\n");
  int* info = realloc( self->rows, 2 * sizeof(int) * (self->arrays_length) );
  if(info != NULL)
    self->rows = info;

  info = realloc( self->cols, 2 * sizeof(int) * (self->arrays_length)  );
  if(info)
    self->cols = info;

  double* info2 = (double*) realloc( self->values, 2 * sizeof(double) * (self->arrays_length)  );
  if(info)
    self->values = info2;

  self->arrays_length *= 2;
}


void append_to_sparseMatrix (SparseMatrix* self, SparseVector* sv){
  // check if sv->dimension == self->dimension
  if (sv->dimension != self->dimension){
    printf("Vector and matrix have different dimensions %d and %d.\n", sv->dimension, self->dimension);
    return;
  }

  if(self->pointer + sv->nnz > self->arrays_length){// we need to double the size
    extend_sparseMatrix(self);
  }

  int i;
  for (i=0; i < sv->nnz ; i++){
    self->rows[self->pointer] = self->nextRow;
    self->cols[self->pointer] = sv->cols[i];
    self->values[self->pointer] = sv->values[i];
    self->pointer ++;
  }
  self->squaredFrob += sv->squaredNorm;
  self->nextRow ++;
  self->current_nnz += sv->nnz;
}


void print_sparseMatrix(SparseMatrix* self){
  int ptr = 0;
  int rowIndex = self->rows[ptr];

  while (ptr != self->pointer){
    if (self->rows[ptr] != rowIndex){
      rowIndex = self->rows[ptr];
      printf("\n");
    }
    printf("(%d,%d: %f)",self->rows[ptr],self->cols[ptr],self->values[ptr]);
    ptr ++;
  }
  printf("\n");    
}

/*
QR decomposition
Q is stored in G in row-wise format
R is not returned */
void qrDecomp(double* G, lapack_int d, lapack_int ell) {
  if(d == 1){ //G is a vector
    normalizeVector(G, ell);
    return;
  }

  int i,j;
  double *tau  = (double *) malloc( ell * sizeof( double ) );
  memset(tau, 0, ell * sizeof(double));

  lapack_int x = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, d, ell, G, ell, tau);
  x = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, d, ell, ell, G, ell, tau);
  free(tau);  
}

/*
computes A^TA.G 
stores the result in G */
void covarianceMult (SparseMatrix* self, int d, int ell, double** G, double* temp, double** product){
  int headptr = 0;
  int ptr = 0;
  int rowIndex = self->rows[headptr];
  int j,t;

  for(j=0; j < d * ell; j++)
    (*product)[j] = 0;
  
  // computation of A^TA.G
  while (ptr != self->pointer){
    headptr = ptr;
    rowIndex = self->rows[headptr];
    for (j=0; j<ell ; j++)
      temp[j] = 0;

    //A*G for one row of A
    while (ptr != self->pointer && self->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++)
	temp[t] += (*G) [self->cols[ptr]*ell + t] * self->values[ptr];
      ptr ++;
    }
    for (j=headptr; j<ptr ; j++){
      for (t=0; t<ell ; t++)
	(*product) [self->cols[j] * ell + t] += temp[t] * self->values[j];
    }
  }
  double* G_addr = *G;
  *G = *product;
  *product = G_addr;
  
}


// computes A*G, product is the output
void leftMult (SparseMatrix *self, int ell, double* G, double* product){
  int ptr = 0;
  int j,t;
  int rowIndex = self->rows[ptr];
  double* temp = (double*) malloc(sizeof(double) * ell);  

  for(j=0; j < (self->nextRow) * ell; j++)
    product[j] = 0;

  while (ptr != self->pointer){
    rowIndex = self->rows[ptr];
    for (j=0; j<ell ; j++)
      temp[j] = 0;

    //A*G for one row of A
    while (ptr != self->pointer && self->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++)
	temp[t] += G [self->cols[ptr]*ell + t] * self->values[ptr];
      ptr ++;
    }
    //product = that row of A*G
    for (t=0; t<ell ; t++)
      product[rowIndex * ell + t] = temp[t];
  }
  free(temp);
}


// computes G*A
double* rightMult (SparseMatrix *self, int ell, double* G){
  int ptr = 0;
  int rowIndex = self->rows[ptr];
  int j,t;

  // product is ell * d
  double* product = (double*) malloc( (self->dimension) * ell * sizeof(double));
  for(j=0; j < (self->dimension) * ell; j++)
    product[j] = 0;


  while (ptr != self->pointer){
    rowIndex = self->rows[ptr];
    while (ptr != self->pointer && self->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++)
	product[t * (self->dimension) +rowIndex] += G[t*(self->nextRow) + rowIndex] * self->values[ptr];

      ptr ++;
    }
  }
  return product;
}

// computes G^T*A
void transposeRightMult (SparseMatrix *self, int ell, double* G, double* product){
  int ptr = 0;
  int rowIndex = self->rows[ptr];
  int j,t;

  // product is ell * d
  for(j=0; j < (self->dimension) * ell; j++)
    product[j] = 0;
  
  while (ptr != self->pointer){
    rowIndex = self->rows[ptr];
    while (ptr != self->pointer && self->rows[ptr] == rowIndex){
      for (t=0; t<ell ; t++){
	product[ t * (self->dimension) +self->cols[ptr] ] += G[rowIndex * ell + t] * self->values[ptr];
      }
      ptr ++;
    }
  }
}


void blockPowerMethod(SparseMatrix *self, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat){
  int iterations = (int) ceil(1 * (log(self->dimension / epsilon) / epsilon));
  int i;
  //printf("numer of iteraiton is %d \n", iterations);

  for(i=0; i < iterations; i++){
    if(i % 3 == 0)
      qrDecomp(G, self->dimension, ell);
    covarianceMult(self, self->dimension, ell, &G, temp_vec, &temp_mat); 
  }

  // approx left singular vectors
  leftMult (self, ell, G, lsv);
  qrDecomp(lsv, self->nextRow, ell);

}


// returns covariance matrix, i.e. AtA
double* getCovariance_sparseMatrix(SparseMatrix* self){
  //printf("here\n");
  double* cov = (double*) malloc(sizeof(double) * self->dimension * self-> dimension);
  memset(cov, 0 , self->dimension * self-> dimension * sizeof(double));

  int headptr = 0, ptr = 0;
  int rowIndex = self-> rows[headptr];
  int i, j, elemIndex;
  double val;

  while(ptr != self-> pointer){
    //printf("in while\n");
    headptr = ptr;
    rowIndex = self-> rows[headptr];

    while(ptr != self-> pointer && self-> rows[ptr] == rowIndex)
      ptr ++;

    for(i = headptr; i < ptr; i++)
      for(j = headptr; j < ptr; j++){
	elemIndex = self-> cols[i] * self-> dimension + self-> cols[j];
	val = (self-> values[i]) * (self-> values[j]);
	cov[elemIndex] += val;
      }
  }
  return cov;
}

void densify_sparseMatrix(SparseMatrix* self, double* output){
  int headptr, ptr = 0;
  int self_rid;
  int output_rid = 0;

  memset(output, 0, sizeof(double) * self->nextRow * self->dimension);

  while(ptr != self->pointer){
    headptr = ptr;
    self_rid = self->rows[headptr];

    while (ptr != self->pointer && self->rows[ptr] == self_rid){
      output[output_rid * self->dimension + self->cols[ptr]] = self->values[ptr];
      ptr ++;
    }
    output_rid ++;
  }
}
