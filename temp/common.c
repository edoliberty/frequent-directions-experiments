#include "common.h"


void print_two_dim(char* desc, double* mat, int m, int n) {
  int i, j;
  printf("%s \n",desc);
  for( i = 0; i < m; i++ ) {
    for(j=0; j < n; j++)
      printf( " %.10e", mat[i*n+j] );
    printf("\n ");
  }
}


void print_one_dim_double(char* desc, double* mat, int length){
  int i;
  printf("%s",desc);
  printf("[");
  for( i = 0; i < length; i++ )
    printf("%f , ",mat[i]);
  printf("]\n");
}

void print_one_dim_int(char* desc, int* mat, int length){
  int i;
  printf("%s",desc);
  printf("[");
  for( i = 0; i < length; i++ )
    printf("%d ,",mat[i]);
  printf("]\n");
}

void dot_product(double* C, int m, int n){
  int i,j;
  double dotproduct = 0;
  for( j = 0; j < n-1; j++ ) {
    dotproduct = 0;
    for( i = 0; i < m; i++ ) 
      dotproduct += C[i*n+j] * C[i*n+j+1];
    printf("dot product of columns %d and %d = %f", j, j+1, dotproduct);
    printf("\n");
  }
}

void normalizeVector(double* vec, int len){
  int i;
  int squaredNorm = 0;
  for(i=0; i<len; i++)
    squaredNorm += pow(vec[i],2);

  for(i=0; i<len; i++)
    vec[i] = vec[i] / sqrt(squaredNorm); 
  
}

void column_norm (double* C, int m, int n){
  double temp = 0;
  int i,j;
  for( j = 0; j < n; j++ ) {
    temp = 0;
    for( i = 0; i < m; i++ ) 
      temp += pow(C[i*n+j],2);
    printf("column %d norm = %f", j, temp);
    printf("\n");
  }
}

void printline(){
  printf("---------------------\n");
}

// computes AtA of dimensions d*d
double* getDenseCovariance(double* mat, int ell, int d){
  double* cov = (double*) malloc(sizeof(double) * d * d);
  memset(cov, 0, sizeof(double) * ell * d);
  
  int i,j,k;
  for(i=0; i<ell; i++)
    for(j=0; j<d; j++)
      for(k=0; k<d; k++){
	cov[j*d+k] += mat[i*d+j] * mat[i*d+k];
      }

  return cov;
}

// returns mat1 - mat2
void subtract(double* mat1, double* mat2, double* res){
  if (sizeof(mat1) != sizeof(mat2)){
    printf("dimensions of two matrices do not match");
    return;
  }
    
  int len = sizeof(res)/sizeof(res[0]);
  int i;
  for(i=0; i<len; i++)
    res[i] = mat1[i] - mat2[i];
}

// computes spectral norm of mat
double getSpectralNorm(double* mat, int ell, int d){
  double S[ell], U[1 * 1], Vt[1 * 1];
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'N', ell, d, mat, ell, S, U, ell, Vt, d);
  return S[0];
}

double computeCovErr(SparseMatrix* A, double* B, int ell, int d){
  double* AtA = getCovariance_sparseMatrix(A);
  double* BtB = getDenseCovariance(B, ell, d);
  subtract(AtA, BtB, AtA);
  double s = getSpectralNorm(AtA, d, d);
  return s;
}

double computeRelCovErr(SparseMatrix* A, double* B, int ell, int d){
  double s = computeCovErr(A,B,ell,d);
  return s / A-> squaredFrob;
}
