#include "common.h"

/*
 * QR decomposition of G
 * Q is returned in G, stored in the row-wise format
 * R is not returned 
*/
void qrDecomp(double* G, lapack_int d, lapack_int ell) {
  /* if G is a vector */
  if(d == 1){ 
    normalizeVector(G, ell);
    return;
  }

  double tau[ell];
  for(int i=0; i<ell; i++)
    tau[i] = 0;

  lapack_int x = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, d, ell, G, ell, tau);
  x = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, d, ell, ell, G, ell, tau);
}


void print_two_dim(char* desc, double* mat, int m, int n) {
  
  printf("%s \n",desc);
  for(int i = 0; i < m; i++ ) {
    for(int j=0; j < n; j++)
      printf( " %.10e", mat[i*n+j] );
    printf("\n ");
  }
}


void print_one_dim_double(char* desc, double* mat, int length){
  printf("%s",desc);
  printf("[");
  for(int i = 0; i < length; i++ )
    if (i < length - 1)
      printf("%f , ",mat[i]);
    else
      printf("%f ",mat[i]);
  printf("],\n");
}

void print_one_dim_int(char* desc, int* mat, int length){
  
  printf("%s",desc);
  printf("[");
  for(int i = 0; i < length; i++ )
    printf("%d ,",mat[i]);
  printf("]\n");
}

void dot_product(double* C, int m, int n){
  double dotproduct = 0;

  for(int j = 0; j < n-1; j++ ) {
    dotproduct = 0;
    for(int i = 0; i < m; i++ ) 
      dotproduct += C[i*n+j] * C[i*n+j+1];
    printf("dot product of columns %d and %d = %f", j, j+1, dotproduct);
    printf("\n");
  }
}

void normalizeVector(double* vec, int len){
  int squaredNorm = 0;

  for(int i=0; i<len; i++)
    squaredNorm += pow(vec[i],2);

  for(int i=0; i<len; i++)
    vec[i] = vec[i] / sqrt(squaredNorm); 
  
}

void column_norm (double* C, int m, int n){
  double temp = 0;
  
  for(int j = 0; j < n; j++ ) {
    temp = 0;
    for(int i = 0; i < m; i++ ) 
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

  for(int i=0; i<d; i++)
    for(int j=0; j<d; j++)
      cov[i*d+j] = 0;

  
 
  for(int i=0; i<ell; i++)
    for(int j=0; j<d; j++)
      for(int k=0; k<d; k++){
	cov[j*d+k] += mat[i*d+j] * mat[i*d+k];
      }

  return cov;
}

// returns mat1 - mat2
void subtract(double* mat1, double* mat2, int n, int d){
  if (sizeof(mat1) != sizeof(mat2)){
    printf("dimensions of two matrices do not match");
    return;
  }


  for(int i=0; i<n; i++)
    for(int j=0; j<d; j++)
      mat1[i*d+j] = mat1[i*d+j] - mat2[i*d+j];
}

// computes spectral norm of mat
double getSpectralNorm(double* mat, int ell, int d){

  double* S = (double*) malloc(sizeof(double) * ell);
  double* U = (double*) malloc(sizeof(double) * ell * 1);
  double* Vt = (double*) malloc(sizeof(double) * d * 1);

  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'N', ell, d, mat, d, S, U, ell, Vt, d);

  free(U); free(Vt);
  double sing = S[0];
  free(S);
  return sing;
}


void write_to_file(double* mat, int n, int d){
  printf("in Write to file \n");
  FILE* fp; 
  fp = fopen("CtC.txt","w");

  for(int i=0; i<n; i++){
    for(int j=0; j<d; j++){
      fprintf(fp, "%f ",mat[i*d +  j]);
    }
    fprintf(fp,"%s","\n");
  }

  fclose(fp);

}

