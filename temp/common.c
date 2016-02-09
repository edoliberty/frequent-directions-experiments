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

  int i,j;
  double *tau  = (double *) malloc( ell * sizeof( double ) );
  memset(tau, 0, ell * sizeof(double));

  lapack_int x = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, d, ell, G, ell, tau);
  x = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, d, ell, ell, G, ell, tau);
  free(tau);  
}


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
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'N', ell, d, mat, d, S, U, ell, Vt, d);
  return S[0];
}


/*
void readfile_transpose(char* filename){
  int n = 7769;
  int d = 26299;
  double* data = (double*) malloc(sizeof(double) * n * d);
  //read file
  

  //transpose
  

}
*/
/*
void readfile(){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen("./dataset/Reuters/1.ssv", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1) {
        printf("Retrieved line of length %zu :\n", read);
        printf("%s", line);
    }

    fclose(fp);
    if (line)
        free(line);
    exit(EXIT_SUCCESS);
}
*/
