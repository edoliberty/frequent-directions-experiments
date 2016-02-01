#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>
#include <time.h>
#include <stdint.h>


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

/*
void write_to_file(double* mat, int m, int n, char* filename){
  FILE* file = fopen(filename, "a" );
  fputs(str,file);
  fclose(file);
}
*/
