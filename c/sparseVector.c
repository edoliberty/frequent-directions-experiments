#include "sparseVector.h"


void init_sparseVector(SparseVector* self, int dim, int cols[], double vals[], int nnz){
  self-> nnz = nnz;
  self-> dimension = dim;
  self-> cols = (int*) malloc(sizeof(int) * self-> nnz);
  self-> values = (double*) malloc(sizeof(double) * self-> nnz);
  self-> squaredNorm = 0;
  
  for (int i=0; i < self-> nnz; i++){
    self-> cols[i] = cols[i];
    self-> values[i] = vals[i];
    self-> squaredNorm += pow(vals[i] , 2);
  }
}


/* it generates a vector of dim dimension, 
   with only nnz non-zeros
   first jlen columns have threshold_prob probability of getting a non-zero
   non-zeros are picked from [-10, 10] uniformly at random
 */
void skew_init_sparseVector(SparseVector* self, int dim, int nnz, int jlen, double threshold_prob){
  self-> dimension = dim;  
  self-> nnz = nnz;
  self-> cols = (int*) malloc(sizeof(int) * self-> nnz);
  self-> values = (double*) malloc(sizeof(double) * self-> nnz);
  self-> squaredNorm = 0;

  double randomVal;  
  int flag, col_id, t;

  for (int i=0; i < self-> nnz; i++){
    randomVal = rand()/(RAND_MAX+1.0); 

    if(randomVal < threshold_prob){ // goes to first "jlen" columns
      col_id = (int) rand() % jlen;
      flag = 1;

      while (flag == 1){
	for (t=0; t < i; t++)
	  if (col_id == self-> cols[t])
	    break;
	if (t == i)
	  flag = 0; 
	else
	  col_id = rand() % jlen;
      }
    }

    else{// goes to the rest of columns
      col_id = jlen + (int) rand() % (dim-jlen);
      flag = 1;

      while (flag == 1){
	for (t=0; t < i; t++)
	  if (col_id == self-> cols[t])
	    break;
	if (t == i)
	  flag = 0; 
	else
	  col_id = jlen + (int) rand() % (dim-jlen);
      }
    }
    self-> cols[i] = col_id;
    int tempr = 2 * (rand()%2) - 1;
    self-> values[i] = tempr * (int)ceil( ((double)rand()/(double)(RAND_MAX)) * 10);
  
    self-> squaredNorm += pow(self-> values[i] , 2);
  }
}


void random_init_sparseVector(SparseVector* self, int dim, int nnz){
  self-> dimension = dim;  
  self-> nnz = nnz;
  self-> cols = (int*) malloc(sizeof(int) * self-> nnz);
  self-> values = (double*) malloc(sizeof(double) * self-> nnz);
  self-> squaredNorm = 0;
  int i;

  for (i=0; i < self-> nnz; i++){
    double newly_gen = rand() % dim;
    int flag = 1;
    int j= 0;

    while (flag == 1){
      for (j=0; j < i; j++)
	if (newly_gen == self-> cols[j])
	  break;
      if (j == i)
	flag = 0; 
      else
	newly_gen = rand() % dim;
    }
    self-> cols[i] = newly_gen;
    self-> values[i] = (int)ceil( ((double)rand()/(double)(RAND_MAX)) * 10);
    self-> squaredNorm += pow(self-> values[i] , 2);
  }
}

void print_sparseVector(SparseVector* self){
  
  for (int i=0; i< self-> nnz; i++)
    printf("(%d, %.2f)", self-> cols[i], self-> values[i]  );
  printf("\n");
}


double* densify_sparseVector(SparseVector* self){
  double* vec = (double*) malloc(sizeof(double) * self->dimension);
  
  for(int i=0; i < self->dimension ; i++)
    vec[i] = 0;
  for(int i=0; i < self->nnz; i++)
    vec[self->cols[i]] = self->values[i];
  return vec;
}

