#include "sparseMatrix.h"


void init_sparseMatrix (SparseMatrix* self, int dim, int len){
  self->current_nnz = 0;
  self->dimension = dim;
  self->nextRow = 0;
  self->squaredFrob = 0;
  self->vectors = (SparseVector*) malloc(sizeof(SparseVector) * len);
}


void append_to_sparseMatrix (SparseMatrix* self, SparseVector* sv){
  self->vectors[self->nextRow] = *sv;
  self->nextRow ++;
  self->squaredFrob += sv->squaredNorm;
  self->current_nnz += sv->nnz;
}


void print_sparseMatrix(SparseMatrix* self){
  int itr = self->nextRow;
  for(int i=0; i < itr; i++)
    print_sparseVector(&(self->vectors[i]));
  
}


/*
 * computes A^TA.G where G is d*ell matrix
 * returns the result in G 
 * temp is 1*ell working memory
 * product is d*ell working memory
*/
void covMultiply_sparseMatrix (SparseMatrix* self, int d, int ell, double** G, double* temp, double** product){
  SparseVector vec;
  int dell = d*ell;
  double val;
  double* g = (double*) malloc(sizeof(double) * ell);
  int rid;

  for(int j=0; j < dell; j++)
    (*product)[j] = 0;
  
  for(int i = 0; i < self-> nextRow; i++){
    memset(temp, 0, sizeof(double) * ell);
    vec = self-> vectors[i];
    
    /* A*G for one row of A */    
    /*
    for (int t=0; t<ell ; t++)
      for(int j=0; j < vec.nnz; j++)
        temp[t] += vec.values[j] * (*G)[vec.cols[j]*ell + t];
    */

    
    for(int j=0; j < vec.nnz; j++){
      val = vec.values[j];
      rid = vec.cols[j] * ell;
      memcpy(g, &((*G)[rid]), sizeof(double) * ell); 
      for (int t=0; t<ell; t++){
	temp[t] += val * g[t];
      }      
    }


    /* At*temp for corresponding column of At */
    for(int j=0; j < vec.nnz; j++){
      rid = vec.cols[j] * ell;
      val = vec.values[j];
      for (int t=0; t<ell; t++)
	(*product) [rid + t] += temp[t] * val;
    }
  }

  double* G_addr = *G;
  *G = *product;
  *product = G_addr;
  free(g);
}


/* computes A*G for G being d*ell matrix
 * output is stored in double* product
 */
void leftMult (SparseMatrix* self, int ell, double* G, double* product){
  memset(product, 0, sizeof(double) * (self->nextRow) * ell);
  SparseVector vec;
  int rid;
  double* g = (double*) malloc(sizeof(double) * ell);
  double val;

  for(int i=0; i < self-> nextRow; i++){
    vec = self-> vectors[i];
    rid = i * ell;
    for (int t=0; t<vec.nnz ; t++){
      memcpy(g, &(G[vec.cols[t]*ell]), sizeof(double) * ell);  // works?
      val = vec.values[t];

      for(int j=0; j < ell; j++){
	product[rid + j] += val * g[j];
      } 
    }
  }
}


/* computes Gt*A
 * G has ell columns
 * output is returned in double* product
 */
void transposeRightMult (SparseMatrix* self, int ell, double* G, double* product){

  memset(product, 0, sizeof(double) * (self->dimension) * ell);
  SparseVector vec;
  double* g = (double*) malloc(sizeof(double) * ell);
  int rid, col;
  double val;

  for(int i=0; i < self-> nextRow; i++){
    vec = self-> vectors[i];
    rid = i*ell;
    memcpy(g, &(G[rid]), sizeof(double) * ell);  // works?

    for(int j=0; j < vec.nnz; j++){
      val = vec.values[j];  
      col = vec.cols[j];
      for (int t=0; t<ell; t++){
	product[t * (self->dimension) + col] += g[t] * val;
      }
    }
  }
  free(g);
}


void blockPowerMethod(SparseMatrix *self, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat){
  int iterations = (int) ceil(1 * (log(self->dimension / epsilon) / epsilon));
  
  for(int i=0; i < iterations; i++){
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

  int elemIndex;
  double val;
  SparseVector temp;

  for(int t=0; t < self->nextRow; t++)
    temp = self->vectors[t];
    for(int i=0; i< temp.nnz; i++)
      for(int j=0; j< temp.nnz; j++){
	elemIndex = temp.cols[i] * self-> dimension + temp.cols[j];
	val = (temp.values[i]) * (temp.values[j]);
	cov[elemIndex] += val;
      }

  return cov;
}

void densify_sparseMatrix(SparseMatrix* self, double* output){

  memset(output, 0, sizeof(double) * self->nextRow * self->dimension);
  SparseVector temp;
  
  for(int t=0; t < self->nextRow; t++){
    temp = self->vectors[t];
    for(int i=0; i < temp.nnz; i++)
      output[ t * self->dimension + temp.cols[i] ] = temp.values[i];
  }
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


double computeRelProjErr(SparseMatrix* A, double* B, int ell, int d, int k){

  double* Adense = (double*) malloc(sizeof(double) * A->nextRow * A->dimension);
  densify_sparseMatrix(A, Adense);

  double* S = (double*) malloc(sizeof(double) * A->nextRow);
  double* U = (double*) malloc(sizeof(double) * A->nextRow * A->nextRow);
  double* Vt = (double*) malloc(sizeof(double) * A->dimension * A->dimension);
  
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', A->nextRow, A->dimension, Adense, A->dimension, S, U, A->nextRow, Vt, A->dimension);

  free(U); free(Adense);
  double tailSquaredFrob = 0;
  
  int itr = min(A->nextRow,A->dimension);
  for(int i = k; i < itr ; i++)
    tailSquaredFrob += pow(S[i],2);
  free(S);

  double projNorm = 0, projErr = 0;
  double* projVec = (double*) malloc(sizeof(double) * k);
  SparseVector vec;

  for(int t=0; t< A->nextRow; t++){
    vec = A->vectors[t];
    projNorm = 0;
    memset(projVec, 0, sizeof(double) * k);
    for(int i=0; i<k; i++){
      for(int j=0; j<vec.nnz; j++){
	projVec[i] += vec.values[j] * Vt[i*A->dimension + vec.cols[j]];
      }
      projNorm += pow(projVec[i],2);
    }
    projErr += vec.squaredNorm - projNorm;
  }

  free(Vt);
  return projErr / tailSquaredFrob;
}
