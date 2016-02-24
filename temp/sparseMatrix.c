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
  for(int i=0; i < self->nextRow; i++)
    print_sparseVector(&(self->vectors[i]));
  
}


/*
 * computes A^TA.G where G is d*ell matrix
 * returns the result in G 
 * temp is 1*ell working memory
 * product is d*ell working memory
*/
void covMultiply_sparseMatrix (SparseMatrix* self, int d, int ell, double** G, double* temp, double** product){

  int rid;
  double val;
  SparseVector vec;

  for(int j=0; j < d * ell; j++)
    (*product)[j] = 0;
  
  for(int i = 0; i < self-> nextRow; i++){
    for(int j=0; j<ell; j++)
      temp[j] = 0;
    vec = self-> vectors[i];

    for(int j=0; j < vec.nnz; j++){
      rid = vec.cols[j] * ell;
      val = vec.values[j];

      for (int t=0; t<ell ; t++)
	temp[t] += (*G)[rid + t] * val;
    }
    
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

}


/* computes A*G for G being d*ell matrix
 * output is stored in double* product
 */
void leftMult (SparseMatrix* self, int ell, double* G, double* product){
  int itr = (self->nextRow) * ell;
  SparseVector vec;
  int rid, gidx;
  double val;

  for(int i=0; i < itr; i++)
    product[i] = 0;

  for(int i=0; i < self-> nextRow; i++){
    vec = self-> vectors[i];
    rid = i * ell;

    for (int t=0; t < vec.nnz ; t++){
      val = vec.values[t];
      gidx = vec.cols[t]*ell;

      for(int j=0; j < ell; j++){
	product[rid + j] += val * G[gidx + j];
      } 
    }
  }

}


/* computes Gt*A
 * G has ell columns
 * output is returned in double* product
 */
void transposeRightMult (SparseMatrix* self, int ell, double* G, double* product){
  int itr = (self->dimension) * ell;
  SparseVector vec;
  int rid, col;
  double val;

  for(int i=0; i<itr; i++)
    product[i] = 0;


  for(int i=0; i < self-> nextRow; i++){
    vec = self-> vectors[i];
    rid = i*ell;

    for(int j=0; j < vec.nnz; j++){
      val = vec.values[j];  
      col = vec.cols[j];

      for (int t=0; t<ell; t++){
	product[t * (self->dimension) + col] += G[rid + t] * val;
      }
    }
  }
}


void blockPowerMethod(SparseMatrix *self, int ell, double epsilon, double* G, double* lsv, double* temp_vec, double* temp_mat){
  int iterations = (int) ceil(1 * (log(self->dimension / epsilon) / epsilon));

  for(int i=0; i < iterations; i++){
    if(i % 10 == 0)
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
  SparseVector vec;

  for(int t=0; t < self->nextRow; t++){
    vec = self->vectors[t];
    for(int i=0; i< vec.nnz; i++){
      for(int j=0; j< vec.nnz; j++){
	elemIndex = vec.cols[i] * self-> dimension + vec.cols[j];
	val = (vec.values[i]) * (vec.values[j]);
	cov[elemIndex] += val;
      }
    }
  }
  
  return cov;
}

void densify_sparseMatrix(SparseMatrix* self, double* output){

  int rid;
  SparseVector vec;

  int itr = self->nextRow * self->dimension;
  for(int i=0; i<itr; i++)
    output[i] = 0;

  
  for(int t=0; t < self->nextRow; t++){
    vec = self->vectors[t];
    rid = t * self->dimension;

    for(int i=0; i < vec.nnz; i++)
      output[ rid + vec.cols[i] ] = vec.values[i];
  }
}


double computeCovErr(SparseMatrix* A, double* B, int ell, int d){
  double* AtA = getCovariance_sparseMatrix(A);
  double* BtB = getDenseCovariance(B, ell, d);
  subtract(AtA, BtB, d, d);
  return getSpectralNorm(AtA, d, d);
}

double computeRelCovErr(SparseMatrix* A, double* B, int ell, int d){
  double s = computeCovErr(A,B,ell,d);
  return s / A-> squaredFrob;

}

double topRank_cov(double* AtA, int d, int k){
  
  double* S = (double*) malloc(sizeof(double) * d);
  double* U = (double*) malloc(sizeof(double) * d * d);
  double* Vt = (double*) malloc(sizeof(double) * d * d);
  
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'N', d, d, AtA, d, S, U, d, Vt, d);

  free(U); free(Vt);
  double tailSquaredFrob = 0;

  for(int i = k; i < d ; i++)
    tailSquaredFrob += S[i];

  free(S);
  return tailSquaredFrob;

}


/* computes top rank k of A, returns it in Vt, returns tail norm of A too */
double topRank(SparseMatrix* A, int k){

  double* Adense = (double*) malloc(sizeof(double) * A->nextRow * A->dimension);
  densify_sparseMatrix(A, Adense);

  double* S = (double*) malloc(sizeof(double) * A->nextRow);
  double* U = (double*) malloc(sizeof(double) * A->nextRow * A->nextRow);
  double* Vt = (double*) malloc(sizeof(double) * A->dimension * A->dimension);
  
  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'N', A->nextRow, A->dimension, Adense, A->dimension, S, U, A->nextRow, Vt, A->dimension);

  free(U); free(Adense); free(Vt);
  int itr = min(A->nextRow, A->dimension);
  double tailSquaredFrob = 0;

  for(int i = k; i < itr ; i++)
    tailSquaredFrob += pow(S[i],2);
  
  free(S);
  return tailSquaredFrob;
}


double computeRelProjErr(SparseMatrix* A, double* B, int ell, int d, int k, double tailSquaredFrob){
 

  double projNorm = 0, projErr = 0;
  double projVec[k];
  SparseVector vec;
  int rid;


  double* S = (double*) malloc(sizeof(double) * 2 * ell);
  double* U = (double*) malloc(sizeof(double) * 4 * ell * ell);
  double* Vt = (double*) malloc(sizeof(double) * d * d);

  int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'A', 2*ell, d, B, d, S, U, 2*ell, Vt, d);
 

  for(int t=0; t< A->nextRow; t++){
    vec = A->vectors[t];
    projNorm = 0;

    for(int i=0; i<k; i++){
      projVec[i] = dotproduct(&vec, Vt, i, d);
      projNorm += pow(projVec[i],2);
    }

    projErr += (vec.squaredNorm - projNorm);
  }

  return projErr / tailSquaredFrob;
}


double dotproduct(SparseVector* sv, double* Vt, int rid, int dim){
  double dp = 0;
  for(int i=0; i < sv->nnz; i++ )
    dp += sv->values[i] * Vt[rid*dim + sv->cols[i]];

  return dp;
}
