#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "randomSum.h"
//#include "randomProjection.h"
#include "sparseSketcher.h"
//#include "rowSampler.h"
#include "common.h"
#include "frequentDirections.h"

int main2(){
  int n = 10000;
  int dim = 2000;
  int ell = 100;
  int k = 10;
  int exp_no = 1;
  int var_set[] = {0.01 * dim, 0.05 * dim, 0.1 * dim, 0.5 * dim}; //sparsity 1% , 10%, ....
  int i, j;
  double start, end, cpu_time_used;
  SparseMatrix A;
  SparseVector arr[n];
  SparseSketcher sfd;
  FrequentDirections fd;
  double sfd_cov_err[exp_no], sfd_proj_err[exp_no], sfd_time[exp_no];
  double fd_cov_err[exp_no], fd_proj_err[exp_no], fd_time[exp_no];

  init_sparseSketcher(&sfd, ell, dim);
  init_fd(&fd, ell, dim);

 for(i=0; i<exp_no; i++){
    sfd_time[i] = 0;
    fd_time[i] = 0;
    init_sparseMatrix(&A, dim, n);

    for (j=0; j < n; j++){
      random_init_sparseVector(&arr[j], dim, var_set[i]);
      append_to_sparseMatrix(&A, &arr[j]);

      // SFD
      start = clock();
      append_to_sparseSketcher(&sfd, &arr[j]);
      sfd_time[i] += (double) (clock() - start);

      //FD
      start = clock();
      append_to_fd(&fd, &arr[j]);
      fd_time[i] += (double) (clock() - start);

    }
    //SFD
    start = clock();
    get_sparseSketch(&sfd);
    sfd_time[i] += (double) (clock() - start);
    sfd_time[i] = sfd_time[i] / CLOCKS_PER_SEC;
    sfd_cov_err[i] = computeRelCovErr(&A, sfd.sketch, ell, dim);
    sfd_proj_err[i] = computeRelProjErr(&A, sfd.sketch, ell, dim, k);

    //FD
    start = clock();
    get_fdSketch(&fd);
    fd_time[i] += (double) (clock() - start);
    fd_time[i] = fd_time[i] / CLOCKS_PER_SEC;
    fd_cov_err[i] = computeRelCovErr(&A, fd.sketch, ell, dim);
    fd_proj_err[i] = computeRelProjErr(&A, fd.sketch, ell, dim, k);
  }

  
  print_one_dim_double("FD COV ERR= ", fd_cov_err, exp_no);
  print_one_dim_double("SFD COV ERR= ", sfd_cov_err, exp_no);

  print_one_dim_double("FD PROJ ERR= ", fd_proj_err, exp_no);
  print_one_dim_double("SFD PROJ ERR= ", sfd_proj_err, exp_no);
  
  print_one_dim_double("FD TIME= ", fd_time, exp_no);
  print_one_dim_double("SFD TIME= ", sfd_time, exp_no);

}


int main_tt(){
  int n = 10000;
  int dim = 1000;
  int ell = 100;
  int k = 10;
  int ell_len = 1;
  int ell_set[] = {100};
  int i, ell_counter;
  double start, end, cpu_time_used;

  SparseMatrix A;
  SparseSketcher sfd;
  SparseVector arr[n];
  //RandomSum hashing;
  //RandomProjection proj;
  //RowSampler sampler;
  FrequentDirections fd;

  init_sparseMatrix(&A, dim, n);
  int nnz =  (int) ceil(0.01*dim);
  printf("%d ", nnz);

  for (i=0; i < n; i++){
    random_init_sparseVector(&arr[i], dim, nnz);
    append_to_sparseMatrix(&A, &arr[i]);
  }

  printf("after initialization of A\n");
  ///print_sparseMatrix(&A);
  //double* dense = (double*) malloc(sizeof(double) * n * dim);
  //densify_sparseMatrix(&A, dense);
  //print_two_dim("dense = ",dense,n,dim);

  // FD
  double fd_proj_err[ell_len];
  double fd_cov_err[ell_len];
  double fd_time[ell_len];
  for(ell_counter = 0; ell_counter < ell_len; ell_counter++){
    printf("FD Itr %d\n", ell_counter);
    ell = ell_set[ell_counter];
    init_fd(&fd, ell, dim);
    
    start = clock();
    for (i=0; i < n; i++){
      append_to_fd(&fd, &arr[i]);
    }
    get_fdSketch(&fd);
    end = clock();

    fd_time[ell_counter] = ((double) (end - start)) / CLOCKS_PER_SEC;
    //fd_cov_err[ell_counter] = computeRelCovErr(&A, fd.sketch, ell, dim);
    //fd_proj_err[ell_counter] = computeRelProjErr(&A, fd.sketch, ell, dim, k);    
  }

  

  
  // SFD
  double sfd_cov_err[ell_len], sfd_proj_err[ell_len];
  double sfd_time[ell_len];
  for(ell_counter = 0; ell_counter < ell_len; ell_counter++){
    printf("SFD Itr %d\n", ell_counter);
    ell = ell_set[ell_counter];
    init_sparseSketcher(&sfd, ell, dim);
    
    start = clock();
    for (i=0; i < n; i++){
      append_to_sparseSketcher(&sfd, &arr[i]);
    }
    get_sparseSketch(&sfd);
    end = clock();
    sfd_time[ell_counter] = ((double) (end - start)) / CLOCKS_PER_SEC;
    //sfd_cov_err[ell_counter] = computeRelCovErr(&A, sfd.sketch, ell, dim);
    //sfd_proj_err[ell_counter] = computeRelProjErr(&A, sfd.sketch, ell, dim, k);
  }

  /*
  print_one_dim_double("FD PROJ ERR= ", fd_proj_err, ell_len);
  print_one_dim_double("SFD PROJ ERR= ", sfd_proj_err, ell_len);

  print_one_dim_double("FD COV ERR= ", fd_cov_err, ell_len);
  print_one_dim_double("SFD COV ERR= ", sfd_cov_err, ell_len);
  */

  print_one_dim_double("FD TIME= ", fd_time, ell_len);
  print_one_dim_double("SFD TIME= ", sfd_time, ell_len);
  

  /*
  // HASHING
  double hashing_cov_err[4];
  double hashing_time[4];
  for(ell_counter = 0; ell_counter < 4; ell_counter++){
    printf("HASHING Itr %d\n", ell_counter);
    ell = ell_set[ell_counter];
    init_randomSum(&hashing, ell, dim);

    start = clock();
    for (i=0; i < n; i++){
      append_randomSum(&hashing, &arr[i]);
    }
    end = clock();
    hashing_time[ell_counter] = ((double) (end - start)) / CLOCKS_PER_SEC;

    hashing_cov_err[ell_counter] = computeRelCovErr(&A, hashing.sketch, ell, dim);
  }
  */
  /*
  // RANDOM PROJECTION
  double proj_cov_err[4];
  double proj_time[4];
  for(ell_counter = 0; ell_counter < 4; ell_counter++){
    printf("PROJ Itr %d\n", ell_counter);
    ell = ell_set[ell_counter];
    init_randomProj(&proj, ell, dim);

    start = clock();
    for (i=0; i < n; i++){
      append_randomProj(&proj, &arr[i]);
    }
    end = clock();
    proj_time[ell_counter] = ((double) (end - start)) / CLOCKS_PER_SEC;

    proj_cov_err[ell_counter] = computeRelCovErr(&A, proj.sketch, ell, dim);
  }
  */
  /*
  // ROW SAMPLER
  double sampler_cov_err[4];
  double sampler_time[4];
  for(ell_counter = 0; ell_counter < 4; ell_counter++){
    printf("SAMPLING Itr %d\n", ell_counter);
    ell = ell_set[ell_counter];
    init_rowSampler(&sampler, ell, dim);

    start = clock();
    for (i=0; i < n; i++){
      append_rowSampler(&sampler, &arr[i]);
    }
    set_rowSamplerSketch(&sampler);
    end = clock();
    sampler_time[ell_counter] = ((double) (end - start)) / CLOCKS_PER_SEC;

    sampler_cov_err[ell_counter] = computeRelCovErr(&A, sampler.sketch, ell, dim);
  }
  */
  /*
  print_one_dim_double("FD COV ERR= ", fd_cov_err, 4);
  print_one_dim_double("SFD COV ERR= ", sfd_cov_err, 4);
  print_one_dim_double("HASHING COV ERR= ", hashing_cov_err, 4);  
  print_one_dim_double("PROJ COV ERR= ", proj_cov_err, 4);  
  print_one_dim_double("SAMPLER COV ERR= ", sampler_cov_err, 4);  

  print_one_dim_double("FD TIME= ", fd_time, 4);
  print_one_dim_double("SFD TIME= ", sfd_time, 4);
  print_one_dim_double("HASHING TIME= ", hashing_time, 4);  
  print_one_dim_double("PROJ TIME= ", proj_time, 4);  
  print_one_dim_double("SAMPLER TIME= ", sampler_time, 4);  
  */
  /*
  RowSampler A;
  init_rowSampler(&A, ell, dim);

  printline();
  start = clock();
  for (i=0; i < n; i++)
    append_rowSampler(&A, &arr[i]);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Done in %f seconds\n",cpu_time_used);

  set_sketch(&A);
  print_two_dim("A=",A.sketch, A.ell, A.dimension);
  */
}



void blockPower_vs_lapacke(){
  printf("in blcokpower_vs_lkapacke\n");
  int n = 10000;//50000;
  int dim = 20000;//2000;
  int ell = 10000;//50000;
  int i;
  double start, end, cpu_time_used;

  SparseMatrix A;
  SparseVector arr[n];

  init_sparseMatrix(&A, dim, dim);
  for (i=0; i < n; i++){
    random_init_sparseVector(&arr[i], dim, 0.1*dim);
    append_to_sparseMatrix(&A, &arr[i]);
  }
  double elapsed_time;
  int info;
  //double* AtA = getCovariance_sparseMatrix(&A);
  //print_two_dim("AtA", AtA, dim, dim);

  
  // lapacke_svd 
  double* Adense = (double*) malloc(sizeof(double) * A.nextRow * A.dimension);
  densify_sparseMatrix(&A, Adense);

  double* S = (double*) malloc(sizeof(double) * A.nextRow);
  double* U = (double*) malloc(sizeof(double) * A.nextRow * A.nextRow);
  double* Vt = (double*) malloc(sizeof(double) * A.dimension * A.dimension);
  
  // row wise
  start = clock();
  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', A.nextRow, A.dimension, Adense, A.dimension, S, U, A.nextRow, Vt, A.dimension);
  end = clock();  
  elapsed_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("elapsed time of ROW-WISE lapacke_svd is %f\n",elapsed_time);
  free(S);free(U);free(Vt);
  /*

  //block power
  double* G = (double*) malloc(sizeof(double) * dim * ell);
  for(i=0; i < dim*ell; i++)
    G[i] = ((float)rand()/(float)(RAND_MAX));

  double* p = (double*) malloc(sizeof(double) * ell * A.nextRow);
  double* temp_vec = (double*) malloc(sizeof(double) * ell);
  double* temp_mat = (double*) malloc(sizeof(double) * ell * A.dimension);

  start = clock();
  blockPowerMethod(&A, ell, 0.25, G, p, temp_vec, temp_mat);
  end = clock();

  elapsed_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("elapsed time of block power is %f\n",elapsed_time);
  */
}

int main(){
  main2();
  //blockPower_vs_lapacke();
}
