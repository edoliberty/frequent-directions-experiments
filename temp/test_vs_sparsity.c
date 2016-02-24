#include "common.h"
#include "frequentDirections.h"
#include "sparseSketcher.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void test_vs_sparsity(){
  int n = 10000;
  int dim = 1000;
  int ell = 50;
  int k = 10;
  int exp_no = 6;
  int var_set[] = {0.01 * dim, 0.05 * dim, 0.1*dim, 0.3*dim, 0.5*dim, 0.7*dim}; 
  //{0.005 * dim, 0.01 * dim, 0.05 * dim, 0.1*dim, 0.3*dim, 0.5*dim};

  double start, end, cpu_time_used;
  SparseMatrix A;
  SparseVector arr[n];
  SparseSketcher sfd;
  FrequentDirections fd;

  double sfd_cov_err[exp_no], sfd_proj_err[exp_no], sfd_time[exp_no];
  double fd_cov_err[exp_no], fd_proj_err[exp_no], fd_time[exp_no];

  double tailSquaredFrob;

  init_sparseSketcher(&sfd, ell, dim);
  init_fd(&fd, ell, dim);
  

  for(int i=0; i<exp_no; i++){
    printf("i = %d\n", i);

    sfd_time[i] = 0;
    fd_time[i] = 0;
    init_sparseMatrix(&A, dim, n);

    for (int j=0; j < n; j++){
      skew_init_sparseVector(&arr[j], dim, var_set[i], (int) (1.5 *var_set[i]), 0.9);
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
    
    tailSquaredFrob = topRank(&A, k);
    
    //SFD
    start = clock();
    get_sparseSketch(&sfd);
    sfd_time[i] += (double) (clock() - start);
    sfd_time[i] = sfd_time[i] / CLOCKS_PER_SEC;
    sfd_cov_err[i] = computeRelCovErr(&A, sfd.sketch, ell, dim);
    sfd_proj_err[i] = computeRelProjErr(&A, sfd.sketch, ell, dim, k, tailSquaredFrob);

    
    //FD
    start = clock();
    get_fdSketch(&fd);
    fd_time[i] += (double) (clock() - start);
    fd_time[i] = fd_time[i] / CLOCKS_PER_SEC;
    fd_cov_err[i] = computeRelCovErr(&A, fd.sketch, ell, dim);
    fd_proj_err[i] = computeRelProjErr(&A, fd.sketch, ell, dim, k, tailSquaredFrob);  
  }  

  printf("SFD:\n");
  print_one_dim_double("\'proj\':", sfd_proj_err, exp_no);
  print_one_dim_double("\'cov\':", sfd_cov_err, exp_no);  
  print_one_dim_double("\'time\':", sfd_time, exp_no);

  print_one_dim_double("\'proj\':", fd_proj_err, exp_no);    
  print_one_dim_double("\'cov\':", fd_cov_err, exp_no);
  print_one_dim_double("\'time\':", fd_time, exp_no);

}
