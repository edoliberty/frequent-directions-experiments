#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "randomSum.h"
#include "randomProjection.h"
#include "sparseSketcher.h"
#include "rowSampler.h"
#include "common.h"
#include "frequentDirections.h"


void test_vs_sparsity(){
  /*
  int n = 10000;
  int dim = 1000;
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
    init_sparseMatrix(&A, dim, var_set[i]*n);

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
  */
}


int main(){
  int d = 3;
  int ell = 2;
  int i;

  double** product = (double**) malloc(sizeof(double*) * d);
  for(i=0; i< d; i++)
    product[i] = (double*) malloc(sizeof(double) * ell);

  for(i=0; i< d; i++)
    memset(product[i], 0 , sizeof(double) * ell);

  //for(j=0; j < d * ell; j++)
  //  (*product)[j] = 0;
  //test_vs_sparsity();
}
