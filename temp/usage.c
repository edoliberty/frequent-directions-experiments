#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "randomSum.h"
#include "randomProjection.h"
#include "sparseSketcher.h"
#include "rowSampler.h"
#include "common.h"
#include "frequentDirections.h"

int main(){
  int n = 1000;
  int dim = 100;
  int ell = 0;
  int ell_set[] = {10,20,30,40};
  int i, ell_counter;
  double start, end, cpu_time_used;

  SparseMatrix A;
  SparseSketcher sfd;
  SparseVector arr[n];
  RandomSum hashing;
  RandomProjection proj;
  RowSampler sampler;
  FrequentDirections fd;

  init_sparseMatrix(&A, dim, dim);
  for (i=0; i < n; i++){
    random_init_sparseVector(&arr[i], dim);
    append_to_sparseMatrix(&A, &arr[i]);
  }

  printf("after initialization of A\n");

  
  // FD
  double fd_cov_err[4];
  double fd_time[4];
  for(ell_counter = 0; ell_counter < 4; ell_counter++){
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
    fd_cov_err[ell_counter] = computeRelCovErr(&A, fd.sketch, ell, dim);
    
  }
  print_one_dim_double("FD COV ERR= ", fd_cov_err, 4);
  print_one_dim_double("FD TIME= ", fd_time, 4);
  

  // SFD
  double sfd_cov_err[4];
  double sfd_time[4];
  for(ell_counter = 0; ell_counter < 4; ell_counter++){
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
    sfd_cov_err[ell_counter] = computeRelCovErr(&A, sfd.sketch, ell, dim);
    
  }
  print_one_dim_double("SFD COV ERR= ", sfd_cov_err, 4);
  print_one_dim_double("SFD TIME= ", sfd_time, 4);


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


