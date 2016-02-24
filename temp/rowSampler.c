#include "rowSampler.h"

void init_rowSampler(RowSampler* self, int ell, int dim ){
  self-> class_name = "RowSampler";
  self-> dimension = dim;
  self-> ell = ell;
  self-> sketch = (double*) malloc(sizeof(double) * ell * dim);
  self-> samplers = (SingleItemSampler*) malloc(sizeof(SingleItemSampler) * ell);
  memset(self-> sketch, 0 , sizeof(double) * ell * dim);
}

void append_to_rowSampler(RowSampler* self, SparseVector* sv){
  int i;
  for(i=0; i < self-> ell; i++)
    add_itemSampler(&(self-> samplers[i]), sv);
}


void get_rowSamplerSketch(RowSampler* self){

  SparseVector* item;
  double item_prob;

  for(int i=0; i < self-> ell; i++){
    item = (self-> samplers[i]).item;
    item_prob = (self-> samplers[i]).item_probability;
    for(int j=0; j< item-> nnz; j++)
      self-> sketch[i * self-> dimension + item-> cols[j]] = (item-> values[j]) / sqrt(item_prob * self-> ell);
  }
}

