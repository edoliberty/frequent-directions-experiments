#include "singleItemSampler.h"

void init_itemSampler(SingleItemSampler* self){
  self-> item = NULL;
  self-> item_weight = 0;
  self-> item_probability = 0;
  self-> sum_w = 0;
  self-> machine_precision = 1e-10;
  srand(time(NULL));
}


void add_itemSampler(SingleItemSampler* self, SparseVector* sv){
  self-> sum_w += sv-> squaredNorm;
  double p = sv-> squaredNorm / max(self-> sum_w , self-> machine_precision);

  double randomVal = rand()/(RAND_MAX+1.0);

  if (randomVal < p){
    self-> item = sv;
    self-> item_weight = sv-> squaredNorm;
    self-> item_probability = p;
  }else{
    self-> item_probability = self-> item_probability * (1.0-p);
  }
}
