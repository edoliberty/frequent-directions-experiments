#ifndef ITEMSAMPLER_H
#define ITEMSAMPLER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "sparseVector.h"
#include "common.h"


typedef struct {
  SparseVector* item;
  double item_weight;
  double item_probability;
  double sum_w;
  double machine_precision;

} SingleItemSampler;


void init_itemSampler(SingleItemSampler* self);
void add_itemSampler(SingleItemSampler* self, SparseVector* sv);

#endif
