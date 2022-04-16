from __future__ import absolute_import
import numpy
from numpy.linalg import norm
from random import random
from .matrixSketcherBase import MatrixSketcherBase


class RowSampler(MatrixSketcherBase):
    def __init__(self, d, ell):
        MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = "RowSampler"
        self.samplers = [singleItemSampler() for i in range(self.ell)]

    def append(self, vector):
        row_norm_square = norm(vector) ** 2
        for i in range(self.ell):
            self.samplers[i].add(vector, row_norm_square)

    def get(self):
        for (i, sampler) in enumerate(self.samplers):
            p = sampler.item_probability
            row = sampler.item
            if row is not None:
                self._sketch[i, :] = row / (numpy.sqrt(p * float(self.ell)))
        return self._sketch


class singleItemSampler:
    def __init__(self):
        self.item = None
        self.item_weight = 0.0
        self.item_probability = 0.0
        self.sum_w = 0.0
        self.machine_precision = 1e-10

    def add(self, item, w=1):
        w = float(w)
        if w <= 0.0:
            return
        self.sum_w += w
        p = w / max(self.sum_w, self.machine_precision)
        if random() < p or self.item is None:
            self.item = item
            self.item_weight = w
            self.item_probability = p
        else:
            self.item_probability = self.item_probability * (1.0 - p)

    def get(self):
        return (self.item, self.item_weight, self.item_probability)
