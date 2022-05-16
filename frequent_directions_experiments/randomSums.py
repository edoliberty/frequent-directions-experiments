from __future__ import absolute_import
from numpy.random import randint, choice
from .matrixSketcherBase import MatrixSketcherBase


class RandomSums(MatrixSketcherBase):
    def __init__(self, d, ell):
        MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = "RandomSums"
        self.signs = [1.0, -1.0]

    def append(self, vector):
        row = randint(self.ell)
        sign = choice(self.signs)
        # v = (sign*vector).tolist()
        # self._sketch[row,:] += v[0]
        self._sketch[row, :] += sign * vector
