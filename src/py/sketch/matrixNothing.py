import numpy
import sys
from matrixSketcherBase import MatrixSketcherBase


class MatrixNothing (matrixSketcherBase):

    def __init__(self,d,ell):
	MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = 'MatrixNothing'

              
    def append(self, vector):
	self.empty_row_index += 1 
