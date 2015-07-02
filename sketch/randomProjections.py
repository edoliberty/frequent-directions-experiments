import sys
from numpy import outer, sqrt
from numpy.random import choice
from matrixSketcherBase import MatrixSketcherBase

class RandomProjections(MatrixSketcherBase):

    def __init__(self , d , ell):
        MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = 'RandomProjections'
	self.rescaled_signs = [-1.0, 1.0]/sqrt(self.ell)
 
    def append(self,vector):        
        randomVector = choice(self.rescaled_signs, self.ell)
        self._sketch += outer(randomVector,vector) 
       
        
