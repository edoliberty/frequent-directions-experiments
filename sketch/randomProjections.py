import sys
from numpy import outer, sqrt
from numpy.random import randn
from matrixSketcherBase import MatrixSketcherBase

class RandomProjections(MatrixSketcherBase):

    def __init__(self , d , ell):
        MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = 'RandomProjections'
 
    def append(self,vector):        
        randomVector = randn(self.ell)/sqrt(self.ell)
        self._sketch += outer(randomVector,vector) 
       
        
