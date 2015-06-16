import sys, numpy
from dataHandler import DataHandler
from matrixMultiplier import MatrixMultiplier
from matrixSketcherBase import MatrixSketcherBase


class MatrixRandomProjections(MatrixSketcherBase):

    def __init__(self , d , ell):
	MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = 'MatrixRandomProjections'
 
            
    def append(self,vector):
        rescaled_vector = vector.flatten()[:self.d] / numpy.sqrt(self.ell)

        for row in xrange(self.ell):
            if numpy.random.randint(2) == 1:
                self._sketch[row,:] += rescaled_vector
            else:
                self._sketch[row,:] -= rescaled_vector   
       
        
