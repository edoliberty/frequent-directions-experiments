import sys, numpy
from matrixSketcherBase import MatrixSketcherBase

class MatrixRandomSums(MatrixSketcherBase):

    def __init__(self , d , ell):
	MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = 'MatrixRandomSums'


    def append(self,vector):
        row = numpy.random.randint(self.ell)
        if numpy.random.randint(2) == 1:
            self._sketch[row,:] += vector.flatten()[:self.d]
        else:
            self._sketch[row,:] -= vector.flatten()[:self.d] 

        
                          
