import sys
import numpy
from matrixSketcherBase import MatrixSketcherBase
from dataHandler import DataHandler
from numpy import zeros

class MatrixSketcher(MatrixSketcherBase):

    def __init__(self , d, ell, c = 0.5):
	MatrixSketcherBase.__init__(self, d, ell)
        self.class_name = 'MatrixSketcher'
        self.shrink_index = self.ell * c

	
            
    def append(self,vector):     
	#MatrixSketcherBase.append(self, vector)
	self._sketch += vector	

        if self.empty_row_index >= self.ell - 1:
            self.__rotate__() 
	    self.empty_row_index = self.shrink_index



    def __rotate__(self):
	self.BBT = numpy.dot(self._sketch , self._sketch.transpose())
        [self.left_sv , self.sv_square , self.left_sv_transpose] = numpy.linalg.svd(self.BBT , full_matrices=False)

        self.D = numpy.sqrt(1 - self.sv_square[self.shrink_index] / self.sv_square[:self.shrink_index])
        self.sv_scaling_dot_left_sv_transpose = numpy.dot(numpy.diag(self.D) , self.left_sv_transpose[:self.shrink_index,:])
        
        self._sketch[:self.shrink_index,:] = numpy.dot(self.sv_scaling_dot_left_sv_transpose , self._sketch)
        self._sketch[self.shrink_index:,:] = 0
    
        
            
        
if __name__ == '__main__':
    import time
    correctness_test=True
    stress_test=False
    
    
    if correctness_test:
        print '---------------------------------------------'
        print "Running correctness test..."   
        N = 300
        dimension = 100
        ell = 50
        data_handler = DataHandler()
        data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=100,\
                                    signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
        
        sketcher = MatrixSketcher(dimension, ell)
        
        for i in xrange(N):
            v = data_handler.makeRow()
            sketcher.append(v)
            #mm.add(v)
            
        #ATA = mm.getATA()
        #ata_trc = numpy.trace(ATA)
        #ata_two = numpy.linalg.norm(ATA,2)
        
        sketch = sketcher.get()
	print sketch.shape
        
    
    
