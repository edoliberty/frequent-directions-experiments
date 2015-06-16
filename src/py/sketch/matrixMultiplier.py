import numpy
import sys
from matrixSketcherBase import MatrixSketcherBase

# bruteforce
#they should give same result
class MatrixMultiplier(MatrixSketcherBase):

    def __init__(self, d, ell):
	self.d = d
        self.class_name = 'MatrixMultiplier'
        self.temp = numpy.zeros((self.d,self.d))
        self.ATA  = numpy.zeros((self.d,self.d))
	self.i = 0
        
                    
    def append(self,vector):
        row = self.i % self.d                
        self.temp[row,:] = vector.flatten()[:self.d]      
        self.i += 1 

        if row >= self.d - 1:
            self.__flush__()


    def __flush__(self):
        self.ATA = self.ATA + numpy.dot(self.temp.transpose(),self.temp)
        self.temp.fill(0)

            
    def get(self):
        self.__flush__()
	Aell = numpy.dot(numpy.diag(numpy.sqrt(S[:self.ell])),V[:self.ell,:])
        return Aell

             
    def getATA(self):
        self.__flush__()
        return self.ATA

             

