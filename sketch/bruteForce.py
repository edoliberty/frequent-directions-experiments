from numpy import zeros, dot, outer, diag, sqrt
from numpy.linalg import svd
import sys
from matrixSketcherBase import MatrixSketcherBase

class BruteForce(MatrixSketcherBase):
    def __init__(self, d, ell):
        self.d = d
        self.ell = ell
        self.class_name = 'BruteForce'
        self.covariance = zeros( (self.d,self.d) )    
                    
    def append(self,vector):
        self.covariance += outer(vector,vector)
                
    def get(self):
	return self.covariance
        #(U,s,Vt) = svd(self.covariance)
        #return dot(diag(sqrt(s[:self.ell])), Vt[:self.ell,:])


             

