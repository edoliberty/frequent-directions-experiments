from numpy import zeros, max, sqrt, dot, diag
from numpy.linalg import svd
from matrixSketcherBase import MatrixSketcherBase

class FrequentDirections(MatrixSketcherBase):

    def __init__(self , d, ell):
        self.class_name = 'FrequentDirections'
        self.d = d
        self.ell = ell
        self.m = 2*self.ell
        self._sketch = zeros( (self.m, self.d) ) 
        self.nextZeroRow = 0
                 
    def append(self,vector):     
    	if self.nextZeroRow >= self.m:
            self.__rotate__()
        self._sketch[self.nextZeroRow,:] = vector 
        self.nextZeroRow += 1
        
    def __rotate__(self):
        [_,s,Vt] = svd(self._sketch , full_matrices=False)
        sShrunk = sqrt(s[:self.ell]**2 - s[self.ell]**2)
        self._sketch[:self.ell,:] = dot(diag(sShrunk), Vt[:self.ell,:])
        self._sketch[self.ell:,:] = 0
        self.nextZeroRow = self.ell
         
    def get(self):
        return self._sketch[:self.ell,:]
    
