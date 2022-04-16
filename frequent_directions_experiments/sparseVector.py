from __future__ import absolute_import
import sys
from numpy import array, sum, float32, uint32, zeros
from six.moves import range


class SparseVector:

    def __init__(self, d, kvList):
        self.d = d
        kvList = [kv for kv in kvList if kv[0] >= 0 and kv[0] < self.d]
        kvList.sort()
        self.cols = array([kv[0] for kv in kvList], dtype=uint32)
        self.values = array([kv[1] for kv in kvList], dtype=float32)
        self.shape = (1,self.d)
        self.nnz = len(self.cols)
        
        self._normSquare = sum(self.values ** 2)
    
    def todense(self):
        v = zeros(self.shape)
        for i in range(self.nnz):
            v[0, self.cols[i]] = self.values[i]
        return v

    def getNnz(self):
        return self.nnz

    def getNormSquare(self):
        return self._normSquare
        
    def distSquare(self,other):
        return self._normSquare + other._normSquare - 2*self.dot(other)
        
        
        
if __name__ == "__main__":
    d = 30
    sv1 = SparseVector(d, [(1,3.1),(23,0.1),(13,13)])
    sv2 = SparseVector(d, [(12,3.1),(23,0.1),(43,-0.4)])
  
    mat = SparseMatrix()
