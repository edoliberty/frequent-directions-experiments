from numpy import zeros
class MatrixSketcherBase:
    def __init__(self,d,ell):
        self.d = d
        self.ell = ell
        self.B = zeros((self.d, self.ell))
        self._machine_precision = 1e-10
        self._next = zeros(self.d)
        
    def __putInNext__(self,v):
        self._next.fill(0)
        self._next[:self.d] = v.flatten()[:self.d]
            
    def sketch(self,vectors):
        for vector in vectors:
            self.add(vector)
            
    def add(self,vector):
        self.__putInNext__(vector)
    
    def get(self):
        return self.B
    
if __name__ == '__main__':
    import sys
    (n,d,ell) = (100,10,3)
    mbs = MatrixSketcherBase(d,ell)
    for i in xrange(n):
        vector = zeros(d)
        mbs.add(vector)
    mbs.get()