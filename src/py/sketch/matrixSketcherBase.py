from numpy import zeros
class MatrixSketcherBase:
    def __init__(self,d,ell):
        self.d = d
        self.ell = ell
        self._sketch = zeros((self.ell, self.d))
        self._machinePrecision = 1e-10
        self._next = zeros(self.d)
        
    # Do NOT overwrite "append". 
    # The right function to overwrite is "__update"         
    def append(self,vector):
        self._next.fill(0)
        self._next[:self.d] = vector.flatten()[:self.d]
        self.__update()
    
    # Convenient looping numpy matrices row by row
    def extend(self,vectors):
        for vector in vectors:
            self.append(vector)
    
    # returns the sketch matrix
    def get(self):
        return self._sketch
            
    # This should take the vector in _next and update the sketch B
    # Should be overwritten by inheriting matrixSketcher classes
    def __update(self):
        pass
    
    # Convenience support for the += operator  append  
    def __iadd__(self,vector):
        self.append(vector)
        return self
         
if __name__ == '__main__':
    import sys
    (n,d,ell) = (100,10,3)
    mbs = MatrixSketcherBase(d,ell)
    for i in xrange(n):
        vector = zeros(d)
        mbs.append(vector)
    mbs.get()