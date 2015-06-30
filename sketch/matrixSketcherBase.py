from numpy import zeros

class MatrixSketcherBase:
    
    def __init__(self, d, ell):
        self.d = d
        self.ell = ell
        self._sketch = zeros((self.ell, self.d))

    # Appending a row vector to sketch
    def append(self,vector):
        pass
    
    # Convenient looping numpy matrices row by row
    def extend(self,vectors):
        for vector in vectors:
            self.append(vector)
    
    # returns the sketch matrix
    def get(self):
        return self._sketch
    
    # Convenience support for the += operator  append  
    def __iadd__(self,vector):
        self.append(vector)
        return self
