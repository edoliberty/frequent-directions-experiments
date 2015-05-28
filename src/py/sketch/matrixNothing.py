import numpy
import sys

class MatrixNothing:
    def __init__(self,d,ell,rank):
        self.class_name = 'MatrixNothing'
        self.d = d
        self.ell = ell
        self.input_placeholder = numpy.zeros(self.d)
        self.sketch = numpy.zeros((self.ell,self.d))
        self.i = 0
        self.rank = rank
        
    def __put_in_input_placeholder__(self,v):
        self.input_placeholder.fill(0)
        self.input_placeholder[:self.d] = v.flatten()[:self.d]
                    
    def add(self,v):
        self.__put_in_input_placeholder__(v)
        self.i+=1 
         
    def get(self):
        if self.rank < self.ell:
            return self.sketch
        else:
            return self.sketch[:self.rank,:]
              
if __name__ == '__main__':    
    d = 10
    ell = 5
    matrix = numpy.diag([.9**i for i in xrange(d)])
    mn = MatrixNothing(d,ell)
    for v in matrix:
        mn.add(v)
    sketch = mn.get()
    print sketch
