import numpy
import sys
from reservoirSampler import ReservoirSampler
from scipy.sparse import dok_matrix
from scipy import float32
from dataHandler import DataHandler


class MatrixEntrySampler:
    def __init__(self,d,ell):
        self.class_name = 'MatrixEntrySampler'
        self.d = d
        self.ell = ell
        self.nnz = d*ell
        self.number_of_columns = 0
        self.sampler = ReservoirSampler(self.nnz)
                    
    def add(self,v):
        self.number_of_columns += 1
        col = self.number_of_columns
        for (row,val) in enumerate(v):
            self.addOneEntry(row,col,val)
        
    def addOneEntry(self,row,col,val):
        weight = abs(val)
        self.sampler.add((row,col,val),abs(val))

    def get(self):
        self.B = dok_matrix((self.d,self.number_of_columns), dtype=float32)
        items = self.sampler.get(with_probabilities=True)
        for item in items:
            ((row,col,val),p) = item
            self.B[row,col] = val/(p*self.nnz)
        return self.B       
                           
if __name__ == '__main__':
    print '---------------------------------------------'
    print "Running correctness test..."   
    N = 30
    dimension = 10
    ell = 20
    data_handler = DataHandler()
    data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=100,\
                                signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
    
    sketcher = MatrixEntrySampler(dimension,3)


    for i in xrange(N):
        v = data_handler.makeRow()
        sketcher.add(v)

    sketch = sketcher.get()    
    
    print sketch
    
    
    
    
