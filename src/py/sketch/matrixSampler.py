import numpy
import sys
from random import random



class singleItemSampler():
    def __init__(self):
        self.item = None
        self.item_weight = 0.0
        self.item_probability = 0.0
        self.sum_w = 0.0
        self.machine_precision = 1e-10
        

    def add(self,item,w=1):
        w = float(w)
        if w <= 0.0:
            return
        self.sum_w += w
        p = w/max(self.sum_w,self.machine_precision)
        if random() < p or item == None:
            self.item = item
            self.item_weight = w
            self.item_probability = p
        else:
            self.item_probability = self.item_probability*(1.0 - p)
            
    def get(self):
        return (self.item ,self.item_weight, item_probability)
        
class MatrixSampler:
    def __init__(self,d,ell,rank):
        self.class_name = 'MatrixSampler'
        self.d = d
        self.ell = ell
        self.samplers = [singleItemSampler() for i in xrange(self.ell)]
        self.B = numpy.zeros((self.ell,self.d))
        self.input_placeholder = numpy.zeros(self.d)
        self.machine_precision = 1e-10
        self.rank = rank

    def __put_in_input_placeholder__(self,v):
        self.input_placeholder.fill(0)
        self.input_placeholder[:self.d] = v.flatten()[:self.d]
                    
    def add(self,v):
        self.__put_in_input_placeholder__(v)
        new_v = self.input_placeholder.copy()
        row_norm_square = numpy.sum(new_v**2)
 
        for i in xrange(self.ell):
            self.samplers[i].add(new_v,row_norm_square)
                    
    def get(self):
        for (i,sampler) in enumerate(self.samplers):
            p = sampler.item_probability
            row = sampler.item
            if p > self.machine_precision and not row == None:
                self.B[i,:] = row/(numpy.sqrt(p*float(self.ell)))

        if self.rank < self.ell:
            [U,S,V] = numpy.linalg.svd(self.B,full_matrices = False)
            return V[:self.rank,:]
        else:
            return self.B
                           
if __name__ == '__main__':
    d = 10
    ell = 5
    matrix = numpy.diag([1 for i in xrange(d)])
    ms = MatrixSampler(d,3)
    for v in matrix:
        ms.add(v)
    sketch = ms.get()
    
    
    print sketch
    
    
    
    
