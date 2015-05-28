import numpy
import sys

from time import time


class MatrixSketcherSVD:
    def __init__(self,d,ell):
        self.class_name = 'MatrixSketcherSVD'
        self.d = d
        self.ell = ell
        self.k = self.ell/2
        
        #self.k = int(max(numpy.ceil(1/eps),1.0))
        #self.ell = ell
        #if self.ell == None:
        #    self.ell = int(numpy.ceil(1/eps))
        #self.ell = min(max(self.ell,self.k),self.d)        
        
        
        self.b = self.ell - self.k
        
        self.B = numpy.zeros((self.ell,self.d))
        self.i = 0 
        #self.BBT = numpy.zeros((self.ell,self.ell))
        #self.machine_precision = 1e-10
        self.input_placeholder = numpy.zeros(self.d)
        
    def __put_in_input_placeholder__(self,v):
        self.input_placeholder.fill(0)
        self.input_placeholder[:self.d] = v.flatten()[:self.d]
            
    def sketch(self,A):
        for a in A:
            self.add(a)
            
    def add(self,v):
        self.__put_in_input_placeholder__(v)
        row = min((self.i - (self.k - 1)) % (self.b+1) + (self.k-1),self.i)        
        self.B[row,:] = self.input_placeholder
        
        self.i+=1 
         
        if row >= self.ell - 1:
            [left_sv,sv,right_sv] = numpy.linalg.svd(self.B,full_matrices=False,compute_uv=True)
            delta = sv[self.k-1]
            sv = numpy.sqrt(numpy.maximum(sv**2-delta**2,0))
            self.B = numpy.dot(numpy.diag(sv),right_sv)

    def get(self):
        return self.B
         
                           
if __name__ == '__main__':   
    N = 30
    dimension=1000
    ell = 30
    data_handler = DataHandler()
    data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=10000,\
                                signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
    
    
    sketcher = MatrixSketcher(dimension,ell)
    mm = MatrixMultiplier(dimension,ell)
    
    for i in xrange(N):
        v = data_handler.makeRow()
        sketcher.add(v)
        mm.add(v)
        
    ATA = mm.getATA()
    ata_trc = numpy.trace(ATA)
    ata_two = numpy.linalg.norm(ATA,2)
    
    sketch = sketcher.get()
    diff = ATA - numpy.dot(sketch.transpose(),sketch)
    diff_trace = numpy.trace(diff) 
    diff_two = numpy.linalg.norm(diff,2)
    
    
    print 'ata_trc=', ata_trc 
    print 'ata_two=', ata_two 
    print 'diff_two=', diff_two 
    print 'ata_trc=', ata_trc/float(ell)     
        
    