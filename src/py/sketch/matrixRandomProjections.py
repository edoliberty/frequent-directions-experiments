import sys, numpy
from dataHandler import DataHandler
from matrixMultiplier import MatrixMultiplier

class MatrixRandomProjections:
    def __init__(self,d,ell,rank):
        self.class_name = 'MatrixRandomProjections'
        self.d = d
        self.ell = ell
        self.B = numpy.zeros((self.ell,self.d))
        self.input_placeholder = numpy.zeros(self.d)
        self.i = 0
        self.rank = rank

    def __put_in_input_placeholder__(self,v):
        self.input_placeholder.fill(0)
        self.input_placeholder[:self.d] = v.flatten()[:self.d]
            
    def sketch(self,A):
        for a in A:
            self.add(a)
            
    def add(self,v):
        self.__put_in_input_placeholder__(v)
        row = numpy.random.randint(self.ell)
        new_v = self.input_placeholder/numpy.sqrt(self.ell)
        for row in xrange(self.ell):
            if numpy.random.randint(2) == 1:
                self.B[row,:] += new_v
            else:
                self.B[row,:] -= new_v    
        self.i+=1
        
    def get(self):
        if self.rank < self.ell:
            [U,S,V] = numpy.linalg.svd(self.B,full_matrices = False)
            return V[:self.rank,:]
        else:
            return self.B
                           
if __name__ == '__main__':   
    N = 300
    dimension = 1000
    ell = 10
    data_handler = DataHandler()
    data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=100,\
                                signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
    
    sketcher = MatrixRandomProjections(dimension,ell)
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

    
    
    
