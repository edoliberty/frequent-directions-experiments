import sys
import numpy

import scipy
import scipy.linalg
import scipy.sparse.linalg
from sparsesvd import sparsesvd
#from scipy import sparse
#from scipy.sparse import linalg
#from sparsesvd import sparsesvd
#from dataHandler import DataHandler
#from matrixMultiplier import MatrixMultiplier
#from matrixSketcher import MatrixSketcher

from scipy.sparse import lil_matrix as sparse_matrix

class SparseMatrixSketcher2:
    def __init__(self,d,ell,target_rank):
        self.class_name = 'SparseMatrixSketcher2'
        self.d = d
        self.ell = ell #it is used to subtract off ell-th sv
        #self.n_sketch_rows = 2 * self.ell
        self.sketch = numpy.zeros((2*self.ell,self.d))
        
        self.buffer_ell = self.d
        self.buffer = sparse_matrix((self.buffer_ell,self.d))
        self.buffer_nnz = 0
        self.buffer_i = 0
        self.buffer_nnz_threshold = self.ell * self.d

        self.total_squared_fro_norm = 0.0
        self.relative_svd_precision = 1e-3
        self.first_empty_row_idx = 0 #this is the index of first empty row in self.sketch
    
    def add(self,v): 
        self.buffer[self.buffer_i,:] = v 
        self.buffer_nnz += v.nnz
        self.buffer_i +=1
        self.total_squared_fro_norm += numpy.linalg.norm(v.data[0]) ** 2 
        
        
        if (self.buffer_i >= self.buffer_ell) or (self.buffer_i >= self.ell and self.buffer_nnz >= self.buffer_nnz_threshold):
            self._shrink_()
            

    '''
    -------------------------------- Shrink --------------------------------------------------------------------------------
    '''
    def _shrink_(self):
        # first shrink the the buffer using a sparse SVD implementation
        #precision_tolerance = numpy.sqrt(self.total_squared_fro_norm)*self.relative_svd_precision/self.ell
       
        rank = numpy.linalg.matrix_rank(self.buffer.todense())
        if rank == 0 : #self.buffer is full zero matrix
            return
        
        rank = min(rank , self.ell)
        [_,s,u] = sparsesvd(scipy.sparse.csc_matrix(self.buffer),rank)

        delta = s.min()**2
        # more efficient than multiplying by a diagonal
        len_s = len(s)
       
        for i in xrange(len_s):
            self.sketch[i + self.first_empty_row_idx] = u[i]*numpy.sqrt(s[i]**2 - delta)

        del self.buffer
        # resetting the buffer matrix
        self.buffer = sparse_matrix((self.buffer_ell,self.d))
        self.buffer_nnz = 0
        self.buffer_i = 0
        

        # a second dense shrink of the sketch matrix
        rank = numpy.linalg.matrix_rank(self.sketch)
        rank = min(rank , self.ell)
        [_,s,u] = sparsesvd(scipy.sparse.csc_matrix(self.sketch),rank)
      
        delta = s.min()**2
        len_s = len(s)
        for i in xrange(len_s):
            self.sketch[i] = u[i]*numpy.sqrt(s[i]**2 - delta)

        self.first_empty_row_idx = len_s - 1


    def get(self):
        self._shrink_()
        return self.sketch[:self.ell].copy()###why copy() here?
    
 
if __name__ == '__main__':
    import time
    import pickle
    correctness_test=True
    stress_test=False
    
    if correctness_test:
        print '---------------------------------------------'
        print "Running correctness test..."   
        N = 800
        dimension = 100
        ell = 10

        A = scipy.sparse.rand(N,dimension,density = 0.1,format = 'lil')        
        sparse_sketcher = SparseMatrixSketcher2(dimension,ell,ell)
       
        for i in xrange(N):
            v = A.getrow(i)
            t = sparse_sketcher.add(v)
            

        sparse_sketch = sparse_sketcher.get()
        print 'shape(sparse_sketch) = ',sparse_sketch.shape
        A = A.todense()
        ATA = A.transpose().dot(A)
        print 'shape(ATA) = ',ATA.shape
        ata_trc = numpy.linalg.norm(A,'fro') ** 2
        diff = ATA - numpy.dot(sparse_sketch.transpose(),sparse_sketch)
        diff_trace = numpy.trace(diff) 
        diff_two = numpy.linalg.norm(diff,2)
        
        #print '------------------'
        #print 'ata_trc=', ata_trc 
        #print 'ata_two=', ata_two 
        #print 'diff_two=', diff_two 
        #print 'ata_trc=', ata_trc
        bound = ata_trc/ell
        if diff_two <= bound:
            print 'Test Passed, accuracy within theoretical guarantee.'
        else:
            print 'Test FAILED'
        
        #print 'eps ||ATA||_2 = ', ata_two
        print 'eps ||ATA - BTB||_2 = ', diff_two 
        print 'eps tr(ATA)/ell', bound , ata_trc
        #print 'FD',numpy.linalg.norm(ATA - numpy.dot(sketch.transpose(),sketch),2)
        print '---------------------------------------------'
        
    
    
    if stress_test:
        print '---------------------------------------------'
        print "Running stress test..."
        t1 = time.time()
        N = 1000000
        dimension = 10000
        ell = 50
        data_handler = DataHandler()
        data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=5,\
                                        signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
        
        sketcher = MatrixSketcher(dimension,ell)


        for i in xrange(N):
            v = data_handler.makeRow()
            sketcher.add(v)
        t2 = time.time()
        time_elapsed = t2-t1
        print "stress_test ended successfully"
        print "ell = %d"%ell
        print "d = %d"%dimension
        print "n = %d"%N
        print "matrix size = %f Gb"%(float(N*dimension*4)/(2**30))
        print "sketch size = %f Gb"%(float(N*ell)*4/(2**30))
        print "time = %f minutes"%(time_elapsed/60.0)
        
        #print "Sketched full %dx%d matrix in %f seconds."%(N,dimension,time_elapsed)
        print '---------------------------------------------'   



