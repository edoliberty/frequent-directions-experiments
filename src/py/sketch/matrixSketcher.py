import sys
import numpy
import scipy
from scipy import sparse
from scipy.sparse import linalg
#import zeros,array,dot,sqrt,diag,trace,apply_along_axis,isfinite,sum
#from numpy.linalg import svd,norm
#from locallinalg import svd,norm
#from scipy.sparse.linalg import eigsh as eigsh
#from time import time
#sys.path.append('/Users/edo/Workspace/cmm/svntrunk/CMM/MLExtraction/ProductExtraction/py-dependencies/numpy/linal')
#from time import time

from dataHandler import DataHandler
from matrixMultiplier import MatrixMultiplier

class MatrixSketcher:
    def __init__(self,d,ell,rank=10,c=0.5):
        print 'in init, dimension = ',d
        self.c = c
        self.class_name = 'MatrixSketcher'
        self.d = d
        self.ell = 2*ell
        self.k = int(0.5*self.ell)
        self.b = self.ell - self.k
        self.rank = rank

        self.B = numpy.zeros((self.ell,self.d))
        self.i = 0 
        self.BBT = numpy.zeros((self.ell,self.ell))
        self.C  = numpy.zeros((self.k,self.d))
        
        
        self.machine_precision = 1e-10
        self.input_placeholder = numpy.zeros(self.d)
        
        #self.A_fro = 0.0
        #self.B_fro = 0.0
        
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
        #self.A_fro += sum(self.input_placeholder**2)
        if row >= self.ell - 1:
            self.__rotate__() 
            #self.__rotate_eigsh__()   

    def get(self):
        if self.rank < self.ell:
            [U,S,V] = numpy.linalg.svd(self.B,full_matrices = False)
            return V[:self.rank,:]
        else:
            return self.B
    
        
    def __svd__(self):
        return numpy.linalg.svd(self.B,full_matrices=False)
    
    def __rotate__(self):
        self.BBT = numpy.dot(self.B,self.B.transpose())
        # should be replaced with eigenvalues
        
        #[self.left_sv,self.sv_square,self.left_sv_transpose] = scipy.sparse.linalg.svds(self.BBT, k = self.k, tol = self.machine_precision)
        [self.left_sv,self.sv_square,self.left_sv_transpose] = numpy.linalg.svd(self.BBT,full_matrices=False)
        
        if self.sv_square[0] < self.machine_precision: # we are done
            print 'not come here'
            return
        
        # the largest index of a significantly positive singular value or k-1
        k_pos = 0
        #print self.sv_square
        for i in xrange(self.k):
            if self.sv_square[i] >= self.machine_precision:
                k_pos = i
        #print 'k_pos = ',k_pos
        self.D = numpy.sqrt(1 - self.sv_square[k_pos]/self.sv_square[:k_pos])
        #print 'shape e D=',self.D.shape
        self.sv_scaling_dot_left_sv_transpose = numpy.dot(numpy.diag(self.D),self.left_sv[:,:k_pos].transpose())        
        self.B[:k_pos,:] = numpy.dot(self.sv_scaling_dot_left_sv_transpose,self.B)
        self.B[k_pos:,:] = 0
        #print 'shape e B = ',self.B.shape
         
    def __assert_not_nan__(self,A):
        if not numpy.isfinite(A).all():
            sys.stderr.write('Matrix contains non finite or non number values\n')
            exit()
        
            
    def __rotate_eigsh__(self):
        self.BBT = numpy.dot(self.B,self.B.transpose())
        
        
        [self.sv_square, self.left_sv] = eigsh(self.BBT, self.k,\
                                              which='LM', v0=None, ncv=None, maxiter=100, tol=1e-10,\
                                              return_eigenvectors=True)#, mode='normal')
         
        if self.sv_square[0] < self.machine_precision: # we are done
            return
        
        self.C[:,:] = numpy.dot(self.left_sv.transpose(),self.B)
    
        delta = min(self.sv_square)
        for i in xrange(self.k):           
            j = self.k - 1 - i
            row_norm = numpy.sqrt(numpy.sum(self.C[j,:]**2))
            new_row_norm = numpy.sqrt(max(row_norm**2 - delta, 0.0))
            self.B[i,:] = self.C[j,:]*new_row_norm/max(row_norm,self.machine_precision)
        self.B[self.k:,:].fill(0)
        
if __name__ == '__main__':
    import time
    correctness_test=True
    stress_test=False
    
    
    if correctness_test:
        print '---------------------------------------------'
        print "Running correctness test..."   
        N = 300
        dimension = 1000
        ell = 50
        data_handler = DataHandler()
        data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=100,\
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
        
        #print '------------------'
        #print 'ata_trc=', ata_trc 
        #print 'ata_two=', ata_two 
        #print 'diff_two=', diff_two 
        #print 'ata_trc=', ata_trc
        bound = ata_trc/float(2*ell)/sketcher.c
        if diff_two <= bound:
            print 'Test Passed, accuracy within theoretical guarantee.'
        else:
            print 'Test FAILED'
        
        print 'eps ||ATA||_2 = ', ata_two
        print 'eps ||ATA - BTB||_2 = ', diff_two 
        print 'eps tr(ATA)/c', bound
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



