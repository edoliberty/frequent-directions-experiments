import sys
import numpy
import scipy
from scipy import sparse
from scipy.sparse import linalg
from sparsesvd import sparsesvd
from dataHandler import DataHandler
from matrixMultiplier import MatrixMultiplier
from matrixSketcher import MatrixSketcher
from scipy.sparse import lil_matrix as sparse_matrix


class SparseMatrixSketcher:
    def __init__(self,d,ell,target_rank):
        self.class_name = 'SparseMatrixSketcher'
        self.d = d
        self.ell = ell #it is used to subtract off ell-th sv
        self.nnz_threshold = 2 * self.d * self.ell
        self.current_non_zero = 0
        self.i = 0 
        self.B = sparse_matrix((self.ell + self.d , self.d)) 
        self.Brank = 0
        self.target_rank = target_rank
        self.machine_precision = 1e-10
        self.total_squared_fro_norm = 0.0
        
            
    def add(self,v):
        self.B[self.i,:] = v
        self.current_non_zero += v.nnz     
        self.i+=1 
        self.total_squared_fro_norm += numpy.linalg.norm(v.data[0]) ** 2 

        if (self.current_non_zero >= self.nnz_threshold) or (self.i >= self.ell + self.d):
            self.__rotate__() 


    def get(self):
        self.__rotate__()
        return self.B[0:self.target_rank,:].copy()
            
    '''
    def get(self):
        if self.target_rank < self.ell:
            [U,S,V] = numpy.linalg.svd(self.B,full_matrices = False)
            return V[:self.target_rank,:]
        else:
            return self.B
    '''

    def __svd__(self):
        V = numpy.random.randn(self.d, 2*self.ell)
        for i in range(int(round(numpy.log(self.d)))):
            #print 'i = ',i
            X = self.B.dot(V)
            #print 'B.shape',self.B.shape
            Y = self.B.transpose().dot(X)
            k = numpy.linalg.matrix_rank(Y)-1
            if k == 0:
                [V,S,U] = numpy.linalg.svd(Y, full_matrices = False)
            else:
                [V,S,U] = scipy.sparse.linalg.svds(Y, k = k, tol = self.machine_precision)               
        return S,V
    

    def __rotate__(self):
        #print 'in rotate'
        [_,S,Vt] = sparsesvd(scipy.sparse.csc_matrix(self.B), self.ell)
       
        delta = S.min() ** 2
        shrinked_s = numpy.sqrt(S**2 - delta)
        temp = numpy.dot(numpy.diag(shrinked_s),Vt)
        self.B = sparse_matrix((self.ell + self.d , self.d))
        self.B[:temp.shape[0],:] = sparse_matrix(temp)
        self.current_non_zero = self.B.nnz
        self.i = len(S) - 1
       
        
if __name__ == '__main__':
    import time
    import pickle
    correctness_test=True
    stress_test=False
    
    if correctness_test:
        print '---------------------------------------------'
        print "Running correctness test..."   
        N = 10000
        dimension = 1000
        ell = 20

        A = scipy.sparse.rand(N,dimension,density = 0.05,format = 'lil')        
        pickle.dump(A, open("dataset/sparse_data_N=10000_d=1000_dens=0.05", "wb"))
        #A = pickle.load(open('sparse_data_N=10000_d=1000','rb'))
        print 'after dumping dataset'
        sparse_sketcher = SparseMatrixSketcher(dimension,ell,ell)
        fd = MatrixSketcher(dimension, ell, ell)

        for i in xrange(N):
            v = A.getrow(i)
            sparse_sketcher.add(v)
            fd.add(v.todense())

        sparse_sketch = sparse_sketcher.get()
        ATA = A.transpose().dot(A)
        ata_trc = sparse_sketcher.total_squared_fro_norm
        sts = sparse_sketch.transpose().dot(sparse_sketch)
        diff = ATA - sts
        diff_two = numpy.linalg.norm(diff.todense(),2)
        print 'diff_two = ',diff_two

        fd_sketch = fd.get()
        fd_sketch = scipy.sparse.lil_matrix(fd_sketch)
        fd_sketch = fd_sketch.todense()
        A2 = A.todense()
        ATA2 = numpy.dot(A2.transpose(), A2)
        sparse_sketch2 = scipy.sparse.lil_matrix(sparse_sketch).todense()
        print 'doff_two in dense form = ', numpy.linalg.norm(ATA2 - numpy.dot(sparse_sketch2.transpose(),sparse_sketch2),2)

        print ata_trc, numpy.linalg.norm(A2,'fro') ** 2       
        print 'second diff_two = ', numpy.linalg.norm(ATA2 - numpy.dot(fd_sketch.transpose(),fd_sketch),2)

        '''
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
        '''

        
    
    
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



