from numpy import zeros, sqrt, dot, diag
from numpy.linalg import norm, svd
from matrixSketcherBase import MatrixSketcherBase
from scipy.sparse import lil_matrix as sparse_matrix
from scipy.sparse import csc_matrix as compressed_matrix
from scipy.sparse import rand
import sparsesvd


class SparseSketcher(MatrixSketcherBase):

    def __init__(self, d, ell):
        self.class_name = 'SparseSketcher'
        self.d = d
        self.ell = ell
        self._sketch = zeros( (2*self.ell, self.d) )
        self.sketch_nextZeroRow = 0 

        self.buffer_ell = self.d
        self.buffer = sparse_matrix( (self.buffer_ell, self.d) )
        self.buffer_nnz = 0
        self.buffer_nextZeroRow = 0
        self.buffer_nnz_threshold = 2 * self.ell * self.d

    
    def append(self, vector): 
        if vector.nnz == 0:
            return

        if (self.buffer_nextZeroRow >= self.buffer_ell) or \
           (self.buffer_nextZeroRow >= self.ell and self.buffer_nnz >= self.buffer_nnz_threshold):
            self.__rotate__()

        self.buffer[self.buffer_nextZeroRow,:] = vector
        self.buffer_nnz += vector.nnz
        self.buffer_nextZeroRow +=1
        
    def __rotate__(self):
        # First rotate the the buffer
        [_,s,Vt] = sparsesvd(compressed_matrix(self.buffer), self.ell)
        sShrunk = sqrt(s[:self.ell]**2 - s[self.ell]**2)
        
        # insert the shrunk part into the sketch
        self._sketch[self.ell:,:] = dot(diag(s), Vt[:self.ell,:])
        
        # resetting the buffer matrix
        del self.buffer
        self.buffer = sparse_matrix( (self.buffer_ell, self.d) )
        self.buffer_nnz = 0
        self.buffer_nextZeroRow = 0

        # A dense shrink of the sketch matrix
        [_,s,Vt] = svd(self._sketch, full_matrices = False)
        sShrunk = sqrt(s[:self.ell]**2 - s[self.ell]**2)
        self._sketch[:self.ell,:] = dot(diag(sShrunk), Vt[:self.ell,:])
        self._sketch[self.ell:,:] = 0

    def get(self):
        return self._sketch[:self.ell,:]
     
if __name__ == '__main__':

    print '---------------------------------------------'
    print "Running correctness test..."   
    N = 800
    dimension = 100
    ell = 10

    A = rand(N,dimension,density = 0.0001,format = 'lil')   
    sparse_sketcher = SparseSketcher(dimension, ell)
       
    for vector in A:
        sparse_sketcher.append(vector)
            

    sparse_sketch = sparse_sketcher.get()
    A = A.todense()
    ATA = A.transpose().dot(A)
    ata_trc = norm(A,'fro') ** 2
    diff = ATA - dot(sparse_sketch.transpose(),sparse_sketch)
    diff_two = norm(diff,2)
    
    print 'eps ||ATA - BTB||_2 = ', diff_two 
    print 'eps tr(ATA)/ell', ata_trc/float(ell)
    print '---------------------------------------------'
        
    
