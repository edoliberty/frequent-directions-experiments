from numpy import zeros, sqrt, dot, diag
from numpy.linalg import norm, svd
from matrixSketcherBase import MatrixSketcherBase
from scipy.sparse import lil_matrix as sparse_matrix
from scipy.sparse import csc_matrix as compressed_matrix
from scipy.sparse import rand
from sparsesvd import sparsesvd

from scipy.sparse import rand
from time import time as timer
from numpy import log
from common import truncateSVD
from numpy import cov as covariance
from simIter import simIter

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

        if (self.buffer_nextZeroRow >= self.buffer_ell or self.buffer_nnz >= self.buffer_nnz_threshold):
            self.__rotate__()


        self.buffer[self.buffer_nextZeroRow,:] = vector
        self.buffer_nnz += vector.nnz
        self.buffer_nextZeroRow +=1
      

    def __rotate__(self):
        # First shrink the buffer
        [Vt, ssquared] = simIter(self.buffer, self.ell)

        # insert the shrunk part into the sketch
        if len(ssquared) >= self.ell:
            sShrunk = sqrt(ssquared[:self.ell] - ssquared[self.ell-1])
            self._sketch[self.ell:,:] = dot(diag(sShrunk), Vt[:self.ell,:])
        else:
            self._sketch[self.ell : self.ell+len(s),:] = dot(sqrt(diag(s)), Vt[:len(s),:])

        # resetting the buffer matrix
        del self.buffer
        self.buffer = sparse_matrix( (self.buffer_ell, self.d) )
        self.buffer_nnz = 0
        self.buffer_nextZeroRow = 0

        # A dense shrink of the sketch
        [_,s,Vt] = svd(self._sketch, full_matrices = False)
        if len(s) >= self.ell:
            sShrunk = sqrt(s[:self.ell]**2 - s[self.ell-1]**2)
            self._sketch[:self.ell,:] = dot(diag(sShrunk), Vt[:self.ell,:])
            self._sketch[self.ell:,:] = 0
        else:
            self._sketch[:len(s),:] = dot(diag(s), Vt[:len(s),:])
            self._sketch[len(s):,:] = 0


    def get(self):
        self.__rotate__()
        return self._sketch[:self.ell,:]
     

    

if __name__ == '__main__':
    n = 1000
    d = 400
    k = 5
    ells = range(10,80,10)
    density = 0.01
    A = rand(n, d, density, format = 'lil')   
    
    B = A.todense()
    ATA = covariance(B.T)
    squared_frob_A = norm(B,'fro') ** 2
    print 'squared frob A = ', squared_frob_A

    A_rank_k = truncateSVD(B, k)
    opt_rank_k_err = norm(B - A_rank_k, 'fro') ** 2

   
    for ell in ells:
        sketcher = SparseSketcher(d, ell)

        t_start = timer()
        for row in A:
            sketcher.append(row)
        
        t_end = timer()
        totalSketchTime = t_end-t_start
        
        sketch = sketcher.get()
        #### cov-error #######
        diff = ATA - dot(sketch.transpose(),sketch)
        relative_cov_err = float(norm(diff,2)) / float(squared_frob_A)

        #### proj-error ######
        [u,s,vt] = svd(sketch, full_matrices = False)
        vt = vt[:k, :] 
        projection = dot(B, dot(vt.transpose(), vt))
        proj_err = norm(B - projection, 'fro') ** 2 
        relative_proj_err = float(proj_err) / float(opt_rank_k_err)

        print 'ell=',ell, 'time=',totalSketchTime, 'cov-err=',relative_cov_err, 'proj-err=',relative_proj_err
        print '###################################################################################'
        
