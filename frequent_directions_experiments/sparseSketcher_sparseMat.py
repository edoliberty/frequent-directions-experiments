from __future__ import absolute_import
from __future__ import print_function
from numpy import zeros, sqrt, dot, diag, ceil, log
from numpy.random import randn
from numpy import cov as covariance
from numpy.linalg import norm, svd, qr, eigh
from scipy.sparse import lil_matrix, csc_matrix, csr_matrix, dok_matrix, rand
from time import time as timer
import pickle

from .matrixSketcherBase import MatrixSketcherBase
from .utils.common import truncateSVD
from .blockPower import blockpower
from .sparseVector import SparseVector
from .frequentDirections import FrequentDirections as FD

from .sparseMatrix import SparseMatrix

# sparse frequent directions sketcher
class SparseSketcher(MatrixSketcherBase):

    def __init__(self, d, ell):
        self.class_name = 'SparseSketcher_sparseMatrix'
        self.d = d
        self.ell = ell
        self._sketch = zeros( (2 * self.ell, self.d) )

        self.buffer_nnz_threshold = 2 * self.ell * self.d
        self.buffer = SparseMatrix (self.buffer_nnz_threshold)

 
    def append(self, vector): 
        if (self.buffer.nnz >= self.buffer_nnz_threshold):
            self.__rotate__()

        self.buffer.append(vector)
       
      

    def __rotate__(self):
        # First shrink the buffer
        [s,vt] = blockpower(self.buffer, self.ell)

        # insert the shrunk part into the sketch
        if len(s) < self.ell:
            self._sketch[self.ell : self.ell+len(s),:] = dot(diag(s), vt[:len(s),:])
        else: #len(s) == self.ell
            sShrunk = sqrt(s[:self.ell]**2 - s[self.ell-1]**2)
            self._sketch[self.ell:,:] = dot(diag(sShrunk), vt[:self.ell,:])


        # resetting the buffer matrix
        del self.buffer
        self.buffer = lil_matrix( (self.d, self.d) )
        self.buffer_nnz = 0
        self.buffer_nextRow = 0


        # A dense shrink of the sketch
        [_,s,vt] = svd(self._sketch, full_matrices = False)
        if len(s) >= self.ell:
            sShrunk = sqrt(s[:self.ell]**2 - s[self.ell-1]**2)
            self._sketch[:self.ell,:] = dot(diag(sShrunk), vt[:self.ell,:])
            self._sketch[self.ell:,:] = 0
        else:
            self._sketch[:len(s),:] = dot(diag(s), vt[:len(s),:])
            self._sketch[len(s):,:] = 0


    def get(self):
        self.__rotate__()
        return self._sketch[:self.ell,:]
        
 
        
if __name__ == '__main__':
    # make input data
    n = 1000
    d = 400
    k = 5
    ells = list(range(10,21,10))
    density = 0.1
    A = rand(n, d, density, format = 'lil')   
    
  
    # error computation
    B = A.todense()
    ATA = covariance(B.T)
    squared_frob_A = norm(B,'fro') ** 2
    A_rank_k = truncateSVD(B, k)
    opt_rank_k_err = norm(B - A_rank_k, 'fro') ** 2

   
    for ell in ells:
        sketcher = SparseSketcher(d, ell)

        t_start = timer()
        for sv in A:
            sketcher.append(sv)
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

        print('sparse: ell=',ell, 'time=',totalSketchTime, 'cov-err=',relative_cov_err, 'proj-err=',relative_proj_err)

        ############### FD ###################################3
        sketcher = FD(d, ell)
        t_start = timer()
        for sv in A:
            sketcher.append(sv)
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

        print('DenseFD: ell=',ell, 'time=',totalSketchTime, 'cov-err=',relative_cov_err, 'proj-err=',relative_proj_err)
