from __future__ import absolute_import
from numpy import zeros, sqrt, dot, diag, ceil, log
from numpy.random import randn
from numpy.linalg import norm, svd, qr, eigh
from scipy.sparse import lil_matrix as sparse_matrix
from scipy.sparse import csc_matrix, rand

from .matrixSketcherBase import MatrixSketcherBase
from six.moves import range


# simultaneous iterations algorithm
# inputs: matrix is input matrix, ell is number of desired right singular vectors
# outputs: transpose of approximated top ell singular vectors, and first ell singular values
def simIter(matrix, ell):
    [m,d] = matrix.shape
    num_of_iter = int(ceil(4 * log(m)))
    init_vectors = randn(m, ell)
    matrix = csc_matrix(matrix)
    matrix_trans = matrix.transpose()

    for i in range(num_of_iter):
        init_vectors = matrix.dot((matrix_trans).dot(init_vectors))

    [Q,_] = qr((matrix_trans).dot(init_vectors))
    M = matrix.dot(Q)

    [_,S,U] = svd(M, full_matrices = False)

    return (U[:,:ell].transpose()).dot(Q.transpose()), S[:ell]

  
# sparse frequent directions sketcher
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
        [Vt, s] = simIter(self.buffer, self.ell)

        # insert the shrunk part into the sketch
        if len(s) >= self.ell:
            sShrunk = sqrt(s[:self.ell]**2 - s[self.ell-1]**2)
            self._sketch[self.ell:,:] = dot(diag(sShrunk), Vt[:self.ell,:])
        else:
            self._sketch[self.ell : self.ell+len(s),:] = dot(diag(s), Vt[:len(s),:])


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
  
