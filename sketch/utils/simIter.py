from numpy import ceil, log, zeros
from numpy.random import randn
from numpy.linalg import qr, svd
from scipy.sparse import rand
from time import time as timer
from scipy.sparse import csc_matrix



def simIter(matrix, ell):
    [m,d] = matrix.shape
    num_of_iter = int(ceil(4 * log(m)))
    init_vectors = randn(m, ell)
    matrix = csc_matrix(matrix)


    for i in xrange(num_of_iter):
        K = (matrix.transpose()).dot(init_vectors)
        init_vectors = matrix.dot(K)

    [Q,_] = qr(K)
    M = matrix.dot(Q)
    [_,S,U] = svd(M, full_matrices = False)
    return (U[:,:ell].transpose()).dot(Q.transpose()), S[:ell]




if __name__ == '__main__':
    N = 1000
    dimension = 400
    density = 0.1
    ell = 20
    
    A = rand(N, dimension, density, format = 'lil')

    st = timer()
    Z,S = simIter(A, ell)
    print timer() - st
    print Z.shape, len(S)










