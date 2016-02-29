from numpy.linalg import svd
from numpy import dot
from numpy import diagflat 

def truncateSVD(A, k):
    U, s, Vt = svd(A, full_matrices = False)
    opt = dot(U[:, 0:k ], dot(diagflat(s[0:k]), Vt[0:k, :]))
    return opt
