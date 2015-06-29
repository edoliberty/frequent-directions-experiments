import sys

from scipy.sparse import dok_matrix
from scipy import dot, diag, sqrt, float32
from scipy.sparse.linalg import svds

from utils.reservoirSampler import ReservoirSampler

class EntrySampler:
    def __init__(self,d,ell):
        self.class_name = 'EntrySampler'
        self.d = d
        self.ell = ell
        self.nnz = d*ell
        self.rows = 0
        self.sampler = ReservoirSampler(self.nnz)
                    
    def append(self,v):
        for (col,val) in enumerate(v):
            self.sampler.add((self.rows,col,val),abs(val))
        self.rows += 1

    def get(self):
        B = dok_matrix((self.rows,self.d), dtype=float32)
        for ((row,col,val),p) in self.sampler.get(with_probabilities=True):
            B[row,col] += val/(p*self.nnz)
        covariance = dot(B.transpose(),B)    
        (_,s,Vt) = svds(covariance, k=self.ell, maxiter=50, return_singular_vectors=True)
        return dot(diag(sqrt(s[:self.ell])), Vt[:self.ell,:])


    
