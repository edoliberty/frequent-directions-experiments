from __future__ import absolute_import
from __future__ import print_function
import sys
from numpy import ceil, log, array, sum, float32, uint32, zeros, empty, arange, concatenate, sqrt,diag
from .sparseVector import SparseVector
from scipy.sparse import rand
from numpy.random import randn
from .utils.common import *
from time import time as timer
from numpy.linalg import qr

class SparseMatrix:
    
    def __init__(self, dim):
        self.rows = empty((1,dim)) 
        self.cols = empty((1,dim))
        self.values = empty((1,dim))
        self.nnz = 0 # number of non-zeros
        self.nextRow = 0
        self.pointer = 0
        self.dimension = dim

    def append(self, vector):
        if(vector.d != self.dimension):
            print("dimension mismatch: can not append this vector to the matrix")
            return

        # extend arrays
        if(self.pointer + vector.nnz > self.rows.shape[1] ):
            z = empty( (1, self.rows.shape[1]) )
            self.rows = concatenate((self.rows,z), axis=1)
            self.cols = concatenate((self.cols,z), axis=1)
            self.values = concatenate((self.values,z), axis=1)
            

        for i in range(vector.nnz):
            self.rows[0, self.pointer] = self.nextRow
            self.cols[0, self.pointer] = vector.cols[i]
            self.values[0, self.pointer] = vector.values[i]
            self.pointer += 1

        self.nextRow += 1
        self.nnz += vector.nnz
        

    def getShape(self):
        return  self.nextRow , self.dimension


    def toDense(self):
        denseMat = zeros(( self.nextRow , self.dimension ))
        rowIndex = self.rows[0,0]
        headPtr = 0

        for ptr in range(self.pointer):
            if ptr != self.pointer-1 and self.rows[0, ptr] != rowIndex :
                for j in range(headPtr, ptr):
                    denseMat[ rowIndex, self.cols[0,j] ] = self.values[0,j]
                # resetting
                rowIndex = self.rows[0,ptr]
                headPtr = ptr

            elif ptr == self.pointer-1:
                if self.rows[0, ptr] == rowIndex :
                    for j in range(headPtr, ptr+1):
                        denseMat[ rowIndex, self.cols[0,j] ] = self.values[0,j]
                elif self.rows[0, ptr] != rowIndex:
                    for j in range(headPtr, ptr):
                        denseMat[ rowIndex, self.cols[0,j] ] = self.values[0,j]
                    rowIndex = self.rows[0,ptr]
                    denseMat[ rowIndex, self.cols[0,ptr] ] = self.values[0,ptr]

        return denseMat



    def sparseShrink(self, ell):
        Z = self.blockpower(ell, 0.25)
        ZtA = self.transposeRightMult(Z)
        [u,s,vt] = svd(ZtA, full_matrices = False)
        for i in range(len(s)):
            s[i] = sqrt(s[i]**2 - s[-1]**2)
        return diag(s).dot(vt)


    def blockpower(self, ell, eps=1):
        n , d = self.getShape()
        init_mat = randn(d, ell)
        num_of_iter = int(10 * ceil(log(d/eps)/eps)) 

        for i in range (num_of_iter):
            [init_mat,_] = qr(init_mat)
            init_mat = self.covarianceMult(init_mat)
    
        K = self.leftMult(init_mat)
        [Q,_] = qr(K)
        del K
        del init_mat
        return Q


    
    ## A^TA * denseMat
    def covarianceMult(self, denseMat):
        rowIndex = self.rows[0,0]
        headPtr = 0
        ptr = 0
        d, ell = denseMat.shape
        temp = zeros((1, ell))
        product = zeros((d, ell))

        while ptr != self.pointer:
            headPtr = ptr
            rowIndex = self.rows[0,headPtr]
            del temp
            temp = zeros((1, ell))

            while ptr != self.pointer and self.rows[0, ptr] == rowIndex:
                temp += denseMat[ self.cols[0,ptr], : ] * self.values[0,ptr]
                ptr += 1

            for j in range(headPtr, ptr):
                product[ self.cols[0,j], : ] += self.values[0,j] * temp[0,:]

        return product


    ## computes G^tA
    def transposeRightMult(self, denseMat):
        ptr = 0
        rowIndex = self.rows[0,ptr]
        m, ell = denseMat.shape
        product = zeros((ell, self.dimension))

        while (ptr != self.pointer):
            rowIndex = self.rows[0,ptr]
            while (ptr != self.pointer and rowIndex == self.rows[0,ptr]):
                for t in range(ell):
                    product[t,self.cols[0,ptr]] += self.values[0,ptr] * denseMat[rowIndex,t]
                ptr += 1

        return product


    ## computes A*G
    def leftMult(self, denseMat):
        rowIndex = self.rows[0,0]
        headPtr = 0
        d, ell = denseMat.shape
        product = zeros((self.nextRow, ell))
        

        for ptr in range(self.pointer):
            # case 1
            if self.rows[0, ptr] != rowIndex and ptr != self.pointer-1: # headPtr -> ptr-1 is one row
                for j in range(headPtr, ptr):
                    product[rowIndex,:] += denseMat[ self.cols[0,j], : ] * self.values[0,j]
                # resetting
                rowIndex = self.rows[0,ptr]
                headPtr = ptr

                # case 2 and 3
            elif ptr == self.pointer-1:
                # case 2
                if self.rows[0, ptr] == rowIndex :
                    for j in range(headPtr, ptr+1):
                        product[rowIndex,:] += denseMat[ self.cols[0,j], : ] * self.values[0,j]
                        
                    # case 3
                elif self.rows[0, ptr] != rowIndex:
                    for j in range(headPtr, ptr):
                        product[rowIndex,:] += denseMat[ self.cols[0,j], : ] * self.values[0,j]

                    # resetting
                    rowIndex = self.rows[0,ptr]
                    product[rowIndex,:] = denseMat[ self.cols[0,ptr], : ] * self.values[0,ptr]

        return product

            

if __name__ == '__main__':
    N = 10000
    dimension = 1000
    density = 0.1
    ell = 10

    A = rand(N, dimension, density, format = 'coo')
    svList, flag = cooToSparseVectorsList(A)    
    sparseMat = SparseMatrix(dimension)
    for sv in svList:
        sparseMat.append(sv)
      
    s = timer()  
    B = sparseMat.sparseShrink(ell)
    e = timer()
    print('elapsed time in python is ',e-s)

    print(B)
    
