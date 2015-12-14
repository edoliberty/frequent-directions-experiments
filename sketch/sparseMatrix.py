import sys
from numpy import ceil, log, array, sum, float32, uint32, zeros, empty, arange
from sparseVector import SparseVector
from scipy.sparse import rand
from numpy.random import randn
from common import *
from time import time as timer

class SparseMatrix:
    
    def __init__(self, nnz, dim):
        self.rows = empty((1,nnz)) #opt
        self.cols = empty((1,nnz))
        self.values = empty((1,nnz))

        self.nnz = 0 # number of non-zeros
        self.nextRow = 0
        self.pointer = 0

        self.columnDim = dim

    def append(self, vector):
        for i in range(vector.nnz):
            self.rows[0, self.pointer] = self.nextRow
            self.cols[0, self.pointer] = vector.cols[i]
            self.values[0, self.pointer] = vector.values[i]
            self.pointer += 1

        self.nextRow += 1
        self.nnz += vector.nnz
        
    def getShape(self):
        return  self.nextRow , self.columnDim

    def toDense(self):
        denseMat = zeros(( self.nextRow , self.columnDim ))
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


        ## multiplication by a dense matrix 
        def mult(self, denseMat):
            rowIndex = self.rows[0,0]
            headPtr = 0
            d, ell = denseMat.shape
            product = zeros((self.nextRow, ell))
            product_itr = 0

            for ptr in range(self.pointer):
                # case 1
                if self.rows[0, ptr] != rowIndex and ptr != self.pointer-1: # headPtr -> ptr-1 is one row
                    for j in range(headPtr, ptr):
                        product[product_itr,:] += denseMat[ self.cols[0,j], : ] * self.values[0,j]
                    product_itr += 1
    
                    # resetting
                    rowIndex = self.rows[0,ptr]
                    headPtr = ptr


                # case 2 and 3
                elif ptr == self.pointer-1:
                    # case 2
                    if self.rows[0, ptr] == rowIndex :
                        for j in range(headPtr, ptr+1):
                            product[product_itr,:] += denseMat[ self.cols[0,j], : ] * self.values[0,j]
                        product_itr += 1

                    # case 3
                    elif self.rows[0, ptr] != rowIndex:
                        for j in range(headPtr, ptr):
                            product[product_itr,:] += denseMat[ self.cols[0,j], : ] * self.values[0,j]
                        product_itr += 1

                        # resetting
                        rowIndex = self.rows[0,ptr]
                        product[product_itr,:] = denseMat[ self.cols[0,ptr], : ] * self.values[0,ptr]

            return product



if __name__ == '__main__':
    N = 10000
    dimension = 1000
    density = 0.01
    ell = 300

    A = rand(N, dimension, density, format = 'coo')
    svList, flag = cooToSparseVectorsList(A)
    sparseMat = SparseMatrix(N * dimension * density, dimension)
    for sv in svList:
        sparseMat.append(sv)
        
    denseMat = randn(dimension, ell)    
    n , d = sparseMat.getShape()
    eps = 1
    num_of_iter = int(10 * ceil(log(n/eps)/eps)) 
        

    s = timer()
    for i in xrange (1):
        sparseMat.covarianceMult(denseMat)
    e = timer()
    print 'elapsed time is ',e-s

    
