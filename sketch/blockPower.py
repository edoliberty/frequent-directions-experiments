#import numpy, scipy
from numpy.random import randn
from numpy import ceil, log, zeros
from numpy.linalg import qr, svd

def blockpower(sparseMat, ell, eps=1):
    n , d = sparseMat.getShape()
    init_mat = randn(d, ell)
    num_of_iter = int(10 * ceil(log(n/eps)/eps)) #constant 10 should be found experimentally based on eps

    for i in xrange (num_of_iter):
        init_mat = sparseMat.covarianceMult(init_mat)
        #K = mat.dot(init_mat)
        #init_mat = (mat.transpose()).dot(K)
    
    K = sparseMat.mult(init_mat)

    [Q,_] = qr(K)
    #M = (Q.transpose()).dot(mat)
    M = (mat.transpose()).dot(Q) #computing transpose of what we need

    [U,S,_] = svd(M, full_matrices = False)
            
    return S, U[:,:ell].transpose() # U is ell*d
    #return (U[:,:ell].transpose()).dot(Q.transpose()), S[:ell] #this step might violate sing val bound we want

        
if __name__ == '__main__':
    A = numpy.random.randn(500,300)
    bpm = BlockPower()
    V = bpm.svds(A,20)
    
    Vnew = numpy.dot(numpy.transpose(A),numpy.dot(A,V))
    
    for j in xrange(20):
        z = numpy.linalg.norm(Vnew[:,j])
        Vnew[:,j] = Vnew[:,j]/z
   
    print numpy.linalg.norm(Vnew - V)**2    

    
    
    #print numpy.dot(numpy.transpose(V),V)
    
    
    
    
