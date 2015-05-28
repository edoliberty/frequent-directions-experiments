import numpy, scipy


class BlockPowerMethod:
    
    def __init__(self):
        self.extra = 5
        self.convergenceEps = 0.00001
    
    def svdsInit(self, mat, V):
        (m,n) = mat.shape
        assert V.shape[0] == n
        assert V.shape[1] <= min(m,n)
        
        [V,_] = numpy.linalg.qr(V)
        phi = 0
        newPhi = 2*self.convergenceEps
        while newPhi > phi + self.convergenceEps:
            phi = newPhi 
            US = numpy.dot(mat,V)
            newPhi = numpy.linalg.norm(US) 
        
            phi = newPhi
            V = numpy.dot(numpy.transpose(mat),US)
            V = numpy.linalg.svd(V,False)[0]
        return V    
         
        
    def svds(self, mat, s):
        (m,n) = mat.shape
        init = numpy.random.randn(n,s)
        res = self.svdsInit(mat,init)
        return res
        
if __name__ == '__main__':
    A = numpy.random.randn(500,300)
    bpm = BlockPowerMethod()
    V = bpm.svds(A,20)
    
    Vnew = numpy.dot(numpy.transpose(A),numpy.dot(A,V))
    
    for j in xrange(20):
        z = numpy.linalg.norm(Vnew[:,j])
        Vnew[:,j] = Vnew[:,j]/z
   
    print numpy.linalg.norm(Vnew - V)**2    

    
    
    #print numpy.dot(numpy.transpose(V),V)
    
    
    
    