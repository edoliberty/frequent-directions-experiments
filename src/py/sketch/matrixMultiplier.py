import numpy
import sys

# bruteforce
#they should give same result
class MatrixMultiplier:
    def __init__(self,d,ell,rank):
        self.class_name = 'MatrixMultiplier'
        self.d = d
        self.ell = ell
        self.temp = numpy.zeros((self.d,self.d))
        self.ATA  = numpy.zeros((self.d,self.d))
        self.input_placeholder = numpy.zeros(self.d)
        self.i = 0
        self.rank = rank
        
    def __put_in_input_placeholder__(self,v):
        self.input_placeholder.fill(0)
        self.input_placeholder[:self.d] = v.flatten()[:self.d]
                    
    def add(self,v):
        self.__put_in_input_placeholder__(v)
        row = self.i % self.d                
        self.temp[row,:] = self.input_placeholder        
        self.i+=1 
         
        if row >= self.d - 1:
            self.__flush__()

    def __flush__(self):
        self.ATA = self.ATA + numpy.dot(self.temp.transpose(),self.temp)
        self.temp.fill(0)
            
    def get(self):
        self.__flush__()
        [U,S,V] = numpy.linalg.svd(self.ATA)
        if self.rank < self.ell:
            return V[:self.rank,:]
        else:
            Aell = numpy.dot(numpy.diag(numpy.sqrt(S[:self.ell])),V[:self.ell,:])
            return Aell

        #Arank = numpy.dot(numpy.diag(numpy.sqrt(S[:self.rank])),V[:self.rank,:])
        #return Arank
             
    def getATA(self):
        self.__flush__()
        return self.ATA
             
if __name__ == '__main__':    
    ell =7
    d = 10
    n = 103
    test_matrix = numpy.random.normal(numpy.zeros((n,d)))
    real_ATA = numpy.dot(test_matrix.transpose(),test_matrix)
    
    mm = MatrixMultiplier(d,ell)
    for v in test_matrix:
        mm.add(v)
    
    mm_ATA = mm.getATA()
    
    print 'TEST: Matrix multiplication'
    if numpy.sum(numpy.abs(real_ATA - mm_ATA)**2) >= 1e-10:
        print "FAIL: matrices are not the same"
    else:
        print "PASS: matrices are the same"
        
        
    print "TEST: svd's are the same as brut force"
    [U,S,V] = numpy.linalg.svd(real_ATA,full_matrices=False)
    real_ell_projection = numpy.dot(V[:ell,:].transpose(),V[:ell,:])
    #why compute like this?
    sketch_SV_ell = mm.get()
    sketch_projected = numpy.dot(sketch_SV_ell, real_ell_projection)
    
    if numpy.sum(numpy.abs(sketch_SV_ell - sketch_projected)) > 1e-10:
        print "FAIL: the projection of the sketch is different from the sketch"
    else:
        print "PASS: the projection of the sketch is like the sketch itself"
      
    
