import numpy
import sys
import time


A = numpy.random.randn(10000,1000)


t_start = time.time()
(Q,R) = numpy.linalg.qr(A)
[Usmall,S,V] = numpy.linalg.svd(R,full_matrices=False,compute_uv=True) 
U = numpy.dot(Q,Usmall)
t_end = time.time()
print 'QR svd took %f seconds'%(t_end - t_start)


t_start = time.time()
[U,S,V] = numpy.linalg.svd(A,full_matrices=False,compute_uv=True)
t_end = time.time()
print 'Direct svd took %f seconds'%(t_end - t_start)


