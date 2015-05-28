
from dataHandler import DataHandler
from matrixSketcher import MatrixSketcher
import time
import numpy

print '---------------------------------------------'
print "Initiating data"
t1 = time.time()
N = 21602
dimension = 1017
ell =48
data_handler = DataHandler()
data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=5,\
                                signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')

sketcher = MatrixSketcher(dimension,ell)


matrix_fro_norm_square = 0.0
print "Sketching the matrix"
for i in xrange(N):
    v = data_handler.makeRow()
    matrix_fro_norm_square += sum(v**2)
    sketcher.add(v)


print "Sketching ocomplete."
final_sketch = sketcher.get()


t2 = time.time()
time_elapsed = t2-t1

print "Test ended successfully"
print "ell = %d"%ell
print "d = %d"%dimension
print "n = %d"%N
print "matrix size = %f Gb"%(float(N*dimension*4)/(2**30))
print "sketch size = %f Gb"%(float(N*ell)*4/(2**30))
print "the sketch is a numpy array of shape:", final_sketch.shape
print "sketching took:", time_elapsed, "seconds"

sketch_fro_norm_square = numpy.linalg.norm(final_sketch,'fro')**2
print "tr(AA^T) =", matrix_fro_norm_square
print "||AA^T - BB^T|| <=", (matrix_fro_norm_square - sketch_fro_norm_square)/(ell*sketcher.c)
