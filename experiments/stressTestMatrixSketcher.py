import sys
import numpy
import time
import json
from dataHandler import DataHandler
from matrixSketcher import MatrixSketcher

NS = range(10000,100001,10000)
dimensions = range(1000,10001,1000)

for N in NS:
    for  dimension in dimensions:             
        ell = min([100,dimension/2,N/2])
        
        data_handler = DataHandler()
        data_handler.initBeforeMake(dimension,signal_dimension=10,signal_to_noise_ratio=5,\
                                        signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
        
        sketcher = MatrixSketcher(dimension,ell)
        
        matrix_fro_square = 0.0
        t1 = time.time()
        for i in xrange(N):
            v = data_handler.makeRow()
            matrix_fro_square += numpy.sum(v**2)
            sketcher.add(v)
        t2 = time.time()
        time_elapsed = t2-t1
        
        sketch = sketcher.get()
        sketch_fro_square = numpy.linalg.norm(sketch,'fro')**2
        accuracy_bound =  2*(matrix_fro_square - sketch_fro_square)/float(ell) 
        sketch_time = float(time_elapsed)
        result = [('n',N),('d',dimension),('ell',ell),('matrix_fro_square',matrix_fro_square),('sketch_fro_square',sketch_fro_square),('sketch_time',sketch_time),('accuracy_bound',accuracy_bound)]
                        
        sys.stdout.write('%s\n'%json.dumps(result))
        