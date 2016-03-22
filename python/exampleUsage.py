import sys
from numpy.linalg import norm
from numpy import dot

from utils.syntheticDataMaker import SyntheticDataMaker
from frequentDirections import FrequentDirections

n = 500
d = 100
ell = 20
k = 5

# this is only needed for generating input vectors
dataMaker = SyntheticDataMaker()
dataMaker.initBeforeMake(d, k, signal_to_noise_ratio=10.0)                                                                                                                                                                                                                                                                                                                         

# This is where the sketching actually happens
sketcher = FrequentDirections(d,ell)
for i in xrange(n):
    row = dataMaker.makeRow()
    sketcher.append(row)
sketch = sketcher.get()

approxCovarianceMatrix = dot(sketch.transpose(),sketch)

print approxCovarianceMatrix.shape 





