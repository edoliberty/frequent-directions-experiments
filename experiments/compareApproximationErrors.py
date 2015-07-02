import sys
from itertools import product
from time import time as timer
from numpy.linalg import svd
from numpy.linalg import norm
from numpy import dot
from numpy import zeros

sys.path.append('../sketch') # needed for imports
from utils.syntheticDataMaker import SyntheticDataMaker
from utils.common import truncateSVD
import bruteForce, frequentDirections, rowSampler, randomProjections, randomSums


if __name__ == '__main__':    
    scketcherClasses = [bruteForce.BruteForce, \
                        rowSampler.RowSampler,\
                        randomProjections.RandomProjections,\
                        randomSums.RandomSums,\
                        frequentDirections.FrequentDirections]
    ns = [500]
    ds = [100]
    ells = range(10,40,10)
    ks = [5]
    rounds = 1
        
    for (n, d, k) in product (ns, ds, ks):
        data_maker = SyntheticDataMaker()
        data_maker.initBeforeMake(d, k, signal_to_noise_ratio=10.0)

        # needed for error computation only
        A = zeros((n, d))
        ATA = zeros((d, d))

        for i in xrange(n):
            vector = data_maker.makeRow()
            A[i,:] = vector 

        squared_frob_A = norm(A,'fro') ** 2
        A_rank_k = truncateSVD(A, k)

        for (sketcherClass, ell, r) in product(scketcherClasses, ells, range(rounds)):
            if ell > d/2:
                continue
                    
            sketcher = sketcherClass(d, ell)
            for i in xrange(n):
                sketcher.append(A[i,:])

            sketch = sketcher.get()    

            if sketcher.class_name == 'BruteForce':
                ATA = sketcher.get()

            diff = ATA - dot(sketch.transpose(),sketch)
            relative_cov_err = norm(diff,2) / squared_frob_A
        
            [u,s,vt] = svd(sketch, full_matrices = False)
            vt = vt[:k, :] 
            projection = dot(A, dot(vt.transpose(), vt))
            proj_err = norm(A - projection, 'fro') ** 2 
            opt_rank_k_err = norm(A - A_rank_k, 'fro') ** 2
            relative_proj_err = float(proj_err) / float(opt_rank_k_err)

            print sketcher.class_name, relative_cov_err, relative_proj_err
             

