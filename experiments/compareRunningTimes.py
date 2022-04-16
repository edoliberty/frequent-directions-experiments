import sys
from itertools import product
from time import time as timer
from numpy import cov as covariance

sys.path.append('../sketch') # needed for imports
from utils.syntheticDataMaker import SyntheticDataMaker
import bruteForce, frequentDirections, rowSampler, randomProjections, randomSums

if __name__ == '__main__':    
    sketcherClasses = [bruteForce.BruteForce, \
                        rowSampler.RowSampler,\
                        randomProjections.RandomProjections,\
                        randomSums.RandomSums,\
                        frequentDirections.FrequentDirections]
    ns = [1000]
    ds = [100]
    ells = range(10,101,10)
    ks = [5]
    rounds = 1
        
    for (n, d, k) in product(ns, ds, ks):
        data_maker = SyntheticDataMaker()
        data_maker.initBeforeMake(d, k, signal_to_noise_ratio=10.0)
        A = data_maker.makeMatrix(n) # n * d matrix

        for (sketcherClass, ell, r) in product(sketcherClasses, ells, range(rounds)):
            if ell > d/2:
                continue
        
            sketcher = sketcherClass(d, ell)
            t_start = timer()
            for row in A:
                sketcher.append(row)
            t_end = timer()
                        
            totalSketchTime = t_end-t_start
            print(sketcher.class_name, totalSketchTime)
