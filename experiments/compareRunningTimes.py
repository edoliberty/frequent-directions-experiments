import sys
from itertools import product
from time import time as timer
 
sys.path.append('../sketch') # needed for imports
from utils.syntheticDataMaker import SyntheticDataMaker
import bruteForce, frequentDirections, rowSampler, randomProjections, randomSums

if __name__ == '__main__':    
    scketcherClasses = [bruteForce.BruteForce, \
                        rowSampler.RowSampler,\
                        randomProjections.RandomProjections,\
                        randomSums.RandomSums,\
                        frequentDirections.FrequentDirections]
    ns = [1000]
    ds = [100]
    ells = range(10,101,10)
    ks = [5]
    rounds = 1
        
    for (sketcherClass, n, d, ell, k, r) in product(scketcherClasses, ns, ds, ells, ks, range(rounds)):
        if ell > d/2:
            continue
        
        data_maker = SyntheticDataMaker()
        data_maker.initBeforeMake(d, k, signal_to_noise_ratio=10.0)
            
        sketcher = sketcherClass(d,ell)
        
        t_start = timer()
        for i in xrange(n):
            vector = data_maker.makeRow()
            sketcher.append(vector)
        t_end = timer()
                        
        totalSketchTime = t_end-t_start
        print sketcher.class_name, totalSketchTime 
