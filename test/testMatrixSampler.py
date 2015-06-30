import sys
import unittest
sys.path.append('../sketch')
from utils.dataHandler import DataHandler
from matrixSampler import MatrixSampler as Sketcher

class testMatrixRandomSums(unittest.TestCase):

  def test_running(self):
    n = 100
    d = 20
    ell = 5
    data_handler = DataHandler()
    data_handler.initBeforeMake(d,signal_dimension=10,signal_to_noise_ratio=5,\
                                    signal_singular_value_decay_factor=1,signal_singular_value_decay_type='lin')
    
    sketcher = Sketcher(d,ell)

    for i in xrange(n):
        v = data_handler.makeRow()
        sketcher.append(v)
    
    sketch = sketcher.get()
    
    self.assertEqual(sketch.shape,(ell,d))


if __name__ == '__main__':
    unittest.main()